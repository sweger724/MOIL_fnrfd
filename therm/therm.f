        program therm_cycle
        implicit none
c
c calculate free energy difference by thermodymamic
c perturbation  method
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONNECT2.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'


c
c local space
c
c **** ALLOCATE SPACE
c define double precision vectors of the chain length
c coor  = coordinate vector
c velo  = velocity vector
c dpot  = forces
c dpold = old forces (previous step)
c dmass = masses
c divms = 1 over the mass
c fact1 = a constant vector used in the Verelt integration
c fact2 =                    
c define vectors for constraints
c  grdcm[x-z] = derivative of the center of mass with respect to [x-z]
c  grdl[x-z]  = derivatives of infitesimal rotations with respect to [x-z]
c pointr pick1 - selection of particle
c shakl   - If true, shake is activated for particles m <= 1.1
c shakb   - If true all bonds are shaked
c
c conversion factor from ps to KAM

        double precision dpold(3,maxpt)
        double precision dmass(maxpt),divms(maxpt)
        double precision fact1(maxpt),fact2(maxpt)
        double precision grdlx(3*maxpt),grdly(3*maxpt)
        double precision grdlz(3*maxpt)

        double precision hot_umbr,beta,partit,exp2,lambda,delta
        double precision temp,tmpr,dt,dt2,hlfdt,tfac


        logical shakl,shakb,shakm

        integer i


        character*4 crdstyl
        character*5 name
        integer namel

        integer uconr,uconp,urcrd,uwcrd,uwvel
        integer nstep,nsvel,newv,npri,ntest,nwcrd,irand,neqstep
        integer level,nlist,neqsvel

        integer geti,of
        logical find,fopen

        double precision winda,windb
        double precision getd
        integer npt3,ndegf
        logical massw

        data tfac/4.2532764d-2/
        data uconr,uconp,urcrd,uwcrd,uwvel/5*99/

c
c General initialization
c
        stdi   = 5
        stdo   = 6
        stderr = 0
        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        nbulk  = 0

c SHAKE ADDITION

        nshak  = 0
c
c inittialization for exact number of COM of new CO model
c

        debug  = .false.
        nocut  = .false.
        shift  = .false.


        shakl  = .false.
        shakb  = .false.
        shakm  = .false.

        name   = 'therm'
        namel  = 5
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c
c chain default parameters
c
        winda  = 0.d0
        windb  = 1.d0
        delta  = 5.d-2
        crdstyl = 'CHAR'
        nstep   = 1
        neqstep = 0
        neqsvel = 10
        nsvel   = 100
        newv    = 1000
        npri    = 1
        ntest   = 500
        nwcrd   = 500
        dt      = 0.001d0
        irand   = -23518284
        temp    = 300.d0


        epshak   = 1.d-6
        epshakv  = 1.d-4
        itershak = 400
c
c energy default parameters
c
        call init_ef()

        nlist  = 1

c
c start interpret the line from stdi
c
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
         if (find('conr')) uconr = of()
         if (find('conp')) uconp = of()
         if (find('rcrd')) urcrd = of()
         if (find('wcrd')) uwcrd = of()
         if (find('wvel')) uwvel = of()
        end if
        temp    = (getd('temp',(temp)))
        nstep   = geti('#ste',nstep)
        neqstep = geti('#equ',neqstep)
        neqsvel = geti('#eqv',neqsvel)
        nsvel   = geti('#sve',nsvel)
        newv    = geti('newv',newv)
        npri    = geti('#pri',npri)
        ntest   = geti('#tes',ntest)
        nwcrd   = geti('#wcr',nwcrd)
        dt      = (getd('step',(dt)))
        irand   = geti('rand',irand)
c
c Energy parameters
c
        nlist   = geti('list',nlist)
          cutvdw2  = (getd('rvmx',(cutvdw2)))
          cutvbig2 = (getd('rvbg',cutvbig2))
          cutele2  = (getd('relx',(cutele2)))
          cutebig2 = (getd('rebg',cutebig2))
          cutmono2 = (getd('cutm',cutmono2))
          rmax     = getd('rmax',rmax)

        eps     = (getd('epsi',(eps)))
        winda   = (getd('firs',(winda)))
        windb   = (getd('fina',(windb)))
        delta   = (getd('sted',(delta)))
        if (find('nocut')) nocut = .true.
        if (shift  .or.  find('shif')) shift   = .true.
        if (find('cdie')) ctrue = .true.
        if (find('rdie')) ctrue = .false.
        if (find('shkb') .and. (.not.shakb)) shakb=.true.

        if (find('cpth')) crdstyl = 'PATH'
        if (find('cchr')) crdstyl = 'CHAR'
        if (find('cini')) crdstyl = 'INIT'
        if (find('cint')) crdstyl = 'INTR'

        if (find('nobo')) ebyes  = .false.
        if (find('noan')) ethyes = .false.
        if (find('noto')) etoyes = .false.
        if (find('noim')) eimyes = .false.
        if (find('novd')) evdyes = .false.
        if (find('noel')) eelyes = .false.
        if (find('symm')) then
         esymyes = .true.
         a = getd('xtra',0.d0)
         b = getd('ytra',0.d0)
         c = getd('ztra',0.d0)
        end if

        if (find('hvdw')) hvdw0 = .false.

c To be added after call rline of the relevant procedure
c
c jmjm
c
c switch for Ewald summation of long range interactions
         if (find('ewald')) then
            ewaldyes = .true.
            dtol = getd('dtol',0.0d0)
            nfft1 = geti('grdx',0)
            nfft2 = geti('grdy',0)
            nfft3 = geti('grdz',0)
            sgridx = getd('sgdx',1.0d0)
            sgridy = getd('sgdy',1.0d0)
            sgridz = getd('sgdz',1.0d0)
            intrpord = geti('iord',4)
            write (stdo,*)
     &  'PME account for long range inter. will be performed'

c if one wants to use PME for vacuum calculations
c one needs a virtual box which is sufficiently large to effectively
c separate image boxes
c in such a case stop should be commented and proper sizes
c of the virtual box  xtra,ytra,ztra should be given

            if (.not.esymyes) then
               write (stdo,*)
     &  'There is no true periodic bound. cond. - symm is missing'
               stop
c              a = getd('xtra',50.0d0)
c              b = getd('ytra',50.0d0)
c              c = getd('ztra',50.0d0)
            end if
         end if
c jmjmend

c
c virtual particles: massless, point charges displaced from the motion
c                    centers (real particles) e.g. TIP4P, 3-site CO model
c
         if (find('vprt')) then
           vp_flag = .true.
           com_flag = .true.
           gcnt_flag = .false.
           if (find('gcnt')) then
              com_flag = .false.
              gcnt_flag = .true.
           end if
         end if
c end of virtual partcles - one may actually think of moving this to con
c
c carlos, add amide plane constraints
         if (find('amid')) call amid()

       if (find('metl')) then
              metalyes   = .true.
              a0_metal   = getd('amtl',50.d0)
              alfa_metal = getd('alfa',1.d0)
                  b_wall     = getd('bwal',b)
                  if (b_wall.gt.b) then
                   call alert('dyna',4,'b_wall must be < b',18,1)
                  end if
                  v_elec     = getd('v_el',0.d0)
c
c compute pre-exponential  factor A0_metal such that
c A0_metal*exp(-alfa_metal*b/2)=a0_metal
c
c              a0_metal = a0_metal*dexp(-0.5d0*alfa_metal*b_wall)
                  write(stdo,*)
                  write(stdo,114)
114               format(1x,' Metal walls perpendicular',
     1          ' to the Y axis will be set!')
                  write(stdo,115)a0_metal
115               format(1x,' Metal repulsive Wall, A0 = ',
     1                  E12.6)
                  write(stdo,116)b_wall/2,v_elec
116               format(1x,' Walls are located at +/- ',f9.3,
     1                  1x,' The Electrode potential is ',f10.4)
       end if

        if (find('acti')) go to 2
        go to 1
2       continue

c set the freezing pointer (zerofreez) to total no freez. At present
c no other options are supported. It is need in shakinit. No shake
c on frozen bonds will be pursued.
c
        do 3 i=1,npt
                zerofrz(i) = 1
3       continue

c rmax is maintained here for old input (with a single
c cutoff to work)
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if


         if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
         if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

c
c check that required files were opened
c

        if (.not.fopen(uconr)) then
         level = 1
         call alert(name,namel,'uconr not opened',15,level)
        else if (.not.fopen(uconp)) then
         level = 1
         call alert(name,namel,'uconp not opened',15,level)
        else if (.not.fopen(urcrd)) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        else if (.not.fopen(uwvel)) then
         level = 1
         call alert(name,namel,'uwvel not opened',16,level)
        end if
c

         cutvdw2  = cutvdw2*cutvdw2
         cutvbig2 = cutvbig2*cutvbig2
         cutele2  = cutele2*cutele2
         cutebig2 = cutebig2*cutebig2

         if (cutmono2.lt.0) then
                cutmono2 = cutebig2*1.44
         else
                cutmono2 = cutmono2*cutmono2
         end if

c *************************************************************
c get connectivity file for reactants (lambda=0)
c *************************************************************
        call rconn(uconr)
c initialize no freeze vector
        inofrz = npt
        do 200 i=1,inofrz
         zerofrz(i) = 1
200     continue
        do 201 i=1,npt
        ptmsA(i)=ptms(i)
        ptchgA(i)=ptchg(i)
        epsgm6A(i)=epsgm6(i)
        epsgm12A(i)=epsgm12(i)
        ptweiA(i)=ptwei(i)
201     continue
        do 202 i=1,nb
        kbondA(i)=kbond(i)
        reqA(i)=req(i)
202     continue
        do 203 i=1,nangl
        kanglA(i)=kangl(i)
        angleqA(i)=angleq(i)
203     continue
        do 204 i=1,ntors
        periodA(i)=period(i)
        ktors1A(i)=ktors1(i)
        ktors2A(i)=ktors2(i)
        ktors3A(i)=ktors3(i)
        phase1A(i)=phase1(i)
        phase2A(i)=phase2(i)
        phase3A(i)=phase3(i)
204     continue
        do 205 i=1,nimp
        kimpA(i)=kimp(i)
        impeqA(i)=impeq(i)
205     continue


c *****************************************************************
c get connectivity file for products(lambda=1)
c *****************************************************************
        call rconn(uconp)
        do 206 i=1,npt
        ptmsB(i)=ptms(i)
        ptchgB(i)=ptchg(i)
        epsgm6B(i)=epsgm6(i)
        epsgm12B(i)=epsgm12(i)
        ptweiB(i)=ptwei(i)
206     continue
        do 207 i=1,nb
        kbondB(i)=kbond(i)
        reqB(i)=req(i)
207     continue
        do 208 i=1,nangl
        kanglB(i)=kangl(i)
        angleqB(i)=angleq(i)
208     continue
        do 209 i=1,ntors
        periodB(i)=period(i)
        ktors1B(i)=ktors1(i)
        ktors2B(i)=ktors2(i)
        ktors3B(i)=ktors3(i)
        phase1B(i)=phase1(i)
        phase2B(i)=phase2(i)
        phase3B(i)=phase3(i)
209     continue
        do 210 i=1,nimp
        kimpB(i)=kimp(i)
        impeqB(i)=impeq(i)
210     continue

        if (urcrd .eq. 99 ) then
                write(*,*)' Sorry, initial coordinate file not provided'
                write(*,*)' missing  - urcrd number '
                write(*,*)' Quitting...., better luck next time'
                level = 1
                call alert(name,namel,'Problem in file',15,level)
        end if
        write(stdo,100)temp,nstep,npri,winda,windb,delta,
     1  stdo,uwcrd,urcrd,crdstyl,nwcrd,nsvel
        write(stdo,101)dt,ntest,neqstep,newv,irand
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR THERMODYNAMIC PERTURBATION :',//,
     1   1x,' temperature: ',f5.1,/,
     2   1x,' number of integration steps: ',i7,/,
     3   1x,' print each ',i7,' steps',/,
     4   1x,' initial value of  LAMBDA : ',f10.5,/,
     5   1x,' final value of  LAMBDA : ',f10.5,/,
     6   1x,' step  of  LAMBDA : ',f10.5,/,
     7   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     8   1x,' initial coordinates are read from unit : ',i5,/,
     9   1x,' Coordinate are read in ',a4,' style ',/,
     1   1x,' coordinates are written at step interval : ',i5,/,
     2   1x,' velocity scaling at step interval of: ',i5)
101     format(1x,' step size: ',e10.3,/,
     1   1x,' period for checking constraints ',i5,/,
     2   1x,' number of thermalization steps  ',i5,/,
     3   1x,' select new velocities each ',i6,' steps',/,
     4   1x,' initial seed for random vel. selection ',i10)
        dt  = dt/tfac


        npt3   = 3*npt

        if ( shakb .or. shakl ) then
          nshak=0
          call shakinit(shakl,shakb,shakm,epshak)
          end if
                ndegf = npt3 - nshak
                write(stdo,106) ndegf
106             format(1x,'Number of degf (exc. shake) ',i7)



c
c do mass weighting
c
        massw=.true.
c
c calculate the constants beta,fact1 and fact2 dmass and 1/m (divms)
c
        beta=1.d0/(temp*1.9878d-3)
        dt2=dt*dt/2.d0
        hlfdt=dt/2.d0

c
c get coordinates
c
        call getcrd(urcrd, 'CHARM')

        do 99 i=1,npt
        coor2(1,i)=coor(1,i)
        coor2(2,i)=coor(2,i)
        coor2(3,i)=coor(3,i)
99      continue

c
c start loop on lambda coordinates
c

        do 2000 lambda=winda,windb,delta
        partit=0.d0
        exp2=0.d0


          do 355 i=1,npt
          ptms(i)=(1.d0-lambda)*ptmsA(i)+ lambda*ptmsB(i)
          dmass(i) =(1.d0-lambda)*ptmsA(i)+lambda*ptmsB(i)
                divms(i) = 1.d0/dmass(i)
                fact1(i) = dt2*divms(i)
                fact2(i) = hlfdt*divms(i)
          ptwei(i)=(1.d0-lambda)*ptweiA(i)+ lambda*ptweiB(i)
          ptchg(i)=(1.d0-lambda)*ptchgA(i)+ lambda*ptchgB(i)
          epsgm6(i)=(1.d0-lambda)*epsgm6A(i)+
     1              lambda*epsgm6B(i)
          epsgm12(i)=(1.d0-lambda)*epsgm12A(i)+
     1              lambda*epsgm12B(i)
355       continue
          do 356 i=1,nb
          kbond(i)=(1.d0-lambda)*kbondA(i)+ lambda*kbondB(i)
          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)
356       continue
          do 357 i=1,nangl
          kangl(i)=(1.d0-lambda)*kanglA(i)+ lambda*kanglB(i)
          angleq(i)=(1.d0-lambda)*angleqA(i)+ lambda*angleqB(i)
357       continue
          do 358 i=1,ntors
          period(i)=(1.d0-lambda)*periodA(i)+lambda*periodB(i)
          ktors1(i)=(1.d0-lambda)*ktors1A(i)+ lambda*ktors1B(i)
          ktors2(i)=(1.d0-lambda)*ktors2A(i)+ lambda*ktors2B(i)
          ktors3(i)=(1.d0-lambda)*ktors3A(i)+ lambda*ktors3B(i)
          phase1(i)=(1.d0-lambda)*phase1A(i) + lambda*phase1B(i)
          phase2(i)=(1.d0-lambda)*phase2A(i) + lambda*phase2B(i)
          phase3(i)=(1.d0-lambda)*phase3A(i) + lambda*phase3B(i)
358       continue
          do 359 i=1,nimp
          kimp(i)=(1.d0-lambda)*kimpA(i)+ lambda*kimpB(i)
          impeq(i)=(1.d0-lambda)*impeqA(i)+ lambda*impeqB(i)
359       continue


c ************************************************
        if ( shakb .or. shakl ) then
          nshak=0
          call shakinit(shakl,shakb,shakm,epshak)
          end if

c ************************************************

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp








c
c assign velocities (initial conditions)
c
         call velinit(temp,1,0)


c factor out center of mass velocity and
c factor out the rigid body rotation by orienting the structure against
c the original configuration and also by calculating the gradient of
c the constraints on rotation and projecting out the components
c of the velocities along these directions
c

         if (.not.esymyes) then
         call prbm(velo,coor,dmass,grdlx,grdly,grdlz,npt)
         end if

              if (debug) then
                write(stdo,*)' after first assignment irand = ',irand
                write(stdo,*)'vx = '
                write(stdo,1000)(velo(1,i),i=1,npt)
                write(stdo,*)'vy = '
                write(stdo,1000)(velo(2,i),i=1,npt)
                write(stdo,*)'vz = '
                write(stdo,1000)(velo(3,i),i=1,npt)
                end if
1000            format(1x,7(f10.5,1x))
c
c get temperature from hot_umbr

        tmpr = hot_umbr(velo,dmass,npt,ndegf,massw)
        write(stdo,102)tmpr
102     format(1x,'current temperature is ',f12.3)


c
c       And now we are ready to call the FREE1 dynamics routine
c
       call dyna_therm(temp,tmpr,dt,dt2,hlfdt,
     1  tfac,beta,partit,exp2,
     1  winda,windb,delta,lambda,npri,nstep,neqsvel,nsvel,
     2  dmass,divms,
     3  grdlx,grdly,grdlz,dpold,
     4  fact1,fact2,stdo,uwcrd,
     5  nwcrd,ntest,neqstep,newv,irand,debug,
     6  shakb,shakl,shakm,
     6  nlist,ndegf,massw)
c
c final printing the results of TP
c
        write(stdo,*)' free energy difference ', partit
        write(stdo,*)' partition function variance ',exp2
        if (exp2.le.0.d0) then
         write(stdo,*)' something fishy in the variance exp2=',exp2
        end if
2000    continue

        stop
        end
