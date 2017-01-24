        program free_ssbp
c
c calculate path using simulation of a chain.
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'

c
c local space
c       Be careful with space allocation which is
c       largely based on "maxpt" - maximum number of particles
c       which is in LENGTH.BLOCK. "maxpt" is usually set to a large
c       number (15,000). For a long chain that will be a LOT of space.
c       You may wish to modify
c       LENGTH.BLOCK for your own purposes. Decision on proper
c       lengths of vectors (for your molecule) should not be
c       difficult, since the connectivity file must be avaliable.
c
c **** ALLOCATE SPACE
c define double precision vectors of the length maxpt*3
c temp  = desired temperature
c dt    = time step
c ***special parameters for free energy calculations along RC ***
c urcrd - a unit number of a file with path binary coordinates
c numpth - the number of path structure of interest
c       *** NOTE THAT IN urcrd NPATH+1 POINTS ARE REQUIRED (INCLUDING
c               PRODUCTS) ***
c nstr  = number of structures along Reaction Coordinate
c nstep = number of sampling  points for a given React.Coord.
c debug = logical variable: true= a lot of information printed out
c rvrs   - logical variable if .true. start the reaction from prod.
c x     = coordinate vector
c y     = coordinate vector
c z     = coordinate vector
c vx   = velocity vector
c vy   = velocity vector
c vz   = velocity vector
c dx    = forces
c dy    = forces
c dz    = forces
c dxold = old forces (previous step)
c dyold = old forces (previous step)
c dzold = old forces (previous step)
c dmass = masses
c divms = 1 over the mass
c fact11 = a constant vector used in the Verlet integration
c fact22 =                    "
c grdp     - slopes of the curvlinear coordinate describing
c               the reaction coordinate at position numpth.
c               Estimated here by finite difference.
c pstep    - a step vector translation from point numpth to
c               numpth+1
c 
c define vectors for constraints
c  grdcm[x-z] = derivative of the center of mass with respect to [x-z]
c  grdl[x-z]  = derivatives of infitesimal rotations with respect to [x-z]
c pointr ipick - selection of particle
c ntest - period to test constraints
c newv  - period for selecting new velocities
c irand- integer seed for random number generator of velocities
c nwcrd - period for writing coordinates in a binary form
c uwcrd  -  unit for writing coordinates
c stdo - where some information on the path is printed out

Cdeb-----------
C pmfyes - logical if true then read a series of structures from
C          path file to find free energy differences and hence pmf
C lmbyes - logical if true then do lambda perturbation for the parameters
C          of the selected atoms 
C nlmbda - No of delta_lambda steps to be used to go from initial to
C          the final state.
C npert  -  Particle number of the ion to be perturbed
c
        double precision dpold(3*maxpt)
        double precision dmass(maxpt),divms(maxpt)
        double precision fact11(maxpt),fact22(maxpt)
        double precision grdcmx(3*maxpt),grdcmy(3*maxpt)
        double precision grdcmz(3*maxpt)
        double precision grdlx(3*maxpt),grdly(3*maxpt)
        double precision grdlz(3*maxpt)
        double precision grdp(3*maxpt),pstep(3*maxpt)

        integer pointr(maxpt),ipick(maxpt)

        double precision hot_umbr,beta,partit,exp2,qrot
        double precision temp,tmpr,dt,dt2,hlfdt,tfac

        double precision cortmp(3,maxpt)

        character*4 crdstyl
        character*5 name
        integer namel

        integer ucon,urcrd,uwcrd,uslop,ufrc
        integer uene
        integer begin,fin
        integer nstr,numpth,nstep,nsvel,newv,npri,ntest,nwcrd
        integer level,irand,npick,neqstep,nlist

        integer geti,of
        double precision getd
        logical find,fopen,rvrs

        double precision scalar(7),sigma(7)
        integer i,j,npick3,npt3,igrid
        logical massw

        logical shakm,freeze,toiyes,nori
        integer ndegf,ii

        logical pmfyes,lmbyes
        integer nlmbda  
 
        integer lc 
        data tfac/4.2532764d-2/
        data ucon,urcrd,uwcrd,uslop,ufrc,uene/6*99/
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

        debug  = .false.
        rvrs   = .false.
        nocut  = .false.
        name   = 'freee'
        namel  = 5

        lc     = 5

        shakm  = .false.
        tolcons= 1.d-5

        lcent  = .false.
        freeze = .false. 
        toiyes = .false.
        nori   = .false.
        qssbp  = .false.
   
        pmfyes = .false.
        lmbyes = .false.

        nlmbda = 10 
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c
c free energy default parameters
c
        igrid   = lgrid
        crdstyl = 'CHAR'
        nstep   = 1
        neqstep = 0
        nsvel   = 100
        newv    = 1000
        npri    = 1
        ntest   = 500
        nwcrd   = 500
        dt      = 0.001d0
        irand   = 1
        temp    = 300.d0
        begin   = 1
        fin     = 1
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
           if (find('conn')) then
               ucon  = of()
               call rconn(ucon)
           endif

           if (find('rcrd')) urcrd = of() 
           if (find('wcrd')) uwcrd = of()
           if (find('wslo')) uslop = of()
           if (find('wfrc')) ufrc  = of()
           if (find('wene')) uene  = of() 
        end if
        if (find('rvrs')) rvrs=.true.
        temp    = getd('temp',temp)
        igrid   = geti('grid',igrid)
        nstr    = geti('#str',nstr)
        nstep   = geti('#ste',nstep)
        neqstep = geti('#equ',neqstep)
        nsvel   = geti('#sve',nsvel)
        newv    = geti('newv',newv)
        npri    = geti('#pri',npri)
        ntest   = geti('#tes',ntest)
        nwcrd   = geti('#wcr',nwcrd)
Cdeb    rvrs    =(geti('rvrs',-1).gt.0)
        dt      = (getd('step',(dt)))
        irand   = geti('rand',irand)
        begin   = geti('bgin',begin)
        fin     = geti('fina',fin)
c Energy parameters
        nlist   = geti('list',nlist)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        eps     = (getd('epsi',(eps)))

        if (find('nocut')) nocut = .true.

        if (find('cdie')) ctrue = .true.
        if (find('rdie')) ctrue = .false.
Cdeb------------
        if (find('nori')) nori = .true.

        if (find('nobo')) ebyes  = .false.
        if (find('noan')) ethyes = .false.
        if (find('noto')) etoyes = .false.
        if (find('noim')) eimyes = .false.
        if (find('novd')) evdyes = .false.
        if (find('noel')) eelyes = .false.


        if (find('mshk')) then
          shakm  = .true.
          tolcons=getd('mtol',1.d-5)
          write(stdo,66)tolcons
66         format (1x,'matrix shake on TIP3/SPCE, tolcons=',d7.1)
          call mshakinit(debug)
          if (nshakm.le.0) then
                write (stdo,*) 'nothing to shake!'
                shakm=.false.
          endif
        endif
c
         if (find('symm')) then
          esymyes = .true.
          a = getd('xtra',0.0d0)
          b = getd('ytra',0.0d0)
          c = getd('ztra',0.0d0)
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
                  write(stdo,104)
104               format(1x,' Metal walls perpendicular',
     1          ' to the Y axis will be set!')
                  write(stdo,105)a0_metal
105               format(1x,' Metal repulsive Wall, A0 = ',
     1                  E12.6)
                  write(stdo,106)b_wall/2,v_elec
106               format(1x,' Walls are located at +/- ',f9.3,
     1                  1x,' The Electrode potential is ',f10.4)
       end if

c
         if (find('cent')) then
           lcent = .true.
           kcenter=getd('kcnt',10.d0)
           xeq = getd('xeqm',0.d0)
           yeq = getd('yeqm',0.d0)
           zeq = getd('zeqm',0.d0)
           call pick(ipick,i)
           icenter=0
           do 35 i=1,npt
             if (ipick(i).ne.0) then
               icenter = icenter+1
               center(icenter)=i
             endif
35         continue
         endif
c
         if (find('nfrz')) then
          freeze = .true.
          call pick(ipick,i)
          inofrz = 0
          do 44 i=1,npt
           if (ipick(i).ne.0) then
            inofrz = inofrz + 1
            nofreez(inofrz) = i
            zerofrz(i)= 1
           else
            zerofrz(i) = 0
           end if
 44       continue
         end if
         if (find('toiy')) toiyes = .true.
         if (find('selc')) then
                call pick(ipick,npick)
                npick=0
                do 3 i=1,npt
                        npick=npick+ipick(i)
3               continue
                j=0
                do 4 i=1,npt
                        if (ipick(i).gt.0) then
                                j=j+1
                                pointr(j)=i
                        end if
4               continue
         endif
Cdeb-------------------------------------------------------
Cdeb+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Spherical Solvent Boundary Potential (SSBP)
C
C Authors:  Dmitrii Beglovd and Benoit Roux  (1994)
C           Department of Chemistry, University of Montreal
C

C Default is include all energy terms for the 
C    solvent boundary potential
C----------------------------------------------------------
       if (find('ssbp')) call init_ssbp(ipick)

Cdeb-------------------------------------------------------
        if (find('acti')) go to 2
        go to 1
2       continue
c
Cdeb-------------------------------------------------------
c check that required files were opened
c
    
        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (.not.fopen(urcrd)) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        end if
Cdeb-------------------------------------------------------
c initialze no-freeze vector
c set the pointer to the selected particles. pointr(i) is 
c the position of selected atom number i in the normal all 
c atom array
        if (.not.freeze) then
        inofrz = npt
        do 21 i=1,inofrz
                nofreez(i) = i
                zerofrz(i) = 1
21      continue
        endif
Cdeb-------------------------------------------------------
c rmax is maintained here for old input (with a single
c cutoff to work)
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if
                
c
c internally the 1-4 scaling factors are defined as the inverse
        if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
        if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0


        cutvdw2  = cutvdw2*cutvdw2
        cutvbig2 = cutvbig2*cutvbig2
        cutele2  = cutele2*cutele2
        cutebig2 = cutebig2*cutebig2

        if (cutmono2.lt.0) then
                cutmono2 = cutebig2*1.44
        else
                cutmono2 = cutmono2*cutmono2
        end if

        if (npick.eq.0) then
         if (debug) write(stdo,*)' ipick = ',(ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
        end if
        write(stdo,100)temp,nstep,npri,igrid,stdo,uwcrd,urcrd
     1  ,nwcrd,nsvel
        write(stdo,101)dt,ntest,neqstep,newv,irand
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR FREE ENERGY SIMULATION:',//,
     1   1x,' temperature: ',f5.1,/,
     2   1x,' number of integration steps: ',i7,/,
     3   1x,' print each ',i7,' steps',/,
     4   1x,' number of path segments : ',i5,/,
     5   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     6   1x,' initial coordinates are read from unit : ',i5,/,
     7   1x,' coordinates are written at step interval : ',i5,/,
     8   1x,' velocity scaling at step interval of: ',i5,/)
101     format(1x,' step size: ',e10.3,/,
     1   1x,' period for checking constraints ',i5,/,
     1   1x,' number of thermalization steps  ',i5,/,
     2   1x,' select new velocities each ',i6,' steps',/,
     3   1x,' initial seed for random vel. selection ',i10,/)
        dt  = dt/tfac

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp



c
c
        npick3 = 3*npick
        npt3   = 3*npt
        if (debug) write(stdo,*) ' npick npt  ',npick,npt
c
c do mass weighting 
c
        massw=.true.
c
c calculate beta - inverse temperature in kcal/mol
c
        beta=1.d0/(temp*1.9878d-3)
c
c calculate the constants fact11 and fact22 dmass and 1/m (divms)
c
        dt2=dt*dt/2.d0
        hlfdt=dt/2.d0

        do 6 i=1,npt
                dmass(i) = ptms(i)
                divms(i) = 1.d0/dmass(i)
                invms(i) = divms(i)
                fact11(i) = dt2*divms(i)
                fact22(i) = hlfdt*divms(i)
6       continue

       IF (QSSBP) call set_ssbp()


c
c start loop on path coordinates
c
        do 2000 numpth=begin,fin
        partit=0.d0
        exp2=0.d0
        qrot=0.d0
c
c get initial path structure,
c path derivatives and the path step, vx vy vz & dxold dyold dzold
c are used as temporary vectors.
c
Cdeb-------------------------------------------

        call getpd_umbr(urcrd,grdp,pstep,cortmp
     1          ,velo,dpold,pointr,npick,numpth
     2          ,igrid,rvrs)
         if (numpth.eq.begin) then
           do 500 i = 1,npt
            coor(1,i) = cortmp(1,i)
            coor(2,i) = cortmp(2,i)
            coor(3,i) = cortmp(3,i)
 500       continue
Cdeb     if (esymyes) call squeeze(a,b,c)
Cdeb     call nbondm()
Cdeb     if (esymyes) call syminit()
         else
           do 510 i = 1,npick
            ii = pointr(i)
            coor(1,ii) = cortmp(1,ii)
            coor(2,ii) = cortmp(2,ii)
            coor(3,ii) = cortmp(3,ii)
 510       continue
c       
         endif 
Cdeb---------------------------------------------
c get center of mass constraints and orthonormalize all constraints
c dpold velo are used here as temporary vectors.
c
        if (.not.qssbp) then
Cdeb-------------Peptide Frozen----------
        call comc_umbr(coor,grdp,divms,npt,pointr,npick,grdcmx,
     1          grdcmy,grdcmz,grdlx,grdly,grdlz,dpold(1),
     2          dpold(1+npt),velo,sigma,debug)
Cdeb    if (debug) then
Cdeb            write(stdo,*)' grdp '
Cdeb            write(stdo,1000)(grdp(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdcmx '
Cdeb            write(stdo,1000)(grdcmx(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdcmy '
Cdeb            write(stdo,1000)(grdcmy(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdcmz '
Cdeb            write(stdo,1000)(grdcmz(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdlx '
Cdeb            write(stdo,1000)(grdlx(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdly '
Cdeb            write(stdo,1000)(grdly(i),i=1,npick3)
Cdeb            write(stdo,*) ' grdlz '
Cdeb            write(stdo,1000)(grdlz(i),i=1,npick3)
Cdeb    end if
1000            format(1x,8(f9.6))
1001            format(20X,3(f15.6))
        end if
c
c get velocities for initial condition
c
        call velinit(temp,1,0)
        if (debug) then
                write(*,*)' after first assignment irand = ',irand
                write(stdo,*)'vx = '
                write(stdo,1000)(velo(1,i),i=1,npt)
                write(stdo,*)'vy = '
                write(stdo,1000)(velo(2,i),i=1,npt)
                write(stdo,*)'vz = '
                write(stdo,1000)(velo(3,i),i=1,npt)
        end if
Cdeb-----------------------------------------------
        if (shakm) then
        ndegf=3*inofrz-nshakm
        else
        ndegf=3*inofrz
        endif
        tmpr=hot_umbr(velo,dmass,npt,ndegf,massw)
        write(stdo,25)tmpr
25      format(1x,'current temperature is ' , f9.3)
c
c       And now we are ready to call the FREE1 dynamics routine
c
        call dyna_free(temp,tmpr,dt,dt2,beta,partit,exp2,qrot,
     1          igrid,numpth,npri,nstep,nsvel,pointr,
     2          npick,dmass,divms,grdp,pstep,grdcmx,grdcmy,grdcmz,
     3          grdlx,grdly,grdlz,
     4          dpold,fact11,fact22,stdo,uwcrd,
     5          nwcrd,ntest,neqstep,newv,irand,debug,scalar,sigma,
     6          nlist,massw,ufrc,shakm,tolcons,toiyes,uene)
c
c final printing the results of TP
c
        write(stdo,*)' partition function ratio ',numpth,' is ',partit
        write(stdo,*)' free energy difference ',-1.d0/beta*dlog(partit)
        write(stdo,*)' partition function variance ',exp2
        if (exp2.le.0.d0) then
         write(stdo,*)' something fishy in the variance exp2=',exp2
        end if
2000    continue


        stop
        end
