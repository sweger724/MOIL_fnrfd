       program umbrella
c
c calculate free energies using the umbrella sampling method
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
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'


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
c upcrd - a unit number of a file with path binary coordinates
c numpth - the number of path structure of interest
c       *** NOTE THAT IN upcrd NPATH+1 POINTS ARE REQUIRED (INCLUDING
c               PRODUCTS) ***
c nstep = number of sampling  points for a given React.Coord.
c debug = logical variable: true= a lot of information printed out
c rvrs   - logical variable if .true. start the reaction from prod.
c coor(i,j)  = coordinate vector i=x,y,z ; j=1,..,npt
c velo(i,j)  = velocity vector i=x,y,z ; j=1,...,npt
c dpot(i,j)  = forces indices as above
c dpold = old forces (previous step)
c dmass = masses
c divms = 1 over the mass
c fact11 = a constant vector used in the Verlet integration
c fact22 =                    "
c grdp     - slopes of the curvlinear coordinate describing
c               the reaction coordinate at position numpth.
c               Estimated here by finite difference, NOT NORMALIZED.
c pstep    - a step vector translation from point numpth to
c               numpth+1
c 
c istart & ifinal - the first and the last structures of the path files
c               to be examined in the present run
c igrid - number of grid point along the reaction path.
c istep - a step along the rection coordinate structures (default is one)
c         useful when adding intermediate points
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
c
        double precision dpold(3,maxpt)
        double precision dmass(maxpt),divms(maxpt)
        double precision fact11(maxpt),fact22(maxpt)
        double precision grdcmx(3,maxpt),grdcmy(3,maxpt)
        double precision grdcmz(3,maxpt)
        double precision grdlx(3,maxpt),grdly(3,maxpt)
        double precision grdlz(3,maxpt)
        double precision grdp(3,maxpt),pstep(3,maxpt)

        integer pointr(maxpt),ipick(maxpt)

        double precision temp,tmpr,dt,dt2,hlfdt,tfac,xfc
        double precision tmp,rcms,rcnorm

        character*5 name
        integer namel

        integer istart,ifinal
        integer ucon,uwcrd,uwvel,upcrd,uqint
        integer numpth,nstep,nsvel,newv,npri,ntest,nwcrd
        integer level,irand,npick,neqstep,nlist

        integer geti,of
        double precision getd
        logical find,fopen,rvrs,effmass

        double precision scalar(7),sigma(7),qq2,q0
        double precision xfc2,eu
        integer i,j,npick3,npt3,igrid,istep
        logical massw

        data tfac/4.2532764d-2/
        data ucon,uwcrd,uwvel/3*99/
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
        effmass= .false.
        name   = 'umbre'
        namel  = 5
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c
c chain default parameters
c
        igrid   = lgrid
        istep   = 1
        istart  = 1
        ifinal  = -1
        nstep   = 1
        neqstep = 0
        nsvel   = 100
        newv    = 1000
        npri    = 1
        ntest   = 500
        nwcrd   = 500
        dt      = 0.001d0
        irand   = -23518284
        temp    = 300.d0
        xfc     = 40.d0
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
         if (find('conn')) ucon = of()
         if (find('rcrd')) upcrd = of()
         if (find('quni')) uqint = of()
         if (find('wcrd')) uwcrd = of()
         if (find('wvel')) uwvel = of()
        end if
        temp    = getd('temp',temp)
        igrid   = geti('grid',igrid)
        istep   = geti('istp',istep)
        istart  = geti('strt',istart)
        ifinal  = geti('finl',ifinal)
        nstep   = geti('#ste',nstep)
        neqstep = geti('#equ',neqstep)
        nsvel   = geti('#sve',nsvel)
        newv    = geti('newv',newv)
        npri    = geti('#pri',npri)
        ntest   = geti('#tes',ntest)
        nwcrd   = geti('#wcr',nwcrd)
        rvrs    =(geti('rvrs',-1).gt.0) .or. rvrs
        effmass = (find('effm').or.effmass)
        dt      = getd('step',dt)
        xfc     = getd('forc',xfc)
        irand   = geti('rand',irand)
c Energy parameters
        nlist   = geti('list',nlist)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        if (find('gbsa')) gbsabool = .true.
        gbsu     = geti('gbsu',gbsu)

        eps     = getd('epsi',eps)
        if (find('nocut')) nocut = .true.

        if (find('cdie')) ctrue = .true.
        if (find('shif')) shift = .true.
        if (find('rdie')) ctrue = .false. 

        if (find('nobo')) ebyes  = .false.
        if (find('noan')) ethyes = .false.
        if (find('noto')) etoyes = .false.
        if (find('noim')) eimyes = .false.
        if (find('novd')) evdyes = .false.
        if (find('noel')) eelyes = .false.

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
        if (find('acti')) go to 2
        go to 1
2       continue

        if (ifinal.lt.0) ifinal  = igrid

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

        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (.not.fopen(upcrd)) then
         level = 1
         call alert(name,namel,'upcrd not opened',16,level)
        else if (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        else if (.not.fopen(uwvel)) then
         level = 1
         call alert(name,namel,'uwvel not opened',16,level)
        end if

        cutvdw2  = cutvdw2*cutvdw2
        cutvbig2 = cutvbig2*cutvbig2
        cutele2  = cutele2*cutele2
        cutebig2 = cutebig2*cutebig2

         if (cutmono2.lt.0) then
                cutmono2 = cutebig2*1.44
         else
                cutmono2 = cutmono2*cutmono2
         end if

c get connectivity
        call rconn(ucon)

        if (gbsabool) call make_rborn
c initialze zero freeze vector
        inofrz = npt
        do 21 i = 1,inofrz
                zerofrz(i) = 1
21      continue

c 
c
        call rline(name,namel,stdi)
        call pick(ipick,npick)
c
c how many atoms where selected ??
c
        npick=0
        do 3 i=1,npt
                npick=npick+ipick(i)
3       continue
        if (npick.eq.0) then
         if (debug) write(stdo,*)' ipick = ',(ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
        end if
        if (upcrd .eq. 99 ) then
                write(*,*)' Sorry, initial coordinate file not provided'
                write(*,*)' missing  - upcrd number '
                write(*,*)' Quitting...., better luck next time'
                level = 1
                call alert(name,namel,'Problem in file',15,level)
        end if
        write(stdo,100)temp,nstep,npri,igrid,stdo,uwcrd,upcrd
     1  ,nwcrd,nsvel,xfc
        write(stdo,101)dt,ntest,neqstep,newv,irand
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR UMBRELLA SIMULATION:',//,
     1   1x,' temperature: ',f5.1,/,
     2   1x,' number of integration steps: ',i7,/,
     3   1x,' print each ',i7,' steps',/,
     4   1x,' number of path segments : ',i5,/,
     5   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     6   1x,' initial coordinates are read from unit : ',i5,/,
     8   1x,' coordinates are written at step interval : ',i5,/,
     9   1x,' velocity scaling at step interval of: ',i5,/,
     9   1x,' force constant for umbrella potential: ',f5.1,/)
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
c set the pointer to the selected particles. pointr(i) is the position
c of selected atom number i in the normal all atom array
c
        j=0
        do 4 i=1,npt
                if (ipick(i).gt.0) then
                        j=j+1
                        pointr(j)=i
                end if
4       continue
        npick3 = 3*npick
        npt3   = 3*npt
        if (debug) write(stdo,*) ' npick npt  ',npick,npt
c
c do mass weighting 
c
        massw=.true.
c
c calculate the constants fact11 and fact22 dmass and 1/m (divms)
c
        dt2=dt*dt/2.d0
        hlfdt=dt/2.d0

        do 6 i=1,npt
                dmass(i) = ptms(i)
                divms(i) = 1.d0/dmass(i)
                fact11(i) = dt2*divms(i)
                fact22(i) = hlfdt*divms(i)
6       continue

        q0=0.d0


        do 2000 numpth=istart,ifinal,istep
c
c get initial path structure,
c path derivatives and the path step, vx vy vz & dxold dyold dzold
c are used as temporary vectors.
c

        call getpd_umbr(upcrd,grdp,pstep,coor
     1          ,velo,dpold,pointr,npick,numpth
     2          ,igrid,rvrs)


        write(*,*)' effmass ',effmass
        if (effmass) then

                rcnorm = 0.d0
                rcms   = 0.d0
                do 11 j=1,npick
                 i   = pointr(j)
                 tmp = grdp(1,j)*grdp(1,j) + grdp(2,j)*grdp(2,j)
     1                  + grdp(3,j)*grdp(3,j)
                 rcnorm = rcnorm + tmp
                 rcms   = rcms + ptms(i)*tmp
11              continue
                rcms = rcms/rcnorm
                write(*,*)numpth,rcms
                go to 2000
        end if
c
c all is required is to calculate the effective mass along the 
c reaction path.
c
                

1000            format(1x,8(f9.6))
                if (debug) then 
                do 9995 j=1,npick
                write(stdo,*) ' grdpx ' , grdp(1,j)
                write(stdo,*) ' grdpy ' , grdp(2,j)
                write(stdo,*) ' grdpz ' , grdp(3,j)
9995            continue

                end if

                do 155 i=1,npt
                coor2(1,i)=coor(1,i)
                coor2(2,i)=coor(2,i)
                coor2(3,i)=coor(3,i)
155             continue

                qq2=0.d0
                do 156 i=1,npick
                qq2=qq2+pstep(1,i)*pstep(1,i)+pstep(2,i)*pstep(2,i)
     1                  +pstep(3,i)*pstep(3,i)
156             continue
                qq2=dsqrt(qq2)

                call nbondm()
                call eforce()

        eu = 0.d0
c the force constant should be normalized to the number of degf
        xfc2 = 2.d0*xfc
        do 311 i=1,npick
                j=pointr(i)
                eu = eu + ((coor(1,j)-coor2(1,j))*grdp(1,i) +
     1          (coor(2,j)-coor2(2,j))*grdp(2,i) +
     2          (coor(3,j)-coor2(3,j))*grdp(3,i)) 
311     continue
        do 312 i=1,npick
                j=pointr(i)
                dpot(1,j) = dpot(1,j) + xfc2 * eu *grdp(1,i)
                dpot(2,j) = dpot(2,j) + xfc2 * eu *grdp(2,i)
                dpot(3,j) = dpot(3,j) + xfc2 * eu *grdp(3,i)
312     continue
        eu = xfc*eu*eu
                write(*,*) ' eu is ' , eu
                e_total=e_total+eu
c
c get center of mass constraints and orthonormalize all constraints
c dxold dyold vx vy vz are used here as temporary vectors.
c
        call comc_umbr(coor2,grdp,divms,npt,pointr,npick,
     1          grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,dpold(1,1),
     1          dpold(1,npt/2),velo,sigma,debug)
        if (debug) then
                write(stdo,*)' grdp '
                write(stdo,1000)((grdp(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdcmx '
                write(stdo,1000)((grdcmx(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdcmy '
                write(stdo,1000)((grdcmy(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdcmz '
                write(stdo,1000)((grdcmz(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdlx '
                write(stdo,1000)((grdlx(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdly '
                write(stdo,1000)((grdly(j,i),j=1,3),i=1,npick)
                write(stdo,*) ' grdlz '
                write(stdo,1000)((grdlz(j,i),j=1,3),i=1,npick)
        end if
c
c get velocities for initial condition
c
c
c       And now we are ready to call the DYNA_UMBR dynamics routine
c
        call dyna_umbr(temp,tmpr,dt,tfac,xfc,q0,
     1          igrid,numpth,npri,nstep,nsvel,pointr,
     2          npick,dmass,divms,grdp,pstep,grdcmx,grdcmy,grdcmz,
     3          grdlx,grdly,grdlz,
     4          dpold,fact11,fact22,uqint,uwcrd,
     5          nwcrd,ntest,neqstep,newv,irand,debug,scalar,sigma,
     6          nlist,massw)

                q0=q0+qq2

2000    continue

        stop
        end
