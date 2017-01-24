       program chmin
       implicit none
c
c calculate path using conjugate gradient minimization
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
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/PDQ.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'
        include 'COMMON/RGYRCST.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'
c        include 'COMMON/PATH.BLOCK'

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
c define double precision vectors of the chain length 
c r     = coordinate vector
c vel   = velocity vector
c dv    = forces
c buffer= temporary space for minimizer

c define vectors for constraints
c  d0 - distances between intermediates (i,i+1 AND i,i+2)
c  e0 - vecotr storing the energies of individual "time" frame
c  e1 - work vector used in the polymer energy routine
c  grdcm[x-z] = derivative of the center of mass with respect to [x-z]
c  grdl[x-z]  = derivatives of infitesimal rotations with respect to [x-z]
c pointr pick1 - selection of particle
c
c total local space (in vectors required):
c               24*maxpt*(lgrid+6) +  4*lgrid - 3 (this includes vectors defined
c               in pwl_chn)
c
        double precision r(3*maxpt*lgrid)
        double precision dv(3*maxpt*lgrid)
        double precision e0(lgrid),e1(2*lgrid),d0(2*lgrid-3)
        double precision grdcmx(3*maxpt),grdcmy(3*maxpt)
        double precision grdcmz(3*maxpt)
        double precision grdlx(3*maxpt),grdly(3*maxpt)
        double precision grdlz(3*maxpt)

        integer pointr(maxpt),ipick(maxpt),irand

        double precision gamma,rho,lambda,rms
        double precision temper

        character*4 crdstyl
        character*3 name
        integer namel

        integer ucon,urcrd,uwcrd
        integer ucon1,ucon2
        integer nstep,npri,ntest,nwcrd,npick
        integer imx1,ix1
        integer level,nlist
        integer sgbboolint

        integer geti,of
        double precision getd
        logical find,fopen,fixend,lap
        logical anneal
        logical select

        double precision scalar(6),sigma(6),sigmav(6)
        double precision tolg,estred
        integer i,j,npick3,npt3,igrid,middle
        integer udata

        data ucon,urcrd,uwcrd/3*99/
c
c General initialization
c
        stdi   = 5
        stdo   = 6
        udata  = stdo
        stderr = 0
        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        nbulk  = 0

        select = .false.
        debug  = .false.
        name   = 'sdp'
        namel  = 3
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c
c chain default parameters
c
        anneal  = .false.
        igrid   = lgrid
        crdstyl = 'CHAR'
        nstep   = 1
        npri    = 1
        ntest   = 500
        nwcrd   = 500
        temper  = 300.d0
        fixend  = .true.
c gamma - the "spring" force constant in the chain describing
c               attraction between nearest neighbors
c               units:  kcal/mol angstrom-2
        sgbboolint = 0
        gamma = 100.d0
c rho   - the repulsion parameter between next nearest neighbors
c               units: kcal/mol
        rho   = 100.d0
c lambda- range parameter in the repulsion exponent - exp(-lambda d^2/<d>^2)
c
        lambda= 2.d0
c tolg - allowed gradient tolerance
c
        tolg  = 1.d-3
        estred = 0.01d0
        sgba = 0
c do an overlap of the structures with respect to each other?
c
        lap = .false.
c
c energy default parameters
c
        call init_ef()
        nocut  = .true.

        nlist  = 1

C default ELastic network model parameters
       Hamilt= 1.d-10

c
c start interpret the line from stdi
c
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
         if (find('conn')) ucon = of()
         if (find('rcrd')) urcrd = of()
         if (find('wcrd')) uwcrd = of()
         if (find('con1')) ucon1 = of()
         if (find('con2')) ucon2 = of()

C   read coarse grained model parameters (associated files)
         call Read_CG_input()
         
        end if

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

        if (find('hvdw')) hvdw0  = .false.
        if (find('anne')) anneal = .true.
        igrid   = geti('grid',igrid)
        if (igrid.gt.lgrid) then
         call alert(name,namel,'Lgrid exceeded',14,1)
        end if
        if (find ('sgbb')) then
           sgbbool =.true.
        end if
        if (find('gbsa')) gbsabool = .true.
        gbsu = geti('gbsu',gbsu)
        nstep   = geti('#ste',nstep)
        npri    = geti('#pri',npri)
        ntest   = geti('#tes',ntest)
        nwcrd   = geti('#wcr',nwcrd)
        temper  = getd('tmpr',temper)
        gamma   = (getd('gama',(gamma)))
        rho     = (getd('repl',(rho)))
        lambda  = (getd('lmbd',(lambda)))
        tolg    = (getd('tolg',(tolg)))
        estred  = (getd('dfpr',(estred)))
        if (find('ovlp')) lap = .true.
        if (find('sele')) select = .true.
c Energy parameters
        if (find('hvdw')) hvdw0 = .false.
        nlist   = geti('list',nlist)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        eps     = (getd('epsi',(eps)))

        if (find('cdie')) ctrue = .true.
        if (find('rdie')) ctrue = .false.


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
        irand   = geti('rand',irand)

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

c rmax is maintained here for old input (with a single
c cutoff to work)
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if

c
c check that required files were opened
c

        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (.not.(fopen(urcrd) .or. urcrd.eq.5)) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        end if

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

c get connectivity
        call rconn(ucon)

c initialize masses of CA atoms
        if (eCGyes) call initCGmass()
        if (eenmyes) then
          do  i=1,npt
            ptms(i)=13
          end do
        end if
        
c initialze no freez vecotr
        inofrz = npt
        do 21 i=1,npt
                zerofrz(i) = 1
21      continue

c Pick atoms (note that chain/chmin is currently NOT working
c with pick, hoping to fix this soon)
c
        if (select) then
                call rline(name,namel,stdi)
                call pick(ipick,npick)
                write(*,*)' after pick '
        j=0
        do 4 i=1,npt
                if (ipick(i).ne.0) then
                 j=j+1
                 pointr(j)=i
                end if
4       continue
        npick = j
        else
                do 22 i=1,npt
                        ipick(i) = 1
                        pointr(i) = i
22              continue
                npick = npt
        end if
        write (6,*) 'select is ',npick
c
        if (npick.eq.0) then
         if (debug) write(stdo,*)' ipick = ',(ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
        end if
        if (urcrd .eq. 99 ) then
                write(*,*)' Sorry, initial coordinate file not provided'
                write(*,*)' missing  - urcrd number '
                write(*,*)' Quitting...., better luck next time'
                level = 1
                call alert(name,namel,'Problem in file',15,level)
        else if (mod(igrid,2) .eq. 0) then
                write(*,*)'  number of grid points better be odd '
                write(*,*)' Input is igrid = ',igrid
                level = 0
                call alert(name,namel,'Problem in grid',15,level)
        end if
        write(stdo,100)nstep,npri,igrid,stdo,uwcrd,urcrd,crdstyl
     1  ,nwcrd
        write(stdo,101)tolg,ntest,gamma,rho,lambda
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR CHAIN MINIMIZATION: ',//,
     1   1x,' number of minimization steps: ',i7,/,
     2   1x,' print each ',i7,' steps',/,
     3   1x,' number of grid points of the chain : ',i5,/,
     4   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     5   1x,' initial coordinates are read from unit : ',i5,/,
     6   1x,' Coordinate are read in ',a4,' style ',/,
     7   1x,' coordinates are written at step interval : ',i5)
101     format(1x,' gradient max. error ',e10.3,/,
     1   1x,' period for checking constraints ',i5,/,
     2   1x,' Monomer - monomer parameters: ',/,
     3   1x,' i,i+1 force constant ',f10.5,/,
     4   1x,' i,i+2 repulsion strength    ',f10.5,/,
     5   1x,' i,i+2 repulsion range       ',f10.5)


       if (eenmyes) then
        write(stdo,102)enm_cutoff,enm_alpha,enm_beta,enm_gamma
102     format(/,1x,' PARAMETERS FOR ELASTIC NETWORK: ',//,
     1   1x,' Cutoff: ',f10.5,/,
     2   1x,' Alpha:  ',f10.5,/,
     3   1x,' Beta : ',f10.5,/,
     4   1x,' Gamma : ',f10.5)
       end if

        call init_wre(stdo)

c
c set the pointer to selected particles. pointr(i) is the position
c of atom number i in the normal all atom array
c
        npick3 = 3*npick
        npt3   = 3*npt
        middle = igrid/2*npt3
c
c initialize sigmav to zero
c
        do 5 i=1,6
                sigmav(i) = 0.d0
5       continue

        if (gbsabool) call make_rborn
c
c get initial coordinates (call rchain)
c depending on "crdstyl", the following options exist
c (i)   DYNAmics files (single precision unformatted compatible with CHARMm)
c (ii)  PATH files (binary, double precision code to maintain 
c        accurate coordinates
c (iii) INIT. Reading formatted coordinates file for reactants and
c             products and generating the rest of the path by
c             linear interpolation
c (iv)  INTRpolate. Given a low resolution path in PATH format, add
c             structures in between to refine the path.
c 
        call rchain(urcrd,r,igrid,crdstyl,pointr,npick)

c
c get center of mass constraints and orthonormalize all constraints
c generating constraints for middle coordinate set that will be used
c by everybody. dv is used here as a temporary vector.
c

c
c@ for chmin calculations all the masses should be set to one
c
        do i=1,npt
                ptms(i) = 1.d0
        end do 

         call comc(r(middle+1),ptms,npt,pointr,npick,grdcmx,
     1          grdcmy,grdcmz,grdlx,grdly,
     2          grdlz,dv(1),dv(1+npt),
     3          dv(1+2*npt),dv(1+3*npt),
     4          dv(1+4*npt),sigma,debug)
        if (debug) then
                write(*,*) ' grdcmx '
                write(*,1000)(grdcmx(i),i=1,npick3)
                write(*,*) ' grdcmy '
                write(*,1000)(grdcmy(i),i=1,npick3)
                write(*,*) ' grdcmz '
                write(*,1000)(grdcmz(i),i=1,npick3)
                write(*,*) ' grdlx '
                write(*,1000)(grdlx(i),i=1,npick3)
                write(*,*) ' grdly '
                write(*,1000)(grdly(i),i=1,npick3)
                write(*,*) ' grdlz '
                write(*,1000)(grdlz(i),i=1,npick3)
1000            format(1x,7(f10.5,1x))
        end if

c overlap all structures with respect to the middle one

        if (lap) then
         imx1 = middle+1
         do i=1,igrid
          ix1 = (i-1)*npt3+1
          if (ix1.ne.imx1) 
     &       call rmsd_weight(npt,r(imx1),r(ix1),rms,.false.,ptms)
         end do
        end if
        tolg = tolg*tolg*(igrid-2)*npt3
c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp
c
c       And now we are ready to call the CHain MIN1 routine
c
        if (.not.anneal) then
         call chmin1(igrid,npri,nstep,pointr,ipick,
     1          npick,d0,e0,e1,grdcmx,grdcmy,grdcmz,grdlx,
     2          grdly,grdlz,r,dv,tolg,udata,uwcrd,
     3          nwcrd,ntest,debug,scalar,sigma,sigmav,gamma,
     4          rho,lambda,fixend,nlist,estred)
        else
         call chnneal(igrid,npri,nstep,temper,pointr,ipick,
     1          npick,d0,e0,e1,grdcmx,grdcmy,grdcmz,grdlx,
     2          grdly,grdlz,r,dv,udata,uwcrd,
     3          nwcrd,ntest,debug,scalar,sigma,sigmav,gamma,
     4          rho,lambda,fixend,nlist)
        end if

        stop
        end
