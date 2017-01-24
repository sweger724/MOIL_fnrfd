        subroutine init_dyna()
c
c initializing as many as possible of the dyna variables
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/COUPLED_DYNA.BLOCK'

c local

c getc - a function  to get a character
c name - name of current routine
c namel - the length (in characters of the name of the 
c               current routine
c crdtyp - the type of coordinates to be read
c urcrd  - unit number to read coordinates
c urtet  - unit number to read tether constraints
c urvel  - unit number to read velocities
c ucon   - unit number for connectivity file
c ucon1 & ucon2 - two special connectivity files
c               to model two electronic energy curves and
c               photodissociation
c geti   - integer function to get an integer
c of     - an integer function to Open a File. It returns the unit
c               number
c getd   - a function returning a double precision number
c find   - a logical function. set to .true. if find expression in line
c fopen  - is set to true if unit is opened
c 

        character*9 name
        integer namel

        integer i,j,k,n
        integer ucon,urvel,urcrd
        integer lc,level
        integer stdin, stdout

        double precision tmp

        logical fopen


        name='init_dyna'
        namel=9

c initialilze some general variables:
        call init_con()

        lc     = 5
c nori - if true DO NOT do reorientation during dynamics
c important when tether
c
        nori   = .false.
        eqms   = .false.

        lcent  = .false.

        debug  = .false.
        shakl  = .false.
        shakb  = .false.
        shakm  = .false.
        freeze = .false.
        boltz  = .true.
        nocut  = .false.
        sdyes  = .false.
        no_scaling = .false.
        hvdw0      = .true.


c DYNAMICS DEFAULT PARAMETERS
        nstep = 1
        neqstep = 1
        ninfo = 1
        ncoor = 0
        newv  = 0
        n_tether = 0
        nrigi    = 10
        nvelo = 0
        nlist = 1
        nscalv  = 0
        ntemp   = 1
        start_dyna = 1
        dt      = 0.001d0
        tempi(1) = 300.d0
        tempf(1) = 300.d0
        tolcons  = 1.d-7
        epshak   = 1.d-7
        epshakv  = 1.d-7
        itershak = 100
        fmax     = -1.d0
        symanneal= .false.
      wfly = .false.
C PRLL default
        first_gather = .true.


C ENERGY DEFAULT PARAMETERS


        call init_ef()

C       Initial data for Brownian Dynamics
C
        LDgamma=0.0d0
        bdcf0=1.0d0
        bdcf1=1.0d0
        bdcf2=1.0d0

C Set Standard input and otput for each paralel run
        stdin=5
        stdout=6
        write(*,*) 'procID: ',procID
        call open_inout(stdin,stdout,procID) 

C READ STANDARD INPUT
        call line_loop(ucon,urcrd,urvel)

c if spherical solvent boundary potential
        if (qssbp) call set_ssbp()

c PROCESS AND TEST SOME INPUT
        call set_cutoffs()

c
c if Ewald summation for periodic boundary conditions
c
c jmjm
c
        if (ewaldyes) then
c          write(stdo,*) 'Calling ewald_init '
           call ewald_init()
        end if
c if virtual particles
        if (vp_flag) then
           call vp_init()
        end if
c
c jmjm end
c

c check symmetry annealing
        if (symanneal.and..not.esymyes) then
          level = 1
          call alert(name,namel,'symm problem!',13,level)
        endif

        if (ncnst.gt.0) write(stdo,118) ncnst
118     format ('Number of TORS constraints is:',i7)


c check that files were opened
        if (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        else if (.not.fopen(uwvel)) then
         level = 1
         call alert(name,namel,'uwvel not opened',16,level)
        end if


c get velocities
        if (.not.boltz) then
c read velcoities like charmm coordinate file
                 call getvel(urvel,'CHARM')
        else
c initialized velocties
c irand - integer to initialize the randomizer
c #### this section can be parallelized too... later....
                call velinit(tempi,ntemp,tpo)
        end if

c factor out center of mass velocity and
c factor out the rigid body rotation by orienting the structure against
c the original configuration and also by calculating the gradient of
c the constraints on rotation and projecting out the components
c of the velocities along these directions
c
        if (.not.(esymyes.or.nori)) then

c ####prbm will be parallelized eventually
         call prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)


        end if

c write an energy head.
        if (my_pe.eq.0) then
         call init_wre(stdo)
c ------------------------------------------------
         if ( emyes0 ) then
          write (stdo,*)'D=',(D(n),n=1,nmb)
          write (stdo,*)'alpha=',(alpha(n),n=1,nmb)
         end if
         if (specl) then
          write (stdo,*)'rcut=',(rcut(n),n=1,nmb)
          write (stdo,*)'lamda=',(lamda(n),n=1,nmb)
         end if
         if (repyes0) then
          write (stdo,*)'Arep=',(Arep(n),n=1,nmb)
          write (stdo,*)'beta =',(beta1(n),n=1,nmb)
         end if
c ---------------------------------------------------------
         if(switch) then
          print*,'delt = ',delt
          print*,'Rcro = ',Rcros
          print*,'dRcr = ',dRcros
          print*,'Force = ',Force
         end if
c --------------------------------------------------------
        end if
c ---------------------------------------------------------
c
c prepare multiplicative constants, first calculate dt in KAM

        dt    = dt/tconv
        twodt = dt*0.5d0
        dt2   = dt*dt*0.5d0
        do 7 i=1,npt
                j = nofreez(i)
                factor1(i) = dt2*invms(j)
                factor2(i) = twodt*invms(j)
7       continue
c
c prepare equilibration rate (if relevant
        if (neqstep.gt.1) then
          tmp = 1.d0/dfloat(neqstep-1)
          do 8 k=1,ntemp
             dtemp(k) = (tempf(k)-tempi(k))*tmp
             tempi(k) = tempi(k) + (start_dyna - 1) * dtemp(k)
8         continue
        else if (neqstep.eq.1) then
          do 81 k=1,ntemp
             dtemp(k) = (tempf(k)-tempi(k))
81        continue
        end if

c
c calculate the box size step to use to squeeze the box if symanneal is true
c the box size will be changed with each call to making new list
c and only during equilibration
c
        if (symanneal) then
                tmp = dfloat(nlist)/dfloat(neqstep)
                dxtra=(syma2-a)*tmp
                dytra=(symb2-b)*tmp
                dztra=(symc2-c)*tmp
                write (stdo,*) 'box size squeezing increments: (x,y,z)'
                write (stdo,*) dxtra,dytra,dztra
        endif

        urst_crd = -1
        cont_yes = .false.

        return

        end

