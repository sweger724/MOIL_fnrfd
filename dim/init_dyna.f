        subroutine init_dyna()
        implicit none
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
        include 'COMMON/FREADY.BLOCK'
        include 'COMMON/LD.BLOCK'

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

        double precision tmp

        logical fopen


        name='init_dyna'
        namel=9


c       checkpointing
        resume = .false.
        nchkpt = 200


c if spherical solvent boundary potential
        if (qssbp) call set_ssbp()

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

c factor out center of mass velocity and
c factor out the rigid body rotation by orienting the structure against
c the original configuration and also by calculating the gradient of
c the constraints on rotation and projecting out the components
c of the velocities along these directions
c
        if (.not.(esymyes.or.nori)) then

         call prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)


        end if

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

c convert the friction constant to KAM units
        LDgamma=LDgamma * tconv

c prepare equilibration rate (if relevant
        if (neqstep.gt.1) then
          tmp = 1.d0/dfloat(neqstep/2-1)
          do 8 k=1,ntemp
             dtemp(k) = (tempf(k)-tempi(k))*tmp
             tempi(k) = tempi(k) + (start_dyna - 1) * dtemp(k)
8         continue
        else if (neqstep.eq.1) then
          do 81 k=1,ntemp
             dtemp(k) = (tempf(k)-tempi(k))
81        continue
        end if


        itershak = 100
        epshak   = 1.d-7
        epshakv  = 1.d-7

        return

        end

