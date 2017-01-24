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
        include 'COMMON/LD.BLOCK'
        include 'COMMON/PT.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PRESS.BLOCK'

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

        integer ith,ia,ja,ka

        double precision tmp

        logical fopen
        logical match1,match2,match3


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
        pdebug = .false.
        shakl  = .false.
        shakb  = .false.
        shakm  = .false.
        freeze = .false.
        boltz  = .true.
        nocut  = .false.
        sdyes  = .false.
        no_scaling = .false.
        hvdw0      = .true.

        pressON = .false.

c DYNAMICS DEFAULT PARAMETERS
        nstep = 1
        LDgamma = 60.d0   ! Langevin dynamics default friction coefficient
        neqstep = 1
        ninfo = 1
        ncoor = 0
        newv  = 0
        n_tether = 0
        nrigi    = 10
        nvelo = 0
        nlist = 1
        nscalv  = 20
        ntemp   = 1
        start_dyna = 1
        irand   = 1
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
c       LDgamma=0.0d0
c       bdcf0=1.0d0
c       bdcf1=1.0d0
c       bdcf2=1.0d0


c default parameters for muta
          muta = .FALSE.
          in_lambda=0.0
          f_lambda=0.0
          n_step_lambda=0
          twait=0
          BandF=.FALSE.
          do i=1,lambda_max
           D_E_L(i)=0.0d0
           sD_E_L(i)=0.0d0
          enddo
c default parameters for urey-bradley bonds
        nurey = 0


c       checkpointing
        resume = .false.
        nchkpt = 200

 

C READ STANDARD INPUT
        call line_loop(ucon,urcrd,urvel)
! Are we running parallel tempering?
        if (nswap .ne. 0 ) PT_run = .true.
        old_step = 0 ! last step that communication accured in PT

C  initialize random numbers
        write(6,*) "irand:", irand
        call RLUXGO(223,irand,0,0)
       if (PT_run) call prepare_random_numbers(0)

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

        if (eCGyes) then
          call CGinit()
        endif

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


c #### this section can be parallelized too... later....
c get velocities
        if (.not.boltz) then
c read velcoities like charmm coordinate file
                 call getvel(urvel,'CHARM')
        else
c initialize velocties
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
          write (stdo,1000)(D(n),n=1,nmb)
          write (stdo,1001)(alpha(n),n=1,nmb)
1000      format(1x,'D=',(1x,f9.3))
1001      format(1x,'alpha=',(1x,f9.3))
         end if
         if (specl) then
          write (stdo,1002)(rcut(n),n=1,nmb)
1002      format(1x,'rcut=',(1x,f9.3))
          write (stdo,1003)(lamda(n),n=1,nmb)
1003      format(1x,'rcut=',(1x,f9.3))
         end if
         if (repyes0) then
          write (stdo,1004)(Arep(n),n=1,nmb)
1004      format(1x,'Arep=',(1x,f9.3))
          write (stdo,1005)(beta1(n),n=1,nmb)
1005      format(1x,'beta=',(1x,f9.3))
         end if
c ---------------------------------------------------------
         if(switch) then
          write(*,1006)delt,Rcros,dRcros,Force
1006      format(1x,' delt Rcro dRcr Force ',4(f9.3))
         end if
c --------------------------------------------------------
        end if
c ---------------------------------------------------------
c
c prepare multiplicative constants, first calculate dt in KAM

        dt    = dt/tconv
        twodt = dt*0.5d0
        dt2   = dt*dt*0.5d0
        do 7 i=1,inofrz
                j = nofreez(i)
                factor1(i) = dt2*invms(j)
                factor2(i) = twodt*invms(j)
7       continue

c       convert friction constant from ps^-1 to KAM
        LDgamma = LDgamma * tconv

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
c@
c                write(stdo,*) ' neqstep nlist syma2 symb2 symc2 a b c '
c                write(stdo,*)neqstep,nlist,syma2,symb2,symc2,a,b,c
c
                write (stdo,*) 'box size squeezing increments: (x,y,z)'
                write (stdo,*) dxtra,dytra,dztra
        endif

        urst_crd = -1
        cont_yes = .false.


        if (.not.fopen(uchx)) then
           level = 1
           call alert(name,namel,'uchx not opened',16,level)
        end if
        if (.not.fopen(uchv)) then
           level = 1
           call alert(name,namel,'uchv not opened',16,level)
        end if


c       resume from checkpoint
        if ( resume ) then
           read(uchx) start_dyna,e_total,((coor(j,i),i=1,npt),j=1,3)
           read(uchv) start_dyna,e_total,((velo(j,i),i=1,npt),j=1,3)
        end if

c if there are urey-bradley bonds...

       DO i=1,nangl
        kangl_urey(i) = 0.d0
       ENDDO

       DO 100 i=1,nurey
        DO ith=1,nangl
         match1 = .false.
         match2 = .false.
         match3 = .false.
         ia = iangl1(ith)
         ja = iangl2(ith)
         ka = iangl3(ith)
         if (ia.eq.urey1(i).or.ia.eq.urey2(i).or.
     1       ia.eq.urey3(i)) match1 = .true.
         if (ja.eq.urey1(i).or.ja.eq.urey2(i).or.
     1       ja.eq.urey3(i)) match2 = .true.
         if (ka.eq.urey1(i).or.ka.eq.urey2(i).or.
     1       ka.eq.urey3(i)) match3 = .true.
         if (match1 .and. match2 .and. match3) then
c          index_urey(i) = ith
          write (stdo,*) 'selected angle',ith
          goto 101
         endif
        ENDDO
         write (stdo,*) i,urey1(i),urey2(i),urey3(i)
         call alert(name,namel,'angle not found',15,1)
101     CONTINUE
c        kkk = kangl(index_urey(i))
c        write (*,*) kkk
c this function returns the proper urey-bradley k and puts kangl=0
        call urey_init(ith)
        write (*,*) kangl(ith),angleq(ith)
        write (*,*) kangl_urey(ith) ,angleq_urey(ith)
c        kangl(index_urey(i)) = kkk
100    CONTINUE

c       do ith=1,nangl
c        write (*,*) 'regular',ith,kangl(ith),angleq(ith)
c        write (*,*) 'urey',ith,kangl_urey(ith) ,angleq_urey(ith)
c       enddo


        return

        end

