       program search
       implicit none
c
c calculate path using SDEL algorithm
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/PATH2.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'
        include 'COMMON/FREADY.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

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
c  d0 - distances between intermediates (i,i+1)
c  grdcm[x-z] = derivative of the center of mass with respect to [x-z]
c  grdl[x-z]  = derivatives of infitesimal rotations with respect to [x-z]
c pointr pick1 - selection of particle
c
c total local space (in vectors required):
c               24*maxpt*(lgrid+6) +  4*lgrid 
c       
c
        double precision d0(lgrid+1)
        double precision grdcmx(3*maxpt),grdcmy(3*maxpt)
        double precision grdcmz(3*maxpt)

        integer pointr(maxpt)
        integer interpolate
        double precision temper,dtopt,cgrd2
        logical constant_tmp

        character*5 name
        integer namel, ierr

        integer urcrd,ucon
        integer npick
        integer imx1,ix1

        integer npri
        logical fopen,lap
        logical select, amid_true

        double precision scalar(6),sigma(6),sigmav(6)
        double precision clo,DD(lgrid),cosine(lgrid),SumD
        integer i,j,k,l,npick3,middle,amount,n
        double precision aa(3,maxpt),bb(3,maxpt)
        double precision normA,normB,a_b

        data urcrd,ucon,uwcrd/3*99/
c
c General initialization
c

        cycleno = 0
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
        irand = -1
        npick = 0

        chmin=.false.
        select = .false.
        paral = .false.
        my_pe = 0
        numprocs =1
        debug  = .false.
        name   = 'path '
        namel  = 5
c ....initialization of the parallel machine
        call init_path_paral()
        write(6,*) 'paral = ',paral,' processors = ',numprocs
        call open_inout(stdi,stdo,my_pe)
c open scratch file for line manipulation
c
        jnkf = 1
        open (unit=jnkf,status='scratch')
c
c chain default parameters
c
        Random_velocities = .true.
        crdstyl = 'CHAR'
        nstep   = 1
        npri    = 1
        nwcrd   = 0
        temper  = 300.d0
        constant_tmp = .false.
        clo=0.0d0
        
        gamma = 100.d0
        dtopt=1.0d-5
c
c initiate parameters for SGB
c
        sgba = 0
        gbsu = 0
        gbsabool = .false.
C     Set our default interpolation mode to 0.
        interpolate = 0
c do an overlap of the structures with respect to each other?
c
        lap = .false.
c
c energy default parameters
c
        call init_ef()
        nocut  = .false.
        amid_true = .false.

        nlist  = 1

C default ELastic network model parameters
       Hamilt= 1.d-10

       enm_cutoff = 7.d0
       enm_alpha = 0.d0
       enm_beta = 1.d0
       enm_gamma = 1.d0 
       
C default READY parameters
        CG_cutoff = 13.5d0 
        HB_cutoff = 8.d0

        cutoff = 0.2   ! conformation RMSD cutoff
        searchLimit = 4.d0   ! maximum deviation from the initial path

        igrid = 1
        
        Stotal = 20    ! structures per expanded node

C     Read inpput file
        call line_loop(dtopt,npri,pointr,npick,
     &       temper,constant_tmp,clo,interpolate,lap,
     &       select,ucon,urcrd,amid_true)

C  change random number seed for each processor       
        irand = irand - my_pe
        call RLUXGO(223,irand,0,0)

        if (nwcrd.eq.0) nwcrd = nlist


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
        if (my_pe.eq.0 .and. .not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (my_pe.eq.0 .and. (.not.(fopen(urcrd)) )) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (my_pe.eq.0 .and. .not.fopen(uwcrd)) then
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
                cutmono2 = cutebig2*1.44d0
        else
                cutmono2 = cutmono2*cutmono2
        end if

c        if (amid_true) call amid()

c initialize masses of CA atoms
        if (eCGyes) call initCGmass()
        if (eenmyes) then
          do  i=1,npt
            ptms(i)=13
          end do
        end if
        
c initialze no freez vector
        inofrz = npt
        do 21 i=1,npt
                zerofrz(i) = 1
21      continue

c Pick atoms (note that chain/chmin is currently NOT working
c with pick, hoping to fix this soon)

        if (npick .eq. 0) then
          do i=1,npt
            ipick(i) = 1
            pointr(i) = i
          end do
          npick = npt
        end if


        if (npick.eq.0) then
         if (debug) write(stdo,*)' ipick = ',(ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
        end if
        write(stdo,100)nstep,npri,stdo,uwcrd,urcrd,crdstyl,nwcrd
        write(stdo,101)gamma
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR SEARCH: ',//,
     1   1x,' number of minimization steps: ',i7,/,
     2   1x,' print each ',i7,' steps',/,
     4   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     5   1x,' initial coordinates are read from unit : ',i5,/,
     6   1x,' Coordinate are read in ',a4,' style ',/,
     7   1x,' coordinates are written at step interval : ',i5)
101     format(1x,' Monomer - monomer parameters: ',/,
     2   1x,' i,i+1 force constant ',f10.5)

        if (eenmyes) then
        write(stdo,102)enm_cutoff,enm_alpha,enm_beta,enm_gamma
102     format(/,1x,' PARAMETERS FOR ELASTIC NETWORK: ',//,
     1   1x,' Cutoff: ',f10.5,/,
     2   1x,' Alpha:  ',f10.5,/,
     3   1x,' Beta : ',f10.5,/,
     4   1x,' Gamma : ',f10.5)
       end if


        call init_wre(stdo)

        npick3 = 3*npick
        npt3   = 3*npt

c
c initialize sigmav to zero
c
        do 5 i=1,6
                sigmav(i) = 0.d0
5       continue

c
c.....output number of processors used
      if (first.and.last) then
         write(stdo,*)' *** ONLY ONE PROCESS IS USED '
      else
         write(stdo,*)' *** ',numprocs,' PROCESSES ARE USED '
      end if

c get initial coordinates (call rchain)
c (iii) INIT. Reading formatted coordinates file for reactants and
c             products and generating the rest of the path by
c             linear interpolation
        call rchain(urcrd,crdstyl)

c
c initiate GB calculations
      if (gbsabool) then
         call make_rborn
       endif
C
c get center of mass constraints and orthonormalize all constraints
c generating constraints for middle coordinate set that will be used
c by everybody. dv is used here as a temporary vector.
c
        if (first) then
      call rmsd_weight(npt,r_initial(1,1),r_final(1,1),rms,.false.,ptms)
          call comc(coor(1,1),ptms,npt,pointr,npick,grdcmx,
     1          grdcmy,grdcmz,grdlx,grdly,
     2          grdlz,sigma,debug)
        end if

       if (ewaldyes) call ewald_init()
       
       if (vp_flag) call vp_init()


c       And now we are ready to call the CHain MIN1 routine
c
        cgrd2 = 0.d0
         call search_2(temper,constant_tmp,pointr,
     1          npick,d0,grdcmx,grdcmy,grdcmz,
     2          stdo,uwcrd,
     3          debug,scalar,sigmav,cgrd2
     5          ,dtopt,clo,npri)
      

c Do analysis of congestion
       write(6,*) 'Processor',my_pe
       
C        call MPI_Buffer_detach( mpi_buffer, 30*maxpt,ierr );
        call comm_exit(ierr)

        stop
        end
