       program dim_prepare
       implicit none
c
c start a simulation in a cell center and terminates on an interface
c with another cell
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
        include 'COMMON/OVERLAP.BLOCK'

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

        integer interpolate
        double precision temper,dtopt,cgrd2, DistAngle
        logical constant_tmp

        character*5 name
        integer namel, ierr

        integer urcrd,ucon,urint,urvel,ufpt
        integer npick
        integer imx1,ix1

        logical fopen,lap
        logical select, amid_true

        double precision scalar(6),sigma(6),sigmav(6)
        double precision clo,DD(lgrid),cosine(lgrid),SumD
        integer i,j,k,l,npick3,middle,amount,n,pointr(maxpt)
        double precision aa(3,maxpt),bb(3,maxpt)
        double precision normA,normB,a_b, CGtorsion

        integer reached(maxnodes), reach, Nreach 
        integer try, i1,i2,i3,i4,tors(maxpt,4)
        double precision dummy(3,maxNodes*maxpt)

        data urcrd,ucon,uwcrd,urint,urvel/5*99/
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

        chmin=.false.
        select = .false.
        paral = .false.
        my_pe = 0
        numprocs =1
        debug  = .false.
        name   = 'path '
        namel  = 5
        
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
        myCell = 1
        
        Stotal = 20    ! structures per expanded node

        call init_paral()

C     Read inpput file
        call line_loop(dtopt,
     &       temper,constant_tmp,clo,interpolate,lap,
     &       select,pointr,ucon,urcrd,urint,urvel,amid_true,ufpt)

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

        if (amid_true) call amid()

c initialize masses of CA atoms
        if (eCGyes) call initCGmass()
        if (eenmyes) then
          do  i=1,npt
            ptms(i)=13
          end do
        end if
        
c initialze no freez vector
        inofrz = npt
        do i=1,npt
          zerofrz(i) = 1
        end do

c        do i=1,npt
c          ipick(i) = 1
c        end do

c
        npick=0
        do 3 i=1,npt
            if(ipick(i).eq.0) then
               npick=npick+1
               mmm(i) = ptms(i)
             else
               mmm(i) = 0.d0
            end if
3       continue
        if (npick.eq.0) npick = npt 
     
         if (npick.eq.0) then
         if (debug) write(stdo,*)' ipick = ',(ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
        end if
        write(stdo,100)nstep,ninfo,stdo,uwcrd,urcrd,crdstyl,nwcrd
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

        write(6,*)"allmass:", allmass

        call init_wre(stdo)

        npick3 = 3*npick
        npt3   = 3*npt

c
c initialize sigmav to zero
c
        do 5 i=1,6
                sigmav(i) = 0.d0
5       continue

c get initial coordinates (call rchain)
c (iii) INIT. Reading formatted coordinates file for reactants and
c             products and generating the rest of the path by
c             linear interpolation
        first = .true.
        call rchain(urcrd,crdstyl)

c
c initiate GB calculations
      if (gbsabool) then
         call make_rborn()
       endif

       if (ewaldyes) call ewald_init()
       
       if (vp_flag) call vp_init()


        cgrd2 = 0.d0


        Ncells = igrid-1
        j = 1
        do i =1, igrid
          if (i .ne. myCell) then
            cell(j) = i
            reached(j) = 0
            j = j + 1
          end if
        end do

      
        do j = 1, Nreduced
          i1 = torsions(j,1)
          i2 = torsions(j,2)
          i3 = torsions(j,3)
          i4 = torsions(j,4)
          do i = 1, Ncells+1
           
          ReducedCoor(j,i)=CGtorsion(i1,i2,i3,i4,centers(1,(i-1)*npt+1))
c         write(6,*)"RRR",ReducedCoor(j,i)*pi180,j,i,(i-1)*npt+1
          end do 
        end do
      
   
        ! align all centers to the 1.st one
c        do i =2, igrid
c         k = (i-1)*npt+1
c         call rmsd_weight(npt,centers(1,1),centers(1,k),rms,.false.,mmm)
c        end do
    
c        dummy coordinates to use in rmsd_weight instead of centers
        do i=1,igrid*npt 
        dummy(1,i)=centers(1,i) 
        dummy(2,i)=centers(2,i) 
        dummy(3,i)=centers(3,i) 
        enddo  


       !normA = 1000.d0
        !do i =1, igrid
        !  k = (i-1)*npt+1
        !  do j = i+1, igrid
        !    if (Nreduced .ne. 0) then
        !     rms = DistAngle(ReducedCoor(1,i),ReducedCoor(1,j),Nreduced)
        !    else
        !      l = (j-1)*npt+1
        !      call Distance(centers(1,k),centers(1,l),rms)
        ! call rmsd_weight(npt,centers(1,l),centers(1,k),rms,.false.,mmm)
        !    end if
        !    write(6,*) "RRR:",i,j, rms
        !    if (rms .lt. normA ) normA = rms
        !  end do
        !end do
       !write(6,*)"Minimal distance:", normA
       !stop
        normA = 1000.d0
        do i =1, igrid
          k = (myCell-1)*npt+1
          if (myCell .ne. i) then
            if (Nreduced .ne. 0) then
            rms = DistAngle(ReducedCoor(1,myCell),
     6                      ReducedCoor(1,i),Nreduced)
            else
              l = (i-1)*npt+1
       call rmsd_weight(npick,centers(1,l),dummy(1,k),rms,.true.,mmm)
              rms = rms * rms
            end if
            write(6,*) "RRR:",i,myCell, rms
            if (rms .lt. normA ) normA = rms
          end if
        end do
        write(6,*)"Minimal distance:", normA
        dx = normA    
 
       !normA = 1000.d0
       !   do j = 1, Ncells
       !     rms = Distance(ReducedCoor(1,i),ReducedCoor(1,j),Nreduced)
       !     write(6,*) "RRR:",i,j, rms
       !     if (rms .lt. normA ) normA = rms
       !   end do
       ! end do
       !write(6,*)"Minimal distance:", normA
       !stop       

 
       ! normA = 1000;
       ! k = (myCell-1)*npt+1
       ! do j = 1, igrid
       !   if (j.ne.myCell) then
       !     l = (j-1)*npt+1
       !  call rmsd_weight(npt,centers(1,k),centers(1,l),rms,.false.,mmm)
       !      write(6,*) "RRR:",j, rms
       !      call Distance(centers(1,k),centers(1,l),rms)
       !      write(6,*) "RRR:",j, rms
       !     if (rms .lt. normA ) normA = rms
       !   end if
       ! end do
       !write(6,*)"Minimal distance:", normA

        do i = 1, npt
          do l = 1,3
             coor2(l,i) = centers(l,(myCell-1)*npt + i)
          end do
        end do
      
        if (eenmyes) call enm_lists(coor(1,1),coor(1,1))
        if (eCGyes)  call CGinit()
        call init_dyna()
     
        try = 0
        Nreach = 0

!        do while (.not.(try .eq. Stotal) .and. (Nreach .lt. Ncells))
        do while (.not.(try .eq. Stotal) .and. (Nreach .lt. 26)) ! 8
         
          do i = 1, npt
            do l = 1,3
             coor(l,i) = coor2(l,i)
            end do
          end do 
          
          call dim_prepare_MD(reach)  
          
          if (reach .ne. -1) then
            if (reached(reach).eq.0) then
              call Write2File(uwcrd,npt,coor(1,1),U(1),0,0)
              reached(reach) = 1
              Nreach = Nreach + 1
              write(6,*)" bounce",myCell, cell(reach)
              !call wener(stdo)
            else
              write(6,*)" again",myCell, cell(reach)
            end if
          end if  
          
          try = try + 1
          
        end do
        
        stop
        end
