       program dim_sample
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
        include 'COMMON/MASSFAC.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
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

        integer urcrd,ucon, urint,urvel,ufpt
        integer npick
        integer imx1,ix1

        logical fopen,lap
        logical select, amid_true

        double precision scalar(6),sigma(6),sigmav(6)
        double precision clo,DD(lgrid),cosine(lgrid),SumD
        integer i,j,k,l,npick3,middle,amount,n
        double precision aa(3,maxpt),bb(3,maxpt)
        double precision normA,normB,a_b, sum, CGtorsion,DistAngle
        real ranx

        integer reached(maxnodes), reach, Nreach
        integer ntry, try, utmpx, i1,i2,i3,i4
        logical drawn(20000)

        data urcrd,ucon,uwcrd,uwvel,urint,urvel/6*99/
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
        MASSWEIGHT=0
        irand = -1

        chmin=.false.
        select = .false.
        paral = .false.
        my_pe = 0
        numprocs =1
        debug  = .false.
        name   = 'path '
        namel  = 5
        
        call init_paral()
        
        !call init_io()
        call open_inout(stdi,stdo,my_pe)
c open scratch file for line manipulation
c
        jnkf = 1
        open (unit=jnkf,status='scratch')

        if (my_pe .eq. 0) then
          call getunit(utmpx)
          open(unit=utmpx,file="tmp_xxx.dcd",form='unformatted'
     &         ,status='unknown')
        end if
     
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

c        conformation RMSD cutoff
        cutoff = 0.2   
c      maximum deviation from the initial path        
        searchLimit = 4.d0  

        igrid = 1
        myCell = 1
        myCell2 = 2
        
c       structures per expanded node        
        Stotal = 20    

C     Read inpput file
        call line_loop(dtopt,
     &       temper,constant_tmp,clo,interpolate,lap,
     &       select,pointr,ucon,urcrd,urint,urvel,amid_true,ufpt)

C  change random number seed for each processor       
        irand = irand 
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
        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if ((.not.(fopen(urcrd)) )) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (my_pe.eq.0 .and. .not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
        else if (my_pe.eq.0 .and. .not.fopen(uwvel)) then
         level = 1
         call alert(name,namel,'uwvel not opened',16,level)
        else if (.not.fopen(urint)) then
         level = 1
         call alert(name,namel,'urint not opened',16,level)
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
        do 21 i=1,npt
                zerofrz(i) = 1
21      continue

c Pick atoms (note that chain/chmin is currently NOT working
c with pick, hoping to fix this soon)
c
        if (.not.select) then
                do 22 i=1,npt
                        mmm(i) = ptms(i)
22              continue
        else
           call rline(name,namel,stdi)
           call pick(ipick,npick)
           do i=1,npt
           if (ipick(i).eq.0) then
            mmm(i) =0
           else
            mmm(i) = ptms(i)
          end if
          enddo
          
        end if

c
        npick=0
        do 3 i=1,npt
            if(ipick(i).gt.0) then   
               npick=npick+1
               rms_pick(npick) = i
            end if 
3       continue
        if (npick.eq.0) npick = npt
        write(stdo,100)nstep,ninfo,stdo,uwcrd,urcrd,crdstyl,nwcrd
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR SEARCH: ',//,
     1   1x,' number of minimization steps: ',i7,/,
     2   1x,' print each ',i7,' steps',/,
     4   1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     5   1x,' initial coordinates are read from unit : ',i5,/,
     6   1x,' Coordinate are read in ',a4,' style ',/,
     7   1x,' coordinates are written at step interval : ',i5)

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
c set the pointer to ALL particles. pointr(i) is the position
c of atom number i in the normal all atom array
c
        j=0
        do 4 i=1,npt
                j=j+1
                pointr(j)=i
4       continue
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

        ! get set of cells
        first = .true.
        call rchain(urcrd,crdstyl)


        ! read starting structure in the interface
        call getcrd(urint,'CHARM')
C
c get center of mass constraints and orthonormalize all constraints
c generating constraints for middle coordinate set that will be used
c by everybody. dv is used here as a temporary vector.
c

          call comc(coor(1,1),ptms,npt,pointr,npick,grdcmx,
     1          grdcmy,grdcmz,grdlx,grdly,grdlz,sigma,debug)

c
c initiate GB calculations
      if (gbsabool) then
         call make_rborn
       endif

       if (ewaldyes) call ewald_init()

       
       if (vp_flag) call vp_init()


        cgrd2 = 0.d0

        Ncells = igrid-2
        j = 1
        do i =1, igrid
          if (i .ne. myCell .and. i .ne. myCell2) then
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
          do i = 1, Ncells+2
          ReducedCoor(j,i)=CGtorsion(i1,i2,i3,i4,centers(1,(i-1)*npt+1))
            write(6,*)"Torsion",j,"in ",i,"structure is"
     &                ,ReducedCoor(j,i)*pi180
          end do
        end do
     
c         align all centers to the 1.st one

c DO NOT ALIGN if esymyes is on. RE
c
        if (.not.esymyes) then
         do i =2, igrid
         k = (i-1)*npt+1
         call rmsd_weight(npick,centers(1,1),centers(1,k),rms,
     & .false.,mmm)
         end do
        end if
c

        normA = 1000.d0
        do i =1, igrid
          k = (myCell-1)*npt+1
          if (myCell .ne. i) then
            if (Nreduced .ne. 0) then
            rms = DistAngle(ReducedCoor(1,myCell),
     6                      ReducedCoor(1,i),Nreduced)
            else
              l = (i-1)*npt+1
         call rmsd_weight(npick,centers(1,l),centers(1,k),rms,
     &     .true.,mmm)
              rms = rms **2
            end if
            write(6,*) "Square of distance between:",i,myCell, rms
            if (rms .lt. normA ) normA = rms
          end if
        end do
        write(6,*)"Minimal squared distance:", normA
        dx = normA

        if (eenmyes) call enm_lists(coor(1,1),coor(1,1))
        if (eCGyes)  call CGinit()
        call init_dyna()
         
        call dim_sample_MD(reach,utmpx)
          
        call parallel_end()
        stop
        end
