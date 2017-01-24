       program path
       implicit none
c
c calculate path using SDEL algorithm
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
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'
        include 'COMMON/FREADY.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'
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
        double precision grdlx(3*maxpt),grdly(3*maxpt)
        double precision grdlz(3*maxpt)

        integer pointr(maxpt),ipick(maxpt)
        integer interpolate
        double precision rms
        double precision temper,dtopt,cgrd2
        logical constant_tmp

        character*5 name
        integer namel, ierr

        integer urcrd,ucon,uwcrd
        integer nstep,npick
        integer imx1,ix1
        integer nlist

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
        MASSWEIGHT=1
        irand = 1

        chmin=.false.
        select = .false.
        paral = .false.
        numprocs =1
        debug  = .false.
        name   = 'path '
        namel  = 5
c ....initialization of the parallel machine
        call init_path_paral()
        num_pes = 1
        my_pe = 0
        prll_on_off = .false.
        write(6,*) 'paral = ',paral,' processors = ',numprocs
        call open_inout(stdi,stdo,procID)
c open scratch file for line manipulation
c
        jnkf = 1
        open (unit=jnkf,status='scratch')
c
c chain default parameters
c
        Random_velocities = .true.
        igrid   = 10
        skpno = 1
        crdstyl = 'CHAR'
        nstep   = 1
        npri    = 1
        nwcrd   = 0
        temper  = 10.d0
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
       Hamilt= 0.d0

C     Read inpput file
        call line_loop(uwcrd,dtopt,npri,
     &       nstep,temper,constant_tmp,clo,interpolate,lap,
     &       select,nlist,ucon,urcrd,amid_true)

C  change random number seed for each processor       
        irand = irand + procID
c       initialize the random number generator
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
        if (first .and. .not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (first .and. (.not.(fopen(urcrd)) )) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
        else if (first .and. .not.fopen(uwcrd)) then
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
        else
                do 22 i=1,npt
                        ipick(i) = 1
22              continue
        end if

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
        write(stdo,100)nstep,npri,igrid,stdo,uwcrd,urcrd,crdstyl,nwcrd
        write(stdo,101)gamma
        write(stdo,*)' debug ? ',debug
100     format(/,1x,' PARAMETERS FOR CHAIN MINIMIZATION: ',//,
     1   1x,' number of minimization steps: ',i7,/,
     2   1x,' print each ',i7,' steps',/,
     3   1x,' number of grid points of the chain : ',i5,/,
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
         pseg=(igrid-2)/numprocs
         if (pseg.lt.1) then
            level=1
        call alert(name,namel,'At least 3 structures required',36,level)
         endif
      else
         write(stdo,*)' *** ',numprocs,' PROCESSES ARE USED '
         if((igrid-2).lt.numprocs) then
            level=1
            call alert(name,namel,'Too many processors!',36,level)
         endif
         if (mod(igrid-2,numprocs).ne.0) then
            level = 1
            if (mod(igrid-2,numprocs).gt.procID ) then
              pseg = (igrid-2)/numprocs + 1
            else
              pseg = (igrid-2)/numprocs
            endif
         else
           pseg=(igrid-2)/numprocs
         end if
         write(stdo,*)' *** ',pseg,' structures at this processor '
      end if


c
c get initial coordinates (call rchain)
c depending on "crdstyl", the following options exist
c (i)   DYNAmics files (for movies using QUANTA)
c (ii)  PATH files (binary, double precision code to maintain 
c        accurate coordinates
c (iii) INIT. Reading formatted coordinates file for reactants and
c             products and generating the rest of the path by
c             linear interpolation
c (iv)  INTRpolate. Given a low resolution path in PATH format, add
c             structures in between to refine the path.
c 
        if (paral) then
          if (pseg+2.gt.lgrid) then
             write (6,*) pseg+2,lgrid
             call alert(name,namel,'Lgrid exceeded',14,1)
           endif
        else
          if (igrid.gt.lgrid) then
            write (6,*) igrid,lgrid
             call alert(name,namel,'Lgrid exceeded',14,1)
          end if
        end if

        call rchain(urcrd,crdstyl,interpolate)
        if (first.and.last) pseg=(igrid-2)/numprocs
        call Get_Dls(pseg,npt,r,d0,ipick)

        if (first.and.last) then
           middle= igrid/2
        else 
           middle= (pseg/2) + 1
        endif
c
c initiate GB calculations
c
      if (gbsabool) then
         call make_rborn
       endif
C
c get center of mass constraints and orthonormalize all constraints
c generating constraints for middle coordinate set that will be used
c by everybody. dv is used here as a temporary vector.
c
         call comc(r(1,npt*(middle-1)+1),ptms,npt,pointr,npick,grdcmx,
     1          grdcmy,grdcmz,grdlx,grdly,
     2          grdlz,sigma,debug)
       if (lap) then 
         if (first) call rmsd_weight(npt,
     $                r(1,middle*npt+1),r_initial(1,1),rms,.true.,ptms)
         if (last) call rmsd_weight(npt,
     $                r(1,middle*npt+1),r_final(1,1),rms,.true.,ptms)
       end if

c overlap all structures with respect to the middle one

C       lap=.true.
        if (lap) then
          imx1 = middle*npt+1
          do 12 i=2,pseg+1
            ix1 = (i-1)*npt+1
            if (ix1.ne.imx1) 
     &          call rmsd_weight(npt,r(1,imx1),r(1,ix1),rms,.true.,ptms)
12        continue
        end if

       if (ewaldyes) call ewald_init()
       
       if (vp_flag) call vp_init()


c       And now we are ready to call the CHain MIN1 routine
c
        cgrd2 = 0.d0
         call path_2(nstep,temper,constant_tmp,pointr,ipick,
     1          npick,d0,grdcmx,grdcmy,grdcmz,grdlx,
     2          grdly,grdlz,stdo,uwcrd,
     3          debug,scalar,sigmav,nlist,cgrd2
     5          ,dtopt,clo,npri)
      

c Do analysis of congestion
       write(6,*) 'Processor',procID
       
       SumD = 0.d0
       do 15 i=1,pseg
         k=(i-1)*npt
         DD(i) = 0.d0

         do j=1,npt
           do l=1,3
             aa(l,j)=r(l,j+k+npt)-r(l,j+k)
             bb(l,j)=r(l,j+k+2*npt)-r(l,j+k+npt)
             DD(i) = DD(i) + (massfac(j)*(r(l,j+k+npt)-r(l,j+k)))**2
           end do
         end do

         normA=0.d0
         normB=0.d0
         a_b = 0.d0

         do j=1,npt
           do l=1,3
              normA=normA + aa(l,j)*aa(l,j)
              normB=normB + bb(l,j)*bb(l,j)
              a_b = a_b + aa(l,j)*bb(l,j)
           end do
         end do

         cosine(i) = a_b / dsqrt(normA*normB)
         write(6,107)i,cosine(i), dsqrt(DD(i)/npt)
         SumD = SumD + dsqrt(DD(i)/npt)
15     continue
107      format(i5,1x,'Cosine',1x,f5.3,1x,'D',1x,f6.4)
       if (last) then
         DD(pseg+1) = 0.d0
         do 27 j=pseg*npt+1,(pseg+1)*npt
           do 27 l=1,3
             DD(pseg+1)=DD(pseg+1)+
     &         (massfac(j-pseg*npt)*(r(l,j+npt)-r(l,j)))**2
27       continue
         write(6,107)igrid-1,3.00,dsqrt(DD(pseg+1)/npt)
         SumD = SumD + dsqrt(DD(pseg+1)/npt)
       endif

       if (paral) then
           if (.not.first) then
              call Send_Double(cosine,pseg,0,procID)
              call Send_Double(DD,pseg,0,1000+procID)
              if (last) call Send_Double(DD(pseg+1),1,0,1001+procID)
           else ! first
             j = pseg + 1

             do n=1,numprocs-1
C               write(6,*) 'Processor ',n
               if (n.lt.mod(igrid-2,numprocs)) then
                 amount = (igrid-2)/numprocs + 1 
               else
                 amount = (igrid-2)/numprocs 
               endif
               call Recv_Double(cosine,amount,n,n)
               call Recv_Double(DD,amount,n,1000+n)
               do i=1,amount
                 write(6,107)j,cosine(i), dsqrt(DD(i)/npt)
                 SumD = SumD + dsqrt(DD(i)/npt)
                 j = j + 1
               enddo
             enddo
             call Recv_Double(DD,1,numprocs-1,1001+numprocs-1)
             write(6,107)j,3.0,dsqrt(DD(1)/npt) 
             SumD = SumD + dsqrt(DD(1)/npt)
           endif ! first
        endif ! paral
        
        write(6,'(a,f10.3)')'Sum of Ds ', SumD


        SumD = 0.d0
        do j=1,npt
          do l=1,3
             aa(l,j)=r_final(l,j)-r_initial(l,j)
             SumD = SumD + aa(l,j)**2
          end do
        end do


        write(6,*)'Processor',procID
        do i=1,pseg+2
          k=(i-1)*npt
          DD(i) = 0.d0

          do j=1,npt
            do l=1,3
              aa(l,j)=r(l,j+k)-r_initial(l,j)
              DD(i) = DD(i) + aa(l,j)*(r_final(l,j)-r_initial(l,j))
            end do
          end do
        end do
108     format(i5,1x,'dproj',1x,f5.3)

        if (first) write(6,108)1, DD(1)/SumD
        do i = 2,pseg+1
           write(6,108)i, DD(i)/SumD
        end do
        if (last) write(6,108)pseg+2, DD(pseg+2)/SumD

        if (paral) then
           if (.not.first) then
              call Send_Double(DD(2),pseg,0,procID)
              if(last) call Send_Double(DD(pseg+2),1,0,1001+procID)
           else ! first
             j = pseg + 2

             do n=1,numprocs-1
C               write(6,*) 'Processor ',n
               if (n.lt.mod(igrid-2,numprocs)) then
                 amount = (igrid-2)/numprocs + 1
               else
                 amount = (igrid-2)/numprocs 
               endif
               call Recv_Double(DD,amount,n,n)
               do i=1,amount
                 write(6,108)j,DD(i)/SumD
                 j = j + 1
               enddo
             enddo
             call Recv_Double(DD,1,numprocs-1,1001+numprocs-1)
             write(6,108)j,DD(1)/SumD
           endif ! first
        endif ! paral

109     format(i5,1x,'|dU/dx|',1x,f10.3)

        write(6,*)'Processor',procID
        SumD = 0.d0

        call UpdateLists(r,npt,pseg)
        do 25 i=1,pseg+2
          k = (i-1)*npt

          do j=1,npt
            do l=1,3
                coor(l,j) = r(l,j+k)*massfac(j)
            end do
          end do

          call eforce()
C          call wener(6)

          DD(i) = 0.d0
          do j=1,npt
            do l=1,3
                DD(i) = DD(i) + dpot(l,j)**2
            end do
          end do
25      continue

        if (first) then
         write(6,109)1,dsqrt(DD(1)/npt3)
         SumD = SumD + dsqrt(DD(1)/npt3)
        endif

        do i=2,pseg+1
         write(6,109)i,dsqrt(DD(i)/npt3)
         SumD = SumD + dsqrt(DD(i)/npt3)
        end do
        if (last) then 
         write(6,109)pseg+2,dsqrt(DD(pseg+2)/npt3)
         SumD = SumD + dsqrt(DD(pseg+2)/npt3)
        endif

        if (paral) then
           if (.not.first) then
              call Send_Double(DD(2),pseg,0,procID)
              if(last) call Send_Double(DD(pseg+2),1,0,1001+procID)
           else ! first
             j = pseg + 2

             do n=1,numprocs-1
C               write(6,*) 'Processor ',n
               if (n.lt.mod(igrid-2,numprocs)) then
                 amount = (igrid-2)/numprocs + 1
               else
                 amount = (igrid-2)/numprocs
               endif
               call Recv_Double(DD,amount,n,n)
               do i=1,amount
                 write(6,109)j, dsqrt(DD(i)/npt3)
                 SumD = SumD + dsqrt(DD(i)/npt3)
                 j = j + 1
               enddo
             enddo
             call Recv_Double(DD,1,numprocs-1,1001+numprocs-1)
             write(6,109)j,dsqrt(DD(1)/npt3)
             SumD = SumD + dsqrt(DD(1)/npt3)
           endif ! first
        endif ! paral

        write(6,'(a,f8.4)')' sqrt(<dU/dx^2>) ',SumD/igrid



C        call MPI_Buffer_detach( mpi_buffer, 30*maxpt,ierr );
        call comm_exit(ierr)

        stop
        end
