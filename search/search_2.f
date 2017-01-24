      subroutine search_2(temper, constant_tmp,
     1     pointr,nselec,d0,grdcmx,grdcmy,grdcmz,
     2     udata,ucrd,
     3     debug,scalar,sigmav,cgrd2,
     5     dtopt,clo,npri)    
c    
      implicit none
      double precision tmp1, tmp2, tmp3
      double precision temper, cgrd2
      integer nselec,ncalls
      integer ucrd,udata,finished,npri
      logical debug, ireq
c     
c     npri   -  print some useful(?) data each NPRI steps
c     nselec -  number of selected particles, on which the chain constraints
c     are imposed.
c     ucrd   -  unit number of file on which the coordinates (path format)
c     are written.
c     udata  -  write info on the run on unit UDATA
c     nwcrd  -  write coordinates on ucrd each NWCRD steps
c     debug  - if .true. print a LOT of debugging info
c     
      double precision grdcmx(*),grdcmy(*),grdcmz(*)
      double precision d0(*)
      double precision scalar(*),sigmav(*)
      double precision estred

      integer n, ji, jf, l, xxx
c     
c     dmass - double precision mass vector (m(i=1,npt)(j),j=1,3*igrid)
c     divms - double precision 1/mass vector(1/m(i=1,npt)(j),j=1,3*igrid)
c     
c     d0         -  vector of distances between i,i+1 & i,i+2 pairs
c     grdcm[x-z] -  gradient of center of mass constraints
c     r - coordinates (rx(i=1,npt),ry(i=1,npt),rz(i=1,npt)(j=1,igrid))
c     dsall- forces      (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c     
      integer pointr(*)
c     
c     pointr - a pointer to the selected particles
c     (which are subject to chain const.
c     


C     
C     next common blocks are used in the energy calculation.
C     list, debugging, coordinate transfer.. etc
C     
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/PATH2.BLOCK'
      include 'COMMON/SEARCH.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'
      
c     
c     local
      integer FindNeighbour
      integer i,j,k,ndegf,imin
      integer ncall, initSearch
      double precision sener,dtopt
      logical constant_tmp
      double precision step,factor,ss(2),tsa1(2*lgrid-2)
      double precision clo,dclo,run_temp
      
      integer closed, Csize, expandNode, count, oldSize,from
      integer missed
      logical found, done
      double precision startS(3,maxpt), Pp, dU, minRms

      ncall  = 0
      npt3  = npt*3
      
      dclo=clo/(nstep+1)

      ireq = 0


C  compute energy of the first stucture
      if (first) then
        do i =1, npt
          do l =1,3
             coor(l,i)=r_initial(l,i)
          end do
        end do

        if (eenmyes) call enm_lists(r_initial(1,1),r_final(1,1))
        if (eCGyes)  call CGinit()
        
        if (esymyes) call squeeze()
        call nbondm()
        if (esymyes) call syminit()
        call eforce()
        call wener(stdo)
        
      endif


C  compute energy of the last stucture
      if (first) then
        do i =1, npt
          do l =1,3
             coor(l,i)=r_final(l,i)
          end do
        end do
        if (esymyes) call squeeze()
        call nbondm()
        if (esymyes) call syminit()
        call eforce()
        call wener(stdo)

      endif        

      allmass = 0.d0
         do i = 1,npt
            allmass = allmass + ptms(i)
         end do

      if (paral) then
C  initial and final structures must be comunicated to each processor
        if(first) then
         do i=1,npt
           coor(1,i) = r_initial(1,i)
           coor(2,i) = r_initial(2,i)
           coor(3,i) = r_initial(3,i)
         enddo

         do i=2,numprocs
           call send_struc(1,(i-1))
         end do
        else
          call recv_struc(1,my_pe)
          do i=1,npt
            r_initial(1,i) = coor(1,i)
            r_initial(2,i) = coor(2,i)
            r_initial(3,i) = coor(3,i)
          enddo

        endif

        if(first) then
         do i=1,npt
           coor(1,i) = r_final(1,i)
           coor(2,i) = r_final(2,i)
           coor(3,i) = r_final(3,i)
         enddo

         do i=2,numprocs
           call send_struc(igrid,(i-1))
         end do
        else
          call recv_struc(igrid,my_pe)
          do i=1,npt
            r_final(1,i) = coor(1,i)
            r_final(2,i) = coor(2,i)
            r_final(3,i) = coor(3,i)
          enddo

! initiate CG force fields according to the initial structure
          if (eenmyes) call enm_lists(r_initial(1,1),r_final(1,1))
          if (eCGyes)  call CGinit()
          
        endif

      endif

C     Start A* algorithm

       closed = 0
       
       if (crdstyl .eq. "INIT") Csize = 2
       if (crdstyl .eq. "PATH") Csize = igrid
       oldSize = 0
       
       if (first) then
         if (crdstyl .eq. "INIT") then
           do i=1,npt
             do l =1,3
               centers(l,i)=r_initial(l,i)
               centers(l,i+npt)=r_final(l,i)
             end do
           end do

           call Distance(centers(1,1),r_final(1,1),rms)
           if (rms .lt. cutoff) rms = cutoff
           
           h(1) = - 50.d0 * (rms-cutoff)/cutoff * dlog(0.15d0);
           g(1) = 0.d0
           f(1) = h(1) + g(1)
           ord(1) = 1
           nodeType(1) = 1

           h(2) = 100000.0
           g(2) = 0.d0
           f(2) = h(2) + g(2)
           ord(2) = 2
           nodeType(2) = 2

         end if
         
         if (crdstyl .eq. "PATH") then
           do i = 1 , Csize
             g(i) = 0.d0
             h(i) = 0.d0
             f(i) = h(i) + g(i)
             ord(i) = i
             nodeType(i) = 1
           end do
         end if
       else

         call init_dyna()

       end if
       initSearch = 998

       do while ( closed .lt. Csize .and. Csize .lt. maxNodes-20 )

         done = .false.
         finished = 0

         if (first) then
           do i = 1,Csize
             if (ord(i).eq.closed+1) then
                expandNode = i
             end if
           end do

           do i = oldSize+1, Csize
             call Write2File(ucrd,npt,centers(1,(i-1)*npt+1),0.0,0,0)
           end do
           
           oldSize = Csize
           count = 0
           missed = 0

           do i=1,npt
             do l =1,3
               startS(l,i) = centers(l,(expandNode-1)*npt+i)
             end do
           end do

           if (crdstyl .eq. "PATH") then
             minRms = 1000.d0
             do i=1, igrid
               call Distance(startS,centers(1,(i-1)*npt+1),rms)
               if (rms .lt. minRms) minRms = rms
             end do
             if (minRms .gt. searchLimit) then
               done = .true.
               write(6,*) "Structure", expandNode, "is too far"
               write(6,*) "from the initial path and"
               write(6,*) "is not going to be expanded."
             end if
           end if

           if (done) startS(1,1)=0.d0

           do i = 1,numprocs-1
C              write(stdo,*)"Sending expand node to process",i
              call Send_Double(startS(1,1),3*npt,i,initSearch)
           end do
         else
C              write(stdo,*)"Receiving expand node at proces",my_pe
              call Recv_Double(startS(1,1),3*npt,0,initSearch)
              if (startS(1,1).eq.0.d0) then
                done = .true.
              else
                call ReceiveFinishSignal(ireq)
              end if
         end if

         if (crdstyl .eq. "INIT") then 
C           call CheckFinal(startS(1,1),r_final(1,1),found)
           if (found) then
             write(stdo,*)"Final structure reached!"
           if (first) then
C             write(stdo,*)"Optimal path:"
C             write(stdo,*)" ",expandNode
C             do while (expandNode.ne.1)
C               do i=1,expandNode-1
C                 k = FindNeighbour(i,expandNode)
C                 if (k .ne. -1) then
C                 if (abs(g(i)-dlog(P(i,k))-g(expandNode)).lt.1.d-6) then
C                     expandNode =i
C                     write(6,*)" ",expandNode
C                     goto 100
C                   end if
                 
C                 endif
C               end do
               
C               write(stdo,*) "ERROR!!!"
C               stop
C100            continue             
C             end do

           end if ! first
             return
           end if
         end if ! crdstyl = INIT

         do while (.not. done)
           
           if (.not. first) then

             do i=1,npt
               do l =1,3
                 coor(l,i) = startS(l,i)
               end do
C             write(6,*)"CCC:",coor(1,i),coor(2,i),coor(3,i)
             end do

C             write(stdo,*)"Starting brownian simulation."
             call Brownian(npri,ireq,done)
           else
             from = MPI_ANY_SOURCE
            call Recv_Double(structures(1,count*npt+1),3*npt+1,from,100)
             xxx = int(structures(1,(count+1)*npt+1))
C             write(stdo,*)"Received result."
             if (structures(1,count*npt+1) .ne. 0.d0 
     &          .and. (count+missed) .lt. Stotal) then
               count = count + 1
               write(stdo,*)"Count:",count, Csize, from, xxx
               call AddCenter(count,Csize,expandNode)
               if ( (count+missed) .eq. Stotal) then
                 do i=1, numprocs-1
                 ! let i-th processor know that we 
                 ! collected all Stotal structures
                   call Send_Integer(3,1,i,155)
                 end do
               end if
             else if (structures(1,count*npt+1) .eq. 0.d0) then
               if (structures(1,count*npt+2) .eq. 0.d0) then
                  finished = finished +1
               else  
                  missed = missed + 1
                  write(stdo,*)"Missed:",missed, Csize, from
                  if ( (count+missed) .eq. Stotal) then
                    do i=1, numprocs-1
                      call Send_Integer(3,1,i,155)
                    end do
                  end if
                  
               end if
             end if ! first
             
             if ((count+missed) .ge. Stotal .and. 
     &           finished .eq. numprocs-1) done = .true.
           end if

         end do

         if (first) then
           write(stdo,*)"Clustering..."
           call ClusterStructures(Csize,oldSize,count,expandNode)
           write(stdo,*)"evaluating..."
           call EvaluateNewCenters(Csize,oldSize,closed,expandNode)
           write(stdo,*)"sorting..."
           call SortCenters(Csize,closed)
           write(stdo,*)"DONE Clustering."
        
           do i=1,numprocs-1
             call Send_Integer(Csize,1,i,1)
           end do
         else
C           write(6,*)"Csize is:",Csize
           call Recv_Integer(Csize,1,0,1)
C           write(6,*)"Csize is:",Csize
         endif
         closed = closed + 1
         if (first) call WriteState(Csize,closed)
       end do

       write(stdo,*)"Iterations finished", Csize, closed

       if (first) then
         do i = oldSize+1, Csize
           call Write2File(ucrd,npt,centers(1,(i-1)*npt+1),0.0,0,0)
         end do
      end if

      return
      end

      
      subroutine CheckFinal(coor3,coor4,found)
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        
           double precision coor3(3,maxpt), coor4(3,maxpt)
           double precision rms
           logical found
        
           call Distance(coor3(1,1),coor4(1,1),rms)
           found = .false.
           if (rms .lt. cutoff) found = .true. 
        return
      end
