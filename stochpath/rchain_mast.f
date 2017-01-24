      subroutine rchain(urcrd,istart,kstru,ratio)
c
c.v1.0 (last changed 24/2/97)
c
c a subroutine to read initial chain coordinates. Three styles
c are supported: 
c (i)   DYNAmics (binary single precision format for compatability
c       with QUANTA)
c (ii)  PATH (double precision format recommend for computations.
c (iii) INIT (initialize a chain using linear interpolation, urcrd
c              contains the names of normal coordinate files for
c              reactants and products).
c (iv)  INTR (interpolation. Missing structures are interpolated
c             (linearly) between existing structures. 
c (v)   TETA (initialize a chain using a theta-function interpolation,
c              first half of the grid at reactants, second half at
c              products. 
c              urcrd contains the names of normal coordinate files for
c              reactants and products).
c (vi)  EXTR (extraction. The amount of structures is reduced by 
c             taking a less dense grid of structures with some skip 
c             factor.
c
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CCRD.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'

      integer urcrd,istart,kstru,ratio
      integer rbin
      logical find
      integer of
c
c local
      integer i,ii,j,istart_file,istart_eff,iend_eff,kstru_eff,namel
      integer istart_sparse,istart_dense,iend_dense,iend_sparse,i_sparse
      integer i_dense, i_dense0,index_file
      integer u1,u2,j1,j2,j3,k,kk,i1,k1,n
      integer nofreez(maxpt)
      double precision div,coor_temp(3,maxpt),fac1
c
      character*6 name
c
      norew = .false.
      lpstr = 1
c
      rbin = 1
c
      name  = 'rchain'
      namel = 6
c
      npt3 = 3*npt
      npt2 = 2*npt
c
      if (crdstyl.eq.'DYNA') then
c........This part looks OK, but it was NEVER CHECKED
         istart_eff=istart
         kstru_eff=kstru
         istart_file=1
         if (bc1) then
c...........in the path file the structures at the extrema are missing 
c...........since they are not used in the program.
            istart_file=2
            if (first.and.last) then
               istart_eff=istart+1
               kstru_eff=kstru-2
            else if (first) then
               istart_eff=istart+1
               kstru_eff=kstru-1
            else if (last) then
               kstru_eff=kstru-1
            end if
         end if
c
         i1 = npt
         do i=istart_eff,istart_eff+kstru_eff-1
            rewind urcrd
            j1 = i-istart_file+1
            j= (i-istart)*npt3
            call rdyncrd(urcrd,j1,i1,nofreez,rbin)
            do k=0,npt3-3,3
               kk = k/3 + 1
               r(j+k+1)  = coor(1,kk)
               r(j+k+2)  = coor(2,kk)
               r(j+k+3)  = coor(3,kk)
            end do
         end do
c     
      else if (crdstyl.eq.'PATH') then
c
       if(first) then
         rewind urcrd
         if(bc1) then
            istart_eff=2
         else
            istart_eff=1
         endif
         iend_eff=pseg+2
c
c reading on master processor  and storing on master processor !!!
c
         do i=istart_eff,iend_eff
            write(stdo,*)' Reading path format i = ',i
            call rpath_seq(urcrd,i)
            write(stdo,*)' Storing path format i = ',i
            j = (i-1)*npt3
            do k=0,npt3-3,3
               kk = k/3 +1
               r(j+k+1)  = coor(1,kk)
               r(j+k+2)  = coor(2,kk)
               r(j+k+3)  = coor(3,kk)
            end do
         end do
c
c reading on master processor  and sending to master processor !!!
c
       
        if(.not.last) then
          do n=2,n_pes
           istart_eff=iend_eff+1
           iend_eff=iend_eff+pseg
           do i=istart_eff,iend_eff
             write(stdo,*)' Reading path format i = ',i
             call rpath_seq(urcrd,i)
             write(stdo,*)' Sending path format i = ',i,
     >                    ' to PROC_ID = ',(n-1)
             call send_struc(i,(n-1))
           end do
          end do
        endif

         istart_eff=iend_eff+1
         if(bc1) then
            iend_eff=istart_eff
         else
            iend_eff=istart_eff+1
         endif

         do i=istart_eff,iend_eff
           write(stdo,*)' Reading path format i = ',i
           call rpath_seq(urcrd,i)
           if(last) then
            write(stdo,*)' Storing path format i = ',i
            j = (i-1)*npt3
            do k=0,npt3-3,3
               kk = k/3 +1
               r(j+k+1)  = coor(1,kk)
               r(j+k+2)  = coor(2,kk)
               r(j+k+3)  = coor(3,kk)
            end do
           else
             write(stdo,*)' Sending path format i = ',i,
     >                    ' to PROC_ID = ',(n_pes-1)
             call send_struc(i,(n_pes-1))
           endif
         end do
c
       else 
c
c    Receiving section (only for paralel code !!!)
c
         istart_eff=my_pe*pseg+3
         iend_eff=istart_eff+pseg-1
         if(last) then
            if(bc1) then
               iend_eff=iend_eff+1
            else
               iend_eff=iend_eff+2
            endif
         endif  
c
         n=my_pe
c
         do i=istart_eff,iend_eff
            write(stdo,*)' Receiving path format i = ',i,
     >                    ' at PROC_ID = ',n
            call recv_struc(i,n)
            j = (i-istart_eff+2)*npt3
            do k=0,npt3-3,3
               kk = k/3 +1
               r(j+k+1)  = coor(1,kk)
               r(j+k+2)  = coor(2,kk)
               r(j+k+3)  = coor(3,kk)
            end do
         end do
c
c     Closing the PATH reading and storing section !!!!
c
       endif
c
      else if (crdstyl.eq.'INIT') then
c
         if (urcrd.ne.stdi) rewind urcrd
c
c........getting unit number of reactants (CHARM format)
         call rline(name,namel,urcrd)
         if (find('file')) then
            u1 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
         end if
c........getting unit number of products (CHARM format)
         call rline(name,namel,urcrd)
         if (find('file')) then
            u2 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
         end if
c
c........reading products
         call getcrd(u2,'CHARM')
c
c........set the position in the path segment of the coordinates of the 
c........last position considered (j). In each path segment only the 
c........positions associated with the segment are kept.
         if (bc1) then
            if (first.and.last) then
               j = (grid-2)*npt3
            else if (last) then
               j = (pseg+2)*npt3
            else
               j = (pseg+3)*npt3
            end if
         else
            if (first.and.last) then
               j = (grid-1)*npt3
            else
               j = (pseg+3)*npt3
            end if
         end if
         do i=0,npt3-3,3
            ii = i/3 + 1
            r(j+i+1)   = coor(1,ii)
            r(j+i+2)   = coor(2,ii)
            r(j+i+3)   = coor(3,ii)
         end do
c
c........set the position in the path segment of the coordinates of the 
c........first position considered (starting from zero).
c........reading reactants
         call getcrd(u1,'CHARM')
         if (bc1) then
            if (first) then
               j1 = npt3
            else
               j1 = 0
            end if
         else
            j1=0
         end if
c
         do i=0,npt3-3,3
            ii = i/3 + 1
            r(j1+i+1)   = coor(1,ii)
            r(j1+i+2)   = coor(2,ii)
            r(j1+i+3)   = coor(3,ii)
         end do
c
c........uses the whole grid to get the linear interpolation coeficient
         if (bc1) then
            div = 1.d0/dble(grid-3)
         else
            div = 1.d0/dble(grid-1)
         end if
         do 5 i=1,npt3
            donsger(i) = (r(j+i)-r(j1+i))*div
5        continue
c
c........Use the linear interpolation coeficient and the actual 
c........position of the segment in the path to calculate the initial
c........structures for the "local processor".
c........i=index in the "time" position
c........k=index of the first coordinate at "time" i locally.
c........j and ii = index of the coordinate.
c
         j1=0
         j2=0
         j3=0 
         if (bc1) then
c...........one deals with the  system without the first two 
c...........structures, which are not used in the program.
            j3 = 1
            if (first.and.last) then
               j1 = 1
               j2 = 1
            else if (first) then
               j1 = 1
            else if (last) then
               j2 = 1
            end if
         end if
         do i=istart+j1,istart+kstru-1-j2
            k = (i-istart)*npt3
            do j=0,npt3-3,3
               ii = j/3 +1
               r(j+k+1) = coor(1,ii) + (i-1-j3)*donsger(j+1)
               r(j+k+2) = coor(2,ii) + (i-1-j3)*donsger(j+2)
               r(j+k+3) = coor(3,ii) + (i-1-j3)*donsger(j+3)
            end do
         end do
c
         close(u1)
         close(u2)
c
      else if (crdstyl.eq.'TETA') then
c.......In this case the form of the boundary condition do not 
c.......affect the form of the initial path (since the extrema are 
c.......not used).
c
         if (urcrd.ne.stdi) rewind urcrd
c
c........getting unit number of reactants (CHARM format)
         call rline(name,namel,urcrd)
         if (find('file')) then
            u1 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
         end if
c........getting unit number of products (CHARM format)
         call rline(name,namel,urcrd)
         if (find('file')) then
            u2 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
         end if
c
c........reading products
         call getcrd(u2,'CHARM')
c
c........set the position in the path segment of the coordinates of the 
c........last position considered (j). In each path segment only the 
c........positions associated with the segment are kept.
c
         k1= (kstru-1)*npt3
c
         do i=0,npt3-3,3
            ii = i/3 + 1
            r(k1+i+1)   = coor(1,ii)
            r(k1+i+2)   = coor(2,ii)
            r(k1+i+3)   = coor(3,ii)
         end do
c
c........set the position in the path segment of the coordinates of the 
c........first position considered (starting from zero).
c........reading reactants
         call getcrd(u1,'CHARM')
         do i=0,npt3-3,3
            ii = i/3 + 1
            r(i+1)   = coor(1,ii)
            r(i+2)   = coor(2,ii)
            r(i+3)   = coor(3,ii)
         end do
c
         if ((istart+kstru-1).le.(grid/2))then
c
            do i=istart+1,istart+kstru-1
               k = (i-istart)*npt3
               do j=1,npt3
                  r(j+k) = r(j)
               end do
            end do
c
         else if (istart.gt.(grid/2))then
c
            do i=istart,istart+kstru-2
               k = (i-istart)*npt3
               do j=1,npt3
                  r(j+k) = r(j+k1)
               end do
            end do
c
         else
c
c........i=index in the "time" position
c........k=index of the first coordinate at "time" i locally.
c........j=index of the coordinate.
c
            do i=istart,grid/2
               k = (i-istart)*npt3
               do j=1,npt3
                  r(j+k) = r(j)
               end do
            end do
c
            do i=grid/2+1,istart+kstru-2
               k = (i-istart)*npt3
               do j=1,npt3
                  r(j+k) = r(j+k1)
               end do
            end do
c
         end if
c.........debug checks
c         do 7 i=1,kstru
c            write(*,*)' *** structure i ',i
c            k = (i-1)*npt3
c            do 7 j=0,npt3-3,3
c               write(*,*)j,r(j+k+1),r(j+k+2),r(j+k+3)
c7        continue
c
         close(u1)
         close(u2)
c
      else if (crdstyl.eq.'INTR') then
c
c
c........Read the path in the PATH format. Ratio has to be such that
c........bc1: ratio=(grid_dense-3)/(grid_sparse-3)(#inputs grid_dense-2)
c........bc2: ratio=(grid_dense-1)/(grid_sparse-1)(#inputs grid_dense) 
c
c........This part is general, also for case in which the first and 
c........last structure in the sparse grid correspond to structures 
c........outside the range of the parallel segment.
c
         rewind urcrd
         istart_eff=istart
         kstru_eff=kstru
         istart_file=1
         if (bc1) then
c...........one reads (and writes) the system without the two 
c...........structures at the extrema, which are not used in the 
c...........program.
            istart_file=2
            if (first.and.last) then
               istart_eff=istart+1
               kstru_eff=kstru-2
            else if (first) then
               istart_eff=istart+1
               kstru_eff=kstru-1
            else if (last) then
               kstru_eff=kstru-1
            end if
         end if
         istart_sparse=(istart_eff-istart_file)/ratio +istart_file
         istart_dense=(istart_sparse-istart_file)*ratio +istart_file
         iend_dense=istart_eff+kstru_eff-1
         iend_sparse=int(dble(iend_dense-istart_file)/dble(ratio)
     &        +0.9999999d0) +istart_file
c
         div = 1.d0/dble(ratio)
c
         do i_sparse=istart_file,istart_sparse-1
            call rpath_seq(urcrd,-i_sparse)
         end do
c
c........read first coordinate
         i_sparse=istart_sparse
         i_dense0=(i_sparse-istart_file)*ratio +istart_file
         write(*,*)' Reading path format i istart kstru ',
     1         i_dense0,istart,kstru
         call rpath_seq(urcrd,i_sparse)
         do k=1,npt
            do kk=1,3
                coor_temp(kk,k) = coor(kk,k)
             end do
         end do
c
         do i_sparse=istart_sparse+1,iend_sparse-1
            i_dense=(i_sparse-istart_file)*ratio+istart_file
            write(*,*)' Reading path format i istart kstru ',
     1           i_dense,istart,kstru
            call rpath_seq(urcrd,i_sparse)
c
c...........interpolate
            do i=max(i_dense0,istart_eff),min(i_dense,iend_dense)-1
               fac1=dble(i-i_dense0)*div
               j = (i-istart_eff)*npt3
               do k=0,npt3-3,3
                  kk = k/3 +1
                  r(j+k+1)= coor_temp(1,kk)+ 
     &                 fac1*(coor(1,kk)-coor_temp(1,kk))
                  r(j+k+2)= coor_temp(2,kk)+ 
     &                 fac1*(coor(2,kk)-coor_temp(2,kk))
                  r(j+k+3)= coor_temp(3,kk)+ 
     &                 fac1*(coor(3,kk)-coor_temp(3,kk)) 
               end do
            end do
c...........update
            i_dense0=i_dense
            do k=1,npt
               do kk=1,3
                  coor_temp(kk,k) = coor(kk,k)
               end do
            end do
c
         end do
c
         i_sparse=iend_sparse
         i_dense=(i_sparse-istart_file)*ratio+istart_file
         write(*,*)' Reading path format i istart kstru ',
     1        i_dense,istart,kstru
         call rpath_seq(urcrd,i_sparse)
c     
c........interpolate
         do i=max(i_dense0,istart_eff),min(i_dense,iend_dense)
            fac1=dble(i-i_dense0)*div
            j = (i-istart_eff)*npt3
            do k=0,npt3-3,3
               kk = k/3 +1
               r(j+k+1)= coor_temp(1,kk)+ 
     &              fac1*(coor(1,kk)-coor_temp(1,kk))
               r(j+k+2)= coor_temp(2,kk)+ 
     &              fac1*(coor(2,kk)-coor_temp(2,kk))
               r(j+k+3)= coor_temp(3,kk)+ 
     &              fac1*(coor(3,kk)-coor_temp(3,kk)) 
            end do
         end do
c
      else if (crdstyl.eq.'EXTR') then
c
c........Read the path in the PATH format. Ratio has to be such that
c........bc1: ratio=(grid_dense-3)/(grid_sparse-3)(#inputs grid_dense-2)
c........bc2: ratio=(grid_dense-1)/(grid_sparse-1)(#inputs grid_dense) 
c
         rewind urcrd
         istart_eff=istart
         kstru_eff=kstru
         istart_file=1
         if (bc1) then
c...........one reads (and writes) the system without the two 
c...........structures at the extrema, which are not used in the 
c...........program.
            istart_file=2
            if (first.and.last) then
               istart_eff=istart+1
               kstru_eff=kstru-2
            else if (first) then
               istart_eff=istart+1
               kstru_eff=kstru-1
            else if (last) then
               kstru_eff=kstru-1
            end if
         end if
c
c........reads the path in a sequencial way until it reaches the desired
c........structure
         index_file=istart_file
         do i=istart_file,istart_eff-1
            do ii=1,ratio
               call rpath_seq(urcrd,-index_file)
               index_file=index_file+1
            end do
         end do
         do i=istart_eff,istart_eff+kstru_eff-2
            write(*,*)' Reading path format i istart kstru ',i,
     1           istart,kstru
c...........rpath_seq reads only the desired structure (the read here 
c...........is sequential anyway). The rpath reads the whole file 
c...........until the structure i (including))
            call rpath_seq(urcrd,index_file)
            index_file=index_file+1
            j = (i-istart)*npt3
            do k=0,npt3-3,3
               kk = k/3 +1
               r(j+k+1)  = coor(1,kk)
               r(j+k+2)  = coor(2,kk)
               r(j+k+3)  = coor(3,kk)
            end do
            do ii=1,ratio-1
               call rpath_seq(urcrd,-index_file)
               index_file=index_file+1
            end do
         end do
c........last structure
         write(*,*)' Reading path format i istart kstru ',i,
     1        istart,kstru
         call rpath_seq(urcrd,index_file)
         j = (i-istart)*npt3
         do k=0,npt3-3,3
            kk = k/3 +1
            r(j+k+1)  = coor(1,kk)
            r(j+k+2)  = coor(2,kk)
            r(j+k+3)  = coor(3,kk)
         end do
c
       end if
c
       return
c
c@       level = 1
c@       call alert(name,namel,'Error reading nstru data',24,level)
c
c@       return
       end
c
c
c  This part should go to the communications !!!!!
c
      subroutine send_struc(i,n)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
c
      include 'mpif.h'
c
c
      integer i,n
      integer rc
      character*80 err_msg
c
      if(my_pe.ne.0) then
        write(stdo,*)
     >      ' Error: data should be sent from master ! my_pe=',my_pe
        stop 101
      endif
c
      if((n.le.0).or.(n.ge.n_pes)) then
        write(stdo,*)
     >      ' Error: out of target processor ! n=',n
        stop 102
      endif   
c   
c

      Call MPI_Send(coor, npt3, MPI_DOUBLE_PRECISION,
     >                         n,i,MPI_COMM_WORLD, rc )

      if (rc.ne.0) then
         write(stdo,*)' rc = ',rc
         write(stdo,*)
     >        ' Error on sending struc i=',i," to proc ID=",n
         call error message(rc,err_msg)
         write(stdo,*)err_msg(1:80)
         stop 103
      end if
c
      return
      end
c
      subroutine recv_struc(i,n)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
c
      include 'mpif.h'
c
c
      integer i,n
      integer rc
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg
c
      if(my_pe.ne.n) then
        write(stdo,*)
     >      ' Error on receiving my_pe=',my_pe,' n=',n
        stop 101
      endif
c

      Call MPI_Recv(coor, npt3, MPI_DOUBLE_PRECISION,
     >                         0,i,MPI_COMM_WORLD, status,rc )

      if (rc.ne.0) then
         write(stdo,*)' rc = ',rc
         write(stdo,*)
     >        ' Error on sending struc i=',i," to proc ID=",n
         call error message(rc,err_msg)
         write(stdo,*)err_msg(1:80)
         stop 103
      end if
c
      return
      end


