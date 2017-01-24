      subroutine excheps_chain()
c
c (last modified by V. Zaloj on Dec. 03, 1999)
c
c  this subroutine updates the forces of chain !!!
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c
      include 'mpif.h'
c

c
c local
      integer i,j
      integer istru
c
c# terra related variables
      integer rc,itemsize
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg
      double precision buffer(3*maxpt)
c
c
c      write(stdo,'(10x,a)') ':==> Start updating chain forces ...'
c
      itemsize = 8
      istru  = pseg + 4
c
c
c.....receive from right: 
c.....last receive from first, which is meaningless here
      if (last) then

         Call MPI_Send(donsger(npt3+1), npt3, MPI_DOUBLE_PRECISION,
     1             (my_pe-1),my_pe, 
     2             MPI_COMM_WORLD, rc )

      else 
         if(first) then
            Call MPI_Recv(buffer, npt3, MPI_DOUBLE_PRECISION,
     1             (my_pe+1),MPI_ANY_TAG, 
     2             MPI_COMM_WORLD, status, rc )
         else
            Call MPI_Sendrecv(donsger(npt3+1),npt3,MPI_DOUBLE_PRECISION,
     1             (my_pe-1),my_pe, 
     2             buffer, npt3, MPI_DOUBLE_PRECISION,
     3             (my_pe+1),MPI_ANY_TAG, 
     4             MPI_COMM_WORLD, status, rc )
         endif
c
         if (rc.ne.0) then
            write(*,*)' rc = ',rc
            write(*,*)' In tcomm lib shift left (even) '
            call error message(rc,err_msg)
            write(*,*)err_msg(1:80)
         end if
c
         j = ndegf+npt3
c
         do 2 i=1,npt3
            donsger(j+i) = donsger(j+i)+buffer(i)
 2       continue 
c
      end if
c     
c
c.....send to right:
      i = ndegf+2*npt3
c     
c.....receive from left
c.....first receive from last, which is meaningless here
      if (first) then

             Call MPI_Send(donsger(i+1), npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe+1),my_pe, 
     2              MPI_COMM_WORLD, rc)
      else
         if(last) then
             Call MPI_Recv(buffer, npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe-1),MPI_ANY_TAG, 
     2               MPI_COMM_WORLD, status, rc)
         else
             Call MPI_Sendrecv(donsger(i+1), npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe+1),my_pe, 
     2              buffer, npt3, MPI_DOUBLE_PRECISION,
     3              (my_pe-1),MPI_ANY_TAG,  
     4              MPI_COMM_WORLD, status, rc)
         endif
         if (rc.ne.0) then
           write(*,*)' rc = ',rc
           write(*,*)' In tcomm lib shift right (even) '
           call error message(rc,err_msg)
           write(*,*)err_msg(1:80)
         end if

         j = 2*npt3
c
         do 3 i=1,npt3
            donsger(j+i) = donsger(j+i)+buffer(i)
 3       continue 
c
      end if
c
c      write(stdo,'(10x,a)') ':==>   End updating chain forces !!!'
c
      return
      end











