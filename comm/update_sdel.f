      subroutine update_chain(upcpucrd)
c
c (last modified by V. Zaloj on Dec. 03, 1999)
c
c  this subroutine updates the coordinates between the processors
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
      integer npt6,istru
c
c# terra related variables
      integer rc,itemsize
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg
      double precision buffer(3*maxpt)
c
      real*8 upcpucrd,elapsedtime
c
      upcpucrd=elapsedtime()
c
c      write(stdo,'(10x,a)') ':==> Start updating chain ...'
c
      npt6  = 6*npt
      itemsize = 8
      istru  = pseg + 4
c
c.....receive from right: 
c.....last receive from first, which is meaningless here
      if (last) then

         Call MPI_Send(r(npt6+1), npt3, MPI_DOUBLE_PRECISION,
     1             (my_pe-1),my_pe, 
     2             MPI_COMM_WORLD, rc )

      else 
         if(first) then
            Call MPI_Recv(buffer, npt3, MPI_DOUBLE_PRECISION,
     1             (my_pe+1),MPI_ANY_TAG, 
     2             MPI_COMM_WORLD, status, rc )
         else
            Call MPI_Sendrecv( r(npt6+1), npt3, MPI_DOUBLE_PRECISION,
     1             (my_pe-1),my_pe, 
     2             buffer, npt6, MPI_DOUBLE_PRECISION,
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
         j = npt6
         do 2 i=1,npt3
            r(ndegf+j+i) = buffer(i)
 2       continue
      end if
c     
c
c.....send to right:
      i = ndegf+npt3
c     
c.....receive from left
c.....first receive from last, which is meaningless here
      if (first) then

             Call MPI_Send( r(i+1), npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe+1),my_pe, 
     2              MPI_COMM_WORLD, rc)
      else
         if(last) then
             Call MPI_Recv(buffer, npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe-1),MPI_ANY_TAG, 
     2               MPI_COMM_WORLD, status, rc)
         else
             Call MPI_Sendrecv( r(i+1), npt3, MPI_DOUBLE_PRECISION,
     1              (my_pe+1),my_pe, 
     2              buffer, npt6, MPI_DOUBLE_PRECISION,
     3              (my_pe-1),MPI_ANY_TAG, 
     4              MPI_COMM_WORLD, status, rc)
         endif
         if (rc.ne.0) then
           write(*,*)' rc = ',rc
           write(*,*)' In tcomm lib shift right (even) '
           call error message(rc,err_msg)
           write(*,*)err_msg(1:80)
         end if
c
         do 3 i=1,npt3
            r(npt3+i) = buffer(i)
 3       continue
      end if
c     
c      write(stdo,'(10x,a)') ':==>   End updating chain !!!'
c
      upcpucrd=elapsedtime()-upcpucrd
c
      return
      end











