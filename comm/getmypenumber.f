      subroutine getmypenumber()
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'
      
      integer ierr
      
      if (proc_max.eq.1)then
         procID = 0
         my_pe = 0
         return
      end if
      
      call MPI_COMM_RANK( MPI_COMM_WORLD, procID, ierr )
      
      my_pe = procID
      
      return
      end









