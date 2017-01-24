      subroutine getmypenumberpath()
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'
       
      integer ierr

      call MPI_COMM_RANK( MPI_COMM_WORLD, procID, ierr )

      return
      end









