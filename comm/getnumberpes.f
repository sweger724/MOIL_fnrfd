      subroutine get number pes()

      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'

      integer ierr
      
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr )
        num_pes = numprocs

      return
      end


