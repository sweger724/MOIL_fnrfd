        subroutine broadcast_state()
        implicit none

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/COORD.BLOCK'
          include 'COMMON/VELOC.BLOCK'
          include 'COMMON/CONNECT.BLOCK'
          include 'COMMON/PARALLEL.BLOCK'
          include 'mpif.h'

          integer ierr

        call MPI_BCAST(coor,3*npt,MPI_DOUBLE_PRECISION,0, MY_COMM, ierr)
        call MPI_BCAST(velo,3*npt,MPI_DOUBLE_PRECISION,0, MY_COMM, ierr)

        end 
