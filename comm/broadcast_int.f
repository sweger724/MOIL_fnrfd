        subroutine broadcast_i(conv,n)
          implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'

        integer conv,n
        integer ierror

        call MPI_bcast(conv,n,MPI_INTEGER,0,MY_COMM,ierror)
c@      write(*,*) ' ierror = ',ierror
c@      write(*,*) ' After mpi_bcast conv = ',conv
        return
        end
