        subroutine gather_int(conv,all)
          implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'
        integer all(maxpe)
        integer conv
        integer ierror
        integer i

        call MPI_GATHER(conv,1,MPI_INTEGER,all(2),1,
     1                  MPI_INTEGER,0,MY_COMM,ierror)

        return
        end
