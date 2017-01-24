        subroutine reduce_1(x)
         implicit none
c this routine adds from all pe-s a single double precision number
c and put the result in pe-0.

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'mpif.h'
c
        double precision x,y

        integer error
        character*80 message

        call MPI_Allreduce (x,y,1,MPI_DOUBLE_PRECISION,
     1             MPI_SUM,  MY_COMM, error )
        x=y
c
        if (error.ne.0) then
                call error message(error,message)
                write(stdo,*)' Error in reduce_1 '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if

        return
        end
        


