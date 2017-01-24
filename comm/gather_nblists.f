        subroutine gather_nblists(x)
          implicit none
c
c Send the total number of neighbours (for this process)
c


        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'
c

        integer x(*), y(4), i


        integer error
        character*80 message

        error  = 0

        do i = 1, 4
          y(i) = x(i)
        end do

        call MPI_Allgather(y,4,MPI_INTEGER,x,4,MPI_INTEGER,
     1    MY_COMM,error)

        if (error.gt.0) then
                call error message(error,message)
                write(stdo,*)' Error in gather_nblists '
                write(stdo,*)' error = ',error
                write(stdo,1)message
        end if
1       format(80a)

        return
        end


