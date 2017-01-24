        subroutine gather_nbmlists(x)
          implicit none
c
c Send the total number of neighbours (for this process)
c

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'mpif.h'
c

        integer x(0:*),y


        integer error
        character*80 message

        error  = 0

        y = x(0)

c       call tcomm lib all gather (x,x,1,4,error)
        call MPI_Allgather(y,1,MPI_INTEGER,x,1,MPI_INTEGER,
     1    MY_COMM,error)

        if (error.gt.0) then
czva            call tcomm lib error message(error,message)
                call error message(error,message)
                write(stdo,*)' Error in gather_nbmlists '
                write(stdo,*)' error = ',error
                write(stdo,1)message
        end if
1       format(80a)

        return
        end


