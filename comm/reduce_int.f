        subroutine reduce_int(x,n)
        implicit none
c
c this routine adds from all pe-s a single integer number
c and put the result in pe-s.
c

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'

        integer x(0:*),y(0:maxpe)
        integer n

        integer i
        integer error

        call MPI_Allreduce(x(0),y(0),n, MPI_INTEGER,
     1          MPI_SUM, MY_COMM,error)

        if (error.ne.0) then
          call alert('reduce_int',10,' reduce failed',14,1)
        end if

        do i=0,n
C           write(6,*)"aaa:",i,x(i),y(i)
           x(i)=y(i)
        end do
        return
        end
