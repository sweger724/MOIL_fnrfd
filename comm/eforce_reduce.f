       subroutine eforce_reduce()
         implicit none
c this routine accumulate all the partial force vectors from
c all of the pe-s and add them up

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'mpif.h'

        integer error, i, j
        character*80 message
        double precision dpotbuf(3,maxpt)

        call MPI_Allreduce(dpot,dpotbuf,3*npt,MPI_DOUBLE_PRECISION,
     1          MPI_SUM,MY_COMM,error)
        do j=1,npt
           do i=1,3
              dpot(i,j)=dpotbuf(i,j)
           end do 
        end do
c
        if (error.ne.0) then
                call error message(error,message)
                write(stdo,*)' Error in eforce_reduce '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if

        return
        end


