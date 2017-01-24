        subroutine reduce_n(x,n)
c
c this routine adds from all pe-s a single double precision number
c and put the result in pe-0.
c

c# terra include file
        include 'COMMON/UNITS.BLOCK'

czva        include '../include/Tcommlib_f.h'
czva    include '../include/Trtlib_f.h'
c
      include 'mpif.h'
c


        double precision x(*),y(10)
        integer n, i

        integer root,count
        integer error
        character*80 message

        root = 0
        count = n

C@
C       write(stdo,*)' root x count TCOMMLIB_OP_SUM TCOMMLIB_D_DOUBLE ' 
C       write(stdo,*)root,x,count,TCOMMLIB_OP_SUM,TCOMMLIB_D_DOUBLE
 
czva        call tcomm lib all reduce (x,x,count,TCOMMLIB_OP_SUM,
czva     1              TCOMMLIB_D_DOUBLE,error)
c
czva            write(stdo,*) " All reduce from : reduce_n"

            call MPI_Allreduce (x,y,count,MPI_DOUBLE_PRECISION, 
     1                MPI_SUM, MPI_COMM_WORLD, error )
        do 333 i=1,n
         x(i)=y(i)
 333    continue
c
        if (error.ne.0) then
czva                call tcomm lib error message(error,message)
                call error message(error,message)
                write(stdo,*)' Error in reduce_n '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if

        return
        end
        


