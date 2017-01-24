        subroutine TransposeData(Q,Da,Db,Dc,start)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'

        double precision send_buffer
        integer start, Da, Db, Dc
        double precision Q(2,Da,Db,Dc)
        integer rem, per_proc,displace(0:maxpe), size(0:maxpe), block
        integer i,ierror

        per_proc = (Dc-1) / num_pes
        rem = mod(Dc-1,num_pes)
        block = 2* Da*Db
        !write(6,*)"DDD:",Dc,per_proc,rem
        !write(6,*)"DDD:",Da,Db
        
        displace(0) = 0
        do i = 0, num_pes-1
          size(i) = per_proc * block
          if (i.le. rem-1) size(i) = size(i) + block
          displace(i+1) = displace(i) + size(i)
        end do

         !write(6,*)"TTT:",size(my_pe),displace(1),block,per_proc
         call MPI_allgatherV(Q(1,1,1,start),size(my_pe),
     1    MPI_DOUBLE_PRECISION,Q(1,1,1,1),size,displace,
     2    MPI_DOUBLE_PRECISION,MY_COMM,ierror)
       
        return
        end
