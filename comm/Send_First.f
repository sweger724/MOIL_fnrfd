        subroutine Send_First(Csize,Cdata)

        implicit none

        integer Csize
        double precision Cdata(*)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'

        integer status(MPI_STATUS_SIZE),rc

C      send to my_pe+1
C      receive from my_pe+1

        if (.not.last) then
         call MPI_send(Cdata(pseg*Csize+1),Csize,MPI_DOUBLE_PRECISION,
     &                 procID+1,1,MPI_COMM_WORLD,rc)

        call MPI_recv(Cdata((pseg+1)*Csize+1),Csize,
     &       MPI_DOUBLE_PRECISION,procID+1,2,MPI_COMM_WORLD,status,rc)
        endif

      return
      end
