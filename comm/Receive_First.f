        subroutine Receive_First(Csize,Cdata)

        implicit none

        integer Csize
        double precision Cdata(*)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'

        integer status(MPI_STATUS_SIZE),rc

C      receive from my_pe-1
C      send to my_pe-1

        if (.not.first) then
           call MPI_recv(Cdata(1),Csize,MPI_DOUBLE_PRECISION,
     &                   procID-1,1,MPI_COMM_WORLD,status,rc)
           call MPI_send(Cdata(1+Csize),Csize,MPI_DOUBLE_PRECISION,
     &                  procID-1,2,MPI_COMM_WORLD,rc)
        endif

        return
        end
