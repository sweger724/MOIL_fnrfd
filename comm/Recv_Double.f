
      subroutine Recv_Double(value,size,from,mes_type)
        implicit none
        include 'mpif.h'

        integer size,from, mes_type
        double precision value(*)


        integer status(MPI_STATUS_SIZE), rc

        call MPI_recv(value,size,MPI_DOUBLE_PRECISION,
     &        from,mes_type,MPI_COMM_WORLD,status,rc)

        if (from .eq. MPI_ANY_SOURCE) then
        
          from = status(MPI_SOURCE)
        end if
C         write(6,*)"message from",from,"with tag",mes_type,"received."

        end
