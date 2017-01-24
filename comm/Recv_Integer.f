      subroutine Recv_Integer(value,size,from,mes_type)
        implicit none
        include 'mpif.h'

        integer size,from, mes_type
        integer value(*)


        integer status(MPI_STATUS_SIZE), rc

        call MPI_recv(value,size,MPI_Integer,
     &        from,mes_type,MPI_COMM_WORLD,status,rc)

        end
