
      subroutine Send_Double(value,size,to,mes_type)
        implicit none
        include 'mpif.h'

        integer size,to, mes_type
        double precision value(*)

        integer rc

        call MPI_send(value,size,MPI_DOUBLE_PRECISION,
     &                 to,mes_type,MPI_COMM_WORLD,rc)

        end
