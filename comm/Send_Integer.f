
      subroutine Send_Integer(value,size,to,mes_type)
        implicit none
        include 'mpif.h'

        integer size,to, mes_type
        integer value(*)


        integer rc

        call MPI_send(value,size,MPI_Integer,
     &                 to,mes_type,MPI_COMM_WORLD,rc)

        end
