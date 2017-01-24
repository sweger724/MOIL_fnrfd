C     
C     send Cdata from processor "from" to processor "to"
C   if store==TRUE the data received are added to the
C   value originaly present at processor "to"
C
      subroutine Transmit(from,to,Cdata,store)

      implicit none

      integer from, to
      double precision Cdata,tmp_data
      logical store

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'

      integer status(MPI_STATUS_SIZE),rc

           if (procID .eq. from) then
             call MPI_send(Cdata,1,MPI_DOUBLE_PRECISION,
     &                     to,to,MPI_COMM_WORLD,rc)
C             write(6,*)"Sending from",from," to",to
           end if

           if (procID .eq. to) then
             call MPI_recv(tmp_data,1,MPI_DOUBLE_PRECISION,
     &                     from,to,MPI_COMM_WORLD,status,rc)
C            write(6,*)"Received message from",from," at",to
             if (store) then
                Cdata = tmp_data
             else
               Cdata = Cdata + tmp_data
             endif

           end if

      return
      end
