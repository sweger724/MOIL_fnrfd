      subroutine recv_struc(i,n)
        implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'mpif.h'


      integer i,n
      integer rc
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg

      if(procID.ne.n) then
        write(stdo,*)
     >      ' Error on receiving procID=',procID,' n=',n
        stop 101
      endif

      Call MPI_Recv(coor, 3*npt, MPI_DOUBLE_PRECISION,
     >                         0,i,MPI_COMM_WORLD, status,rc )

      if (rc.ne.0) then
         write(stdo,*)' rc = ',rc
         write(stdo,*)
     >        ' Error on sending struc i=',i," to proc ID=",n
         call error message(rc,err_msg)
         write(stdo,*)err_msg(1:80)
         stop 103
      end if

      return
      end
