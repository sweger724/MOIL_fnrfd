      subroutine send_struc(i,n)
        implicit none
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
c
      include 'mpif.h'
c
c
      integer i,n
      integer rc
      character*80 err_msg
c
      if(procID.ne.0) then
        write(stdo,*)
     >      ' Error: data should be sent from master ! procID=',procID
        stop 101
      endif
c
      if((n.le.0).or.(n.ge.numprocs)) then
        write(stdo,*)' Error: out of target processor ! n=',n
        stop 102
      endif
c
c
      Call MPI_Send(coor, npt3, MPI_DOUBLE_PRECISION,
     >                         n,i,MPI_COMM_WORLD, rc )

      if (rc.ne.0) then
         write(stdo,*)' rc = ',rc
         write(stdo,*)' Error on sending struc i=',i,' to proc ID=',n
         call error message(rc,err_msg)
         write(stdo,*)err_msg(1:80)
         stop 103
      end if
c
      return
      end
