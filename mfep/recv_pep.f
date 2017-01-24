      subroutine recv_pep(n,i,outdata,nselec)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'mpif.h'

      double precision outdata(3,nselec),avg_pep(3,nselec)
      integer i,n,ierr
      integer rc
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg


c      CALL MPI_COMM_RANK(MPI_COMM_WORLD,procID,ierr)
c
c      if(procID.eq.n) then
c      CALL MPI_RECV(avg_pep,10, MPI_DOUBLE_PRECISION
c     &,i, 122, MPI_COMM_WORLD, status, rc)
c      CALL MPI_SEND(avg_pep,10, MPI_DOUBLE_PRECISION
c     &,i, 123, MPI_COMM_WORLD, rc)
c      elseif (procID.eq.i) then
c      CALL MPI_SEND(avg_pep,10, MPI_DOUBLE_PRECISION
c     &,n, 122, MPI_COMM_WORLD, rc)
c      CALL MPI_RECV(avg_pep,10, MPI_DOUBLE_PRECISION
c     &,n, 123, MPI_COMM_WORLD, status, rc)
c      endif



c      if(procID.ne.n) then
c        write(stdo,*)
c     >      ' Error on receiving procID=',procID,' n=',n
c        stop 101
c      end if
c      write(stdo,*) 'recv1',avg_pep(1,1),i,n

       if(procID.eq.n) then

c         Call MPI_Recv(outdata, 3*nselec, MPI_DOUBLE_PRECISION,
c     >                        i,1,MPI_COMM_WORLD, status,rc )
c         write(68,*) n,avg_pep(1,1)
 

      endif

       if(procID.eq.i) then
c          Call MPI_Send(avg_pep, 3*nselec, MPI_DOUBLE_PRECISION,
c     >                         n,1,MPI_COMM_WORLD, rc )
c        write(67,*) i,avg_pep(1,1)

       endif
c




      if (rc.ne.0) then
         write(stdo,*)' rc = ',rc
         write(stdo,*)
     >        ' Error on receiving struc from proc i=',i
c         call error message(rc,err_msg)
         write(stdo,*)err_msg(1:80)
         stop 103
      end if
       
        call MPI_FINALIZE(ierr)

      return
      end
