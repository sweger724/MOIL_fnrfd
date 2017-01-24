      subroutine comm_exit(ierr)
c
      include 'COMMON/LENGTH.BLOCK'
      integer ierr
c
      include 'mpif.h'
c
c
      Call MPI_FINALIZE(ierr)
c
      ierr=0
c
      return
      end


