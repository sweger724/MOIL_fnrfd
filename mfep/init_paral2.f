      subroutine init_paral2(my_pe,num_pes)
      implicit none
      include 'mpif.h'
c
c (dec 31 2003 -- RE)
c
c Subroutine for MPI initialization of straightforward MD
c
      include 'COMMON/LENGTH.BLOCK'
c      include 'COMMON/PARALLEL.BLOCK'
c     
      character*10 name
      integer namel,my_pe,num_pes
      integer ierr
c
c
      name   = 'init_paral'
      namel  = 10
c     
      call MPI_INIT( ierr )
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_pes,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,my_pe,ierr)
     

c
      return
      end





