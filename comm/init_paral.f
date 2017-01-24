      subroutine init_paral()
      implicit none
c Subroutine for MPI initialization of straightforward MD
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'
c     
      character*10 name
      integer namel
      integer ierr
c
c
      name   = 'init_paral'
      namel  = 10
c     
      call MPI_INIT( ierr )
      call get my pe number()
      write(6,*)' my processing element ',my_pe
      call get number pes (num_pes)

      MY_COMM = MPI_COMM_WORLD
      prll_on_off = .true.
c
      return
      end





