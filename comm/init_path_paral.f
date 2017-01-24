      subroutine init_path_paral()
      implicit none
c
c Subroutine for the initialization of the parallel machine and the 
c variables related to it.
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'mpif.h'
     
      integer ierr
      character*10 name
      integer namel


      call MPI_INIT( ierr )

      name   = 'init_paral'
      namel  = 10
c
c.....subroutines for processor identification in a parallel machine
      call get my pe number path()
      write(6,*)' my processing element ',procID
c
      call get number pes()
      write(6,*)' number of processes', numprocs
c
      if (numprocs .eq. 1) then
         first = .true.
         last = .true.
         if (procID .ne. 0)then
            level = 1
            call alert(name,namel,'wrong number of processors',
     1           26,level)
         end if
      else
         if (procID.eq.0) then
            first = .true.
            last = .false.
         else if (procID+1.eq.numprocs) then 
            first = .false.
            last = .true.
         else
            first = .false.
            last  = .false.
         end if
      end if
      write(*,*) " # proc=",procID," : first ",first,"  last ",last
c
c
c.....redefine istart and istru (defined in init_sto) in case of 
c.....running in parallel.
c
      paral=.not.(first.and.last)

      return
      end
