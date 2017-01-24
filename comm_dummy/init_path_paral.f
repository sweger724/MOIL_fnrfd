      subroutine init_path_paral()
      implicit none
c
c Subroutine for the initialization of the parallel machine and the 
c variables related to it.
c (last modified in 30/5/2008 by Peter)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
c     
      integer ierr
c
      character*10 name
      integer namel
c
c     
      name   = 'init_paral'
      namel  = 10
c
c.....subroutines for processor identification in a parallel machine
      my_pe = 0
      write(*,*)' my processing element ',my_pe
c
      numprocs = 1
      write(*,*)' number of processes', numprocs
c
      if (numprocs.eq.1) then
         first = .true.
         last = .true.
         if (my_pe.ne.0)then
            level = 1
            call alert(name,namel,'wrong number of processors',
     1           26,level)
         end if
      else
         if (my_pe.eq.0) then
            first = .true.
            last = .false.
         else if (my_pe+1.eq.numprocs) then 
            first = .false.
            last = .true.
         else
            first = .false.
            last  = .false.
         end if
      end if
      write(*,*) " # proc=",my_pe," : first ",first,"  last ",last
c
c
c.....redefine istart and istru (defined in init_sto) in case of 
c.....running in parallel.
c
      paral=.not.(first.and.last)

      return
      end
