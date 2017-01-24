      subroutine init_paral(istart,istru)
c
c.v1.0 (last changed 24/2/97)
c
c Subroutine for the initialization of the parallel machine and the 
c variables related to it.
c (last modified in 30/4/96 by RO)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/PVM_LOCAL.BLOCK'
      include 'COMMON/COMM_LOCAL.BLOCK'
c     
      integer istart,istru,ierr
c
      character*10 name
      integer namel
c
      include 'mpif.h'
c
c     
      call MPI_INIT( ierr )
c
      name   = 'init_paral'
      namel  = 10
c
c.....subroutines for processor identification in a parallel machine
      pvm_first_call=.true.
      comm_first_call=.true.
      call get my pe number(my_pe)
      write(*,*)' my processing element ',my_pe
      proc_id = my_pe + 1
c
      call get number pes(n_pes)
      write(*,*)' number of processes', n_pes
c
c.....first and last are defined so that one can by pass the 
c.....definition of proc, i.e, if nothing declared the default is 
c.....to run with a single processor. (modified for using pvm)
c.....some consistency checks:
c
      if (n_pes.ne.proc_max) then 
c........this cannot happen in pvm since it starts exactly proc 
c........processes, or kill everybody if fails. But in MPI or in general
c........it is possible.
         level = 1
         call alert(name,namel,'wrong number of processors',
     1        26,level)
      end if
c
      if (n_pes.eq.1) then
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
         else if (my_pe+1.eq.n_pes) then 
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
c
      if (paral) then
         istart = my_pe*pseg + 1
         istru  = pseg + 4
      else
         istart  = 1
         istru   = grid
      end if

c
ct......calls related to the talk() subroutine for comunication
ct......between workstations (tcpip protocol)  
ct
ct      proj_id = 1
ct      max_proc =  1
ct      proc_id = 1
c
ct      lme = 6
ct      lup = 6
ct      low = 6
ct      me(1:lme)    = 'noname'
ct      upper(1:lup) = 'noname'
ct      lower(1:low) = 'noname'
c      
      return
      end





