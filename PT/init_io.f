      subroutine init_io()
      
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/PT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'

      character*80 filenamei, filename1
      character name*7
      integer len, namel,ierr
      
      namel=7
      name = "init_io"
      
      ! read number of replicas from standard input
      if (procID .eq. 0) read(stdi,*,err=7,end=7)Nreplicas
      
      ! let everybody know what the number of replicas is
      call mpi_bcast(Nreplicas,1,MPI_INTEGER,
     &                  0,MY_COMM, MPIerr)

      if (Nreplicas .eq. numprocs) then
         prll_on_off = .false.
      else
         prll_on_off = .true.
      end if

      ! every processor sets its replicaID and my_pe (based on procID)
      if (mod(numprocs,Nreplicas).ne. 0) then
        call alert(name,namel,'Number of processors have to be divisible
     & by number of replicas!',25,1)
      end if
      
      num_pes = numprocs / Nreplicas 
       
      my_pe = mod(procID,num_pes)
      replicaID = ( procID - my_pe ) / num_pes

        !write(6,*)"num_pes,Nreplica,my_pe:",num_pes,Nreplicas,my_pe
        !write(6,*)"procID,:",procID,"replicaID:",replicaID

        ! create new groups of processors
      call MPI_Comm_split(MPI_COMM_WORLD,replicaID,my_pe,MY_COMM,ierr)
      !write(6,*)"My new group:",MY_COMM
     

      ! initialize input/output files

      len=13
      write(filenamei,'(a5,i4.4,a4)') 'dyna_',replicaID,'.inp'
      open (unit=stdi,status='unknown',file=filenamei(1:len))
      rewind stdi
     
      if (my_pe .eq. 0) then
        ! only leading processor of each replica writes output
        len=13
        write(filename1,'(a5,i4.4,a4)') 'dyna_',replicaID,'.out'
        open (unit=stdo,status='unknown',file=filename1(1:len))
      else
        len=13
        write(filename1,'(a5,i4.4,a4)') 'dyma_',procID,'.out'
        open (unit=stdo,status='unknown',file=filename1(1:len))
      end if

c     initialzation of acceptance ratios for parallel tempering 
        sumswap1=0.
        sum_accept1=0
        sumswap2=0.
        sum_accept2=0
 
      return
      
7     continue     
        write(6,*)"You need to put a single number on standard input"
        write(6,*)"that specifies number of replicas run!"
        write(6,*)"This can be done by placing this number"
        write(6,*)"to a temporary file and redirecting the "
        write(6,*)"standard input from this file: "
        write(6,*)" i.e.: dynapt < replicaNumber.file"
        call alert(name,namel,'Error during read',17,1)
        return
      end
