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
      Nreplicas = 1 

       ! read the replica ID from standard input
       if (procID .eq. 0) read(stdi,*,err=7,end=7)replicaID

       ! let everybody know what the number of our replica
       call mpi_bcast(replicaID,1,MPI_INTEGER,
     &                  0,MY_COMM, MPIerr)
      
      ! initialize input/output files (according to the replicaID)

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
        len=18
      write(filename1,'(a5,i4.4,a1,i4.4,a4)') 
     &   'dyma_',replicaID,'_',my_pe,'.out'
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
        write(6,*)"that specifies the replica ID for this run!"
        write(6,*)"This can be done by placing this number"
        write(6,*)"to a temporary file and redirecting the "
        write(6,*)"standard input from this file: "
        write(6,*)" i.e.: dynapt < replicaNumber.file"
        call alert(name,namel,'Error during read',17,1)
        return
      end
