        subroutine prepare_random_numbers()

          implicit none

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/PARALLEL.BLOCK'
          include 'COMMON/PT.BLOCK'
          include 'mpif.h'

          integer i

          ! draw sequence of random numbers
          ! it is important that all processors in replicaID
          ! draw randomSize random numbers, because random number 
          ! generators must be synchronized within a replica
          if (replicaID .eq. 0) call RANLUX(RANDOM_ARRAY,randomSize)
          
          ! broadcast this sequence to all processors/replicas
          call mpi_bcast(RANDOM_ARRAY,randomSize,MPI_REAL,
     &                   0,MPI_COMM_WORLD, MPIerr)
          
          ! initialize position in the sequence
          random_index = 0
          
          return
        end
