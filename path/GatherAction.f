      subroutine GatherAction(S)

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      double precision S
      integer logN,log2
        
      if(.not. paral) return
      
C Parallel algorithm here:

C    compute log2(proc)
      logN=log2(numprocs)

      call Gather_Data(logN,S,numprocs)
C   redistribute computed average
C   call Distribute_Data(logN,S,numprocs)

      return
      end
