      subroutine init_paral()
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
         
         my_pe = 0
         num_pes = 1
         prll_on_off = .false.
      
         return
      end





