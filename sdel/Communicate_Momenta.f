      subroutine Communicate_Momenta(pp)
        
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/SDEL.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      double precision pp(*)

C      get the momentum of the previous frame
      if (first) pp(1)=p0_initial
      
      if (last)  pp(pseg+2)=p0_final
      
      if (mod(my_pe,2) .eq. 0) then
        call Send_first(1,pp)
        call Receive_first(1,pp)
      else
        call Receive_first(1,pp)
        call Send_first(1,pp)
      endif
 
      return
      end
