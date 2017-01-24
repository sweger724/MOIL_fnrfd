      subroutine Communicate_dStmp(ss)
        
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/SDP.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      double precision ss(*)

C      get the dSmtp of the previous frame
      if (first) ss(1)=dStmp_initial
      
      if (last)  ss(pseg+2)=dStmp_final
      
      if (mod(my_pe,2) .eq. 0) then
        call Send_first(1,ss)
        call Receive_first(1,ss)
      else
        call Receive_first(1,ss)
        call Send_first(1,ss)
      endif
 
      return
      end
