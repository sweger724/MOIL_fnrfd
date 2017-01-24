      Subroutine inittime

      implicit none

      return
      end

      Subroutine printcpu(stdo)

      implicit none
      integer stdo
      real*8 dummytime
c
      dummytime=0.0d0    
c
      write(stdo,'(/5x,a,f12.2,a/)')
     $     'Dummy (disabled) CPU time=',dummytime,' (sec)'
      return
      end
