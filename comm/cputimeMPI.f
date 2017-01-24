      real*8 function elapsedtime()
c     ===========================
c     
c     calculation of the elapsed time in the NT system
c     
c     declaration part
c     ================
c     
      implicit none
      real*8 stwtime
      COMMON/MOILTIME/ stwtime
c
      include 'mpif.h'
c
c     execution part
c     ==============
      elapsedtime=MPI_Wtime()-stwtime
      return 
      end
c
      Subroutine inittime

      implicit none
      real*8 stwtime,elapsedtime
      COMMON/MOILTIME/ stwtime
      stwtime=0.0d0
c
      stwtime=elapsedtime()
      return
      end

      Subroutine printcpu(stdo)

      implicit none
      integer stdo
      real*8 elapsedtime      
c
      write(stdo,'(/5x,a,f12.2,a/)')
     $     'Elapsed CPU time=',elapsedtime(),' (sec)'
      return
      end

      Subroutine printrtc(stdo)

      implicit none 
      integer stdo 
c      real*8 RTC      
c
c      write(stdo,'(/5x,a,f15.2,a/)')
c     $     'Elapsed Greenwich mean time=',RTC(),' (sec)'
      return
      end



