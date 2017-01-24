      subroutine Distance(c1,c2,rms)

      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/SEARCH.BLOCK'

      double precision c1(3,maxpt), c2(3,maxpt),rms

      integer i

        rms = 0.d0

        do i = 1,npt
            rms = rms + mmm(i) * (c1(1,i) - c2(1,i))**2
            rms = rms + mmm(i) * (c1(2,i) - c2(2,i))**2
            rms = rms + mmm(i) * (c1(3,i) - c2(3,i))**2
        end do

        rms = rms / allmass

      return
      end
