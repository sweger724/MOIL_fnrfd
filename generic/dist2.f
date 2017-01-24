      subroutine dist2(x1, x2, n, d)
      implicit none
c     IN
c     x1, x2: equal-length vectors
c     n: length of x1, x2
c
c     OUT
c     d: square of Euclidean distance between x1, x2

      integer i, n
      double precision d, x1(*), x2(*)

      d = 0.d0
      do i = 1, n
         d = d + (x1(i) - x2(i))**2
      end do 

      end
