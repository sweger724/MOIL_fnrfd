      subroutine vecmin(x, y, n, xmy)

c     IN
c     x, y: equal-length vectors
c     n: length of x & y
c
c     OUT
c     xmy: x - y

      integer i, n
      double precision x(*), y(*), xmy(*)


      do 100 i = 1, n
         xmy(i) = x(i) - y(i)
 100  continue

      end
