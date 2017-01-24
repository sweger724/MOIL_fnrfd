      subroutine radgyr(coor, npt)

c     IN
c     coor: coordinates
c     n: length of x (e.g., 3 * npt)

c     OUT
c     xn: normalized output
c     norm: norm of x

      integer n
      double precision norm, normr, x(*), xn(*)

      norm = 0.d0
      do 100 i = 1, n
         norm = norm + x(i)**2
 100  continue
      norm = sqrt(norm)
      normr = 1.0/norm
      do 200 i = 1, n
         xn(i) = x(i) * normr
 200  continue


      end
