      subroutine vecmin_fp(x, y, xmy, nselec, pointr)

      implicit none

      integer i, j, k, nselec, pointr(*)
      double precision x(3,*), y(3,*), xmy(3,*)

      do 100 j = 1, nselec
         i = pointr(j)
         do 1 k = 1,3
            xmy(k,j) = x(k,i) - y(k,j)
 1       continue
 100  continue

      end
