      double precision function dist2_fp(x1, x2, nselec, pointr)

c     IN
c     x1, x2: equal-length vectors (3,*)
c     nselec, pointr: particle selection
c
c     OUT
c     dist2: square of Euclidean distance between x1, x2
      
      integer i, j, pointr(*)
      double precision c1, c2, c3, x1(3,*), x2(3,*)
      
      dist2_fp = 0.d0
      
      do 100 j = 1, nselec
         i = pointr(j)
         do 1 k = 1, 3
            dist2_fp = dist2_fp + (x1(k,i) - x2(k,i))**2
 1       continue
 100  continue

      return
      end
