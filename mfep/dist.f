      double precision function dist2(x1, x2, nselec)

c     IN
c     x1, x2: equal-length vectors (3,*)
c     nselec particle selection
c
c     OUT
c     dist2: square of Euclidean distance between x1, x2
      
      integer i, j, pointr(*)
      double precision c1, c2, c3, x1(3,*), x2(3,*)
      
      dist2 = 0.d0
      
      do 100 j = 1, nselec
    
         do 1 k = 1, 3
            dist2 = dist2 + (x1(k,j) - x2(k,j))**2
 1       continue
 100  continue

      return
      end
