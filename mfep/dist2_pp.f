      double precision function dist2_pp(x1, x2, p1,nselec,pointr)

c     IN
c     x1 point on the plane 
c     x2 point distant from the plane 
c     p1 plane normal 
c
c     OUT
c     dist2_pp : distance between plane p1 and x2
      implicit none
      integer i, j,k,pointr(*),nselec
      double precision  x1(3,*), x2(3,*),x2mx1(3,nselec),p1(3,*)

      dist2_pp = 0.d0

      do 100 j = 1, nselec
         do 1 k = 1, 3
            x2mx1(k,j) = x2(k,j) - x1(k,j)
 1       continue
 100  continue


      do 200 j = 1, nselec
         do 2 k = 1, 3
            dist2_pp = dist2_pp+(x2mx1(k,j)*p1(k,j))
 2       continue
 200  continue


      return
      end
