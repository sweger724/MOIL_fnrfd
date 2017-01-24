      subroutine half_harm()

c In  SSBP  simulations with  fixed radius, water oxygens are 
c restrained  with  a half harmonic potential which opearates
c when water molecules cross the fixed boundary of the sphere.

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

c LOCAL
      integer i,j
      double precision delr       
   
      do 10 i=1,nharm
       j = harm(i)
       delr = radist(j)-RMAXM    
       enharm = enharm + 50.d0*delr*delr
       dpot(1,j) = dpot(1,j) + 100.d0*delr*coor(1,j)*rrdist(j)
       dpot(2,j) = dpot(2,j) + 100.d0*delr*coor(2,j)*rrdist(j)
       dpot(3,j) = dpot(3,j) + 100.d0*delr*coor(3,j)*rrdist(j)
 10   continue

      return
      end
