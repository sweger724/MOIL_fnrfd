      subroutine vec2matav
c
czva: Subroutine vec2matav should compute sum of all forces that act on slow
c     subsystem and sum of all secon derivatives (numerical computation) !!!
c
c  virtual shift of coordinates to compute the second derivatives !!!
c
      double precision dlts2d
      parameter (dlts2d=1.0d-5)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c    
      character*9 name
      integer namel
      integer i,j,k,l,m,n
c
      double precision d2x(4),d2y(4),d2z(4)
      double precision wktmp
c
c      data d2x/1.0d-5,-1.0d-5,0.0d0,0.0d0/
c      data d2y/0.0d0,1.0d-5,-1.0d-5,0.0d0/
c      data d2z/0.0d0,0.0d0,1.0d-5,-1.0d-5/
c
      data d2x/1.0d-5, 0.0d0, 0.0d0, 0.0d0/
      data d2y/ 0.0d0,1.0d-5, 0.0d0, 0.0d0/
      data d2z/ 0.0d0, 0.0d0,1.0d-5, 0.0d0/
c
      data name/'vec2matav'/
      data namel/9/
c
c     average forces and second derivatives
c     
c the structure of the second order derivatives are
c   
c    dFx/dx     dFx/dy     dFx/dz    
c
c    dFy/dx     dFy/dy     dFy/dz
c
c    dFz/dx     dFz/dy     dFz/dz
c
c  the matrix should be diagonal but because of numerical computations
c it may not obey the symmetry !!!!
c
      do 10 l=1,npts
         coor2(1,l) = coor(1,l)
         coor2(2,l) = coor(2,l)
         coor2(3,l) = coor(3,l)
 10   continue
c
      do 300 m=1,4
c
c translate all slow part for m=1,3. The slow part comes to initial
c position for m=4
c
         do 100 l=1,npts
            coor(1,l) = coor2(1,l)+d2x(m)
            coor(2,l) = coor2(2,l)+d2y(m)
            coor(3,l) = coor2(3,l)+d2z(m)
 100     continue
c
         call eforce()
c
         if(m.eq.4) goto 300
c        
         k=3*(m-1) 
c
c  store the second derivatives
c
         do 200 l=1,npts          
            d2pt_ave(k+1,l)=d2pt_ave(k+1,l)+(dpot(1,l)/dlts2d)        
            d2pt_ave(k+2,l)=d2pt_ave(k+2,l)+(dpot(2,l)/dlts2d)        
            d2pt_ave(k+3,l)=d2pt_ave(k+3,l)+(dpot(3,l)/dlts2d)
 200     continue
c
 300  continue

      do 500 l=1,npts
         do 400 j=1,3
c
c    average forces
c
            dpot_ave(j,l) = dpot_ave(j,l) + dpot(j,l) 
c
c   substract the base position part for second derivatives
c 
            wktmp=dpot(j,l)/dlts2d  
            d2pt_ave(j,l)=d2pt_ave(j,l)-wktmp
            d2pt_ave(j+3,l)=d2pt_ave(j+3,l)-wktmp
            d2pt_ave(j+6,l)=d2pt_ave(j+6,l)-wktmp
 400     continue
 500  continue
c
      return
c
      end







