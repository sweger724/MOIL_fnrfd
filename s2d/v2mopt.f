      subroutine vec2matav
c
czva: Subroutine vec2matav should compute sum of all forces that act on slow
c     subsystem and sum of all secon derivatives (numerical computation) !!!
c
c  virtual shift of coordinates to compute the second derivatives !!!
c
      double precision dlts2d
      parameter (dlts2d=1.0d-3)
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
      integer i,j,k,l,m,n,npt6
c
      double precision wktmp
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
      npt6=6*npt
c
      do 200 l=1,npt6
            diag(l) = 0.0d0 
 200  continue

      call eforce3()
c
      k=0
c
      do 500 l=1,npts
c
c    average forces
c
            dpot_ave(1,l) = dpot_ave(1,l) + dpot(1,l) 
            dpot_ave(2,l) = dpot_ave(2,l) + dpot(2,l) 
            dpot_ave(3,l) = dpot_ave(3,l) + dpot(3,l) 
c
c   substract the base position part for second derivatives
c 
            k=k+1
c  xx
            d2pt_ave(1,l)=d2pt_ave(1,l) + diag(k)
            k=k+1
c xy & yx
            d2pt_ave(2,l)=d2pt_ave(2,l) + diag(k)
            d2pt_ave(4,l)=d2pt_ave(4,l) + diag(k)
            k=k+1
c xz & zx
            d2pt_ave(3,l)=d2pt_ave(3,l) + diag(k)
            d2pt_ave(7,l)=d2pt_ave(7,l) + diag(k)
            k=k+1
c yy
            d2pt_ave(5,l)=d2pt_ave(5,l) + diag(k)
            k=k+1
c yz & zy
            d2pt_ave(6,l)=d2pt_ave(6,l) + diag(k)
            d2pt_ave(8,l)=d2pt_ave(8,l) + diag(k)
            k=k+1
c zz
            d2pt_ave(9,l)=d2pt_ave(9,l) + diag(k)
 400     continue
 500  continue
c
      return
c
      end







