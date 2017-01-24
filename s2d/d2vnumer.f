      subroutine d2vnumer(d2v)
c
c Calculation of the second derivatives by taking the  
c numerical derivatives of the force.
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/GETVEC.BLOCK'


c      character*5 name
c      integer namel,level
      integer i,j,kx,ky,kz,jj
      
      double precision delta,factor
      double precision d2v(3*maxpt,3*maxpt)
c
c.....initializations
      call nbondm()
      if (esymyes) call syminit()
c         
      delta = 1.d-5
      factor = 1.d0/(2.d0*delta)
c
c
      do 22 i=1,npt
c
c........F(i_x-dx)
         coor(1,i) = coor(1,i) - delta
         call eforce()
         do 23 j=1,npt
            coor2(1,j) = dpot(1,j)
            coor2(2,j) = dpot(2,j)
            coor2(3,j) = dpot(3,j)
 23      continue
c
c........F(i_x+dx)
         coor(1,i) = coor(1,i) + 2.d0*delta
         call eforce()
c
c........back to i_x
         coor(1,i) = coor(1,i) - delta
c
c........pass to the 3*n representation
         kx = 3*(i-1)+1
c
c........dF/di_x
         do 24 j=1,npt
            jj = 3*(j-1)+1
            d2v(jj,kx)   = factor*(dpot(1,j)-coor2(1,j))
            d2v(jj+1,kx) = factor*(dpot(2,j)-coor2(2,j))
            d2v(jj+2,kx) = factor*(dpot(3,j)-coor2(3,j))
 24      continue
c
c........F(i_y-dy)
         coor(2,i) = coor(2,i) - delta
         call eforce()
         do 25 j=1,npt
            coor2(1,j) = dpot(1,j)
            coor2(2,j) = dpot(2,j)
            coor2(3,j) = dpot(3,j)
 25      continue
c
c........F(i_y+dy)
         coor(2,i) = coor(2,i) + 2.d0*delta
         call eforce()
c
c........back to i_y
         coor(2,i) = coor(2,i) - delta
c
c........pass to the 3*n representation
         ky = 3*(i-1)+2
c
c........dF/di_y
         do 26 j=1,npt
            jj = 3*(j-1)+1
            d2v(jj,ky)   = factor*(dpot(1,j)-coor2(1,j))
            d2v(jj+1,ky) = factor*(dpot(2,j)-coor2(2,j))
            d2v(jj+2,ky) = factor*(dpot(3,j)-coor2(3,j))
 26      continue
c
c........F(i_z-dz)
         coor(3,i) = coor(3,i) - delta
         call eforce()
         do 27 j=1,npt
            coor2(1,j) = dpot(1,j)
            coor2(2,j) = dpot(2,j)
            coor2(3,j) = dpot(3,j)
 27      continue
c
c........F(i_z+dz)
         coor(3,i) = coor(3,i) + 2.d0*delta
         call eforce()
c
c........back to i_z
         coor(3,i) = coor(3,i) - delta
c
c........pass to the 3*n representation
         kz = i*3
c
c........dF/di_z
         do 28 j=1,npt
            jj = 3*(j-1)+1
            d2v(jj,kz)   = factor*(dpot(1,j)-coor2(1,j))
            d2v(jj+1,kz) = factor*(dpot(2,j)-coor2(2,j))
            d2v(jj+2,kz) = factor*(dpot(3,j)-coor2(3,j))
 28      continue
 22   continue

c      call house(d2v,3*npt,3*npt,eigenv,work,i)
c      write(stdo,*)' Eigenvalues '
c      write(stdo,100)(eigenv(j),j=1,3*npt)
c      write(stdo,*)' Error = ',i

      return
      end



