      subroutine poschi(i1,i2,i3,i4,req12,angleq123,chi1234)
c
c Based on poschi.f V.3 from Ora Schueler (23/12/97)
c
c Placing sidechain atom 1:
c based on chi-angle 1-2-3-4, angle 1-2-3, distance 1-2
c        ichi: number of angle in list of chi-angles

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
c
      integer level
      integer i,j
      integer i1,i2,i3,i4
      double precision req12,angleq123,chi1234
      double precision sin1,cos1
      double precision e34(3),e32(3),e21(3),e3(3)
      double precision pi,phi
c
      character*6 name
      integer namel
c
      pi     = 4.d0*datan(1.d0)
      name   = 'poschi'
      namel  = 6
c
c
c build prtc 4: 3 + e34
c
      do i=1,3
         e34(i) = coor(i,i4)-coor(i,i3)
         e32(i) = coor(i,i2)-coor(i,i3)
      end do
c
c sin1 is used here as a buffer for normalization
      sin1  = e32(1)*e32(1) + e32(2)*e32(2) + e32(3)*e32(3)
      sin1  = 1.d0/dsqrt(sin1)
      do i=1,3
         e32(i) = e32(i)*sin1
      end do
c orthonormalize e34 with respect to e32
      sin1  = e32(1)*e34(1) + e32(2)*e34(2) + e32(3)*e34(3)
      do i=1,3
         e34(i) = e34(i) - sin1*e32(i)
      end do
c normalize the new e34
      sin1  = e34(1)*e34(1) + e34(2)*e34(2) + e34(3)*e34(3)
      sin1  = 1.d0/dsqrt(sin1)
      do i=1,3
         e34(i) = e34(i)*sin1
      end do
c generate a third unit vector e3 by a vector products of e32 and e34
      e3(1) = e32(2)*e34(3) - e32(3)*e34(2)
      e3(2) = e32(3)*e34(1) - e32(1)*e34(3)
      e3(3) = e32(1)*e34(2) - e32(2)*e34(1)
c check that e3 is normalized
      sin1 = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
      if (dabs(sin1-1.d0).gt.1.d-6) then
         level = 1
         call alert(name,namel,'Fishy vector product',20,level)
      end if
c rotate e34 about chi degrees (axis = e32): 
c       e21 = cos(chi) * e34 + sin(chi) * e3
      cos1= dcos(chi1234)
      sin1= dsin(chi1234)
      do i=1,3
         e21(i) = cos1*e34(i) + sin1*e3(i)
      end do
c check that e21 is normalized
      sin1  = e21(1)*e21(1) + e21(2)*e21(2) + e21(3)*e21(3)
      if (dabs(sin1-1.d0).gt.1.d-6) then
         level = 1
         call alert(name,namel,'e21 not normalized 1',20,level)
      end if
      sin1 = e21(1)*e34(1)+e21(2)*e34(2)+e21(3)*e34(3)
c check if angle ok
      if (dabs(sin1 - cos1) .gt. 1.d-6) then
         write(stdo,*)'ichi1234,chi exp,cos obs',
     1        chi1234,cos1,sin1
         call alert(name,namel,'chi not ok',10,1)
      end if
c rotate e21 about (angleq123-pi/2) degrees (axis = e21*e32): 
c     e21 = cos(angleq123-pi/2)*e21 + sin(angleq123-pi/2) * e32
      sin1= dsin(angleq123-pi/2.d0)
      cos1= dcos(angleq123-pi/2.d0)
      do i=1,3
         e21(i) = cos1*e21(i) + sin1*e32(i)
      end do
c check that e21 is normalized
      sin1  = e21(1)*e21(1) + e21(2)*e21(2) + e21(3)*e21(3)
      if (dabs(sin1-1.d0).gt.1.d-6) then
         write(stdo,*)sin1,sin1-1.d0
         level = 1
         call alert(name,namel,'e21 not normalized 2',20,level)
      end if
c build coordinates of prtc 1:
      do i=1,3
         coor(i,i1)=coor(i,i2)+req12*e21(i)
      end do
c
cdeb      call gettors(i1,i2,i3,i4,phi)
c
cdeb      write(*,*)i1,i2,i3,i4,chi1234/pi180
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine gettors(i1,i2,i3,i4,phi)
c
c get the value of the torsional diedral formed by atoms i[1-4]
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
c
      integer i1,i2,i3,i4
      double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
      double precision ux,uy,uz,vx,vy,vz,uv,uu,vv,phi,pi
      parameter(pi=3.14159265358979323844d0)
c
c      
      dx1 = coor(1,i2) - coor(1,i1)
      dy1 = coor(2,i2) - coor(2,i1)
      dz1 = coor(3,i2) - coor(3,i1)
      
      dx2 = coor(1,i3) - coor(1,i2)
      dy2 = coor(2,i3) - coor(2,i2)
      dz2 = coor(3,i3) - coor(3,i2)
      
      dx3 = coor(1,i4) - coor(1,i3)
      dy3 = coor(2,i4) - coor(2,i3)
      dz3 = coor(3,i4) - coor(3,i3)
      
      ux  = dy1*dz2 - dz1*dy2
      uy  = dz1*dx2 - dx1*dz2
      uz  = dx1*dy2 - dy1*dx2
      
      vx  = dy2*dz3 - dz2*dy3
      vy  = dz2*dx3 - dx2*dz3
      vz  = dx2*dy3 - dy2*dx3
      
      uu  = (ux*ux+uy*uy+uz*uz)
      vv  = (vx*vx+vy*vy+vz*vz)
      uv  = (ux*vx+uy*vy+uz*vz)/dsqrt(uu*vv)
c
      if (abs(uv).gt.1.d0)uv=sign(1.d0,uv)
      
      phi = dacos(uv)
      
      dx1 = uy*vz - uz*vy
      dy1 = uz*vx - ux*vz
      dz1 = ux*vy - uy*vx
      
      if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) phi = - phi
c     phi = phi - 180.d0
      if (phi.gt.pi)  phi = phi - 2*pi 
      if (phi.lt.-pi) phi = phi + 2*pi
c
c YP compatible with my convert block
c      write(*,*)i1,i2,i3,i4,phi/pi180
      write(*,*)i1,i2,i3,i4,phi*pi180
c 
      return
      end
c
