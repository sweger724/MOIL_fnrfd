      subroutine torsion ( i, j, k, l, ang, dang, doderiv )

c     IMPORTANT: Assumes arith is off (false).

      implicit none


      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'

      integer i, j, k, l
      double precision ax,ay,az,bx,by,bz,cx,cy,cz,ab,bc,ac,aa,bb,cc
      double precision uu,vv,uv,co,ux,uy,uz,vx,vy,vz,dx1,dy1,dz1
      double precision ang,den,df,yy,pi
      double precision dot,d1,d2,d3,e1,e2,e3
      double precision co1,a0x,a0y,a0z,b0x,b0y,b0z,c0x,c0y,c0z
      double precision uu2,a1x,a1y,a1z,b1x,b1y,b1z
      double precision vv2,a2x,a2y,a2z,b2x,b2y,b2z
      logical doderiv

c     gradient
      double precision dang(3,4)

c     MACRO
      dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3

      pi=4.d0*datan(1.d0)

c     calculate distances
      ax = coor(1,j) - coor(1,i)
      ay = coor(2,j) - coor(2,i)
      az = coor(3,j) - coor(3,i)

      bx = coor(1,k) - coor(1,j)
      by = coor(2,k) - coor(2,j)
      bz = coor(3,k) - coor(3,j)

      cx = coor(1,l) - coor(1,k)
      cy = coor(2,l) - coor(2,k)
      cz = coor(3,l) - coor(3,k)

c     calculate dot products
      ab = dot(ax,ay,az,bx,by,bz)
      bc = dot(bx,by,bz,cx,cy,cz)
      ac = dot(ax,ay,az,cx,cy,cz)
      aa = dot(ax,ay,az,ax,ay,az)
      bb = dot(bx,by,bz,bx,by,bz)
      cc = dot(cx,cy,cz,cx,cy,cz)


c     calculate cos(ang)	
      uu = (aa * bb) - (ab * ab)
      vv = (bb * cc) - (bc * bc)
      uv = (ab * bc) - (ac * bb)
      den= 1.d0 / dsqrt(uu * vv)
      co = uv * den

      if (co.gt.1.d0) co = 1.d0 
      if (co.lt.-1.d0) co = -1.d0 

      
c     now calculate sin(ang) because cos(ang) is symmetric, so
c     we can decide between +-ang.
      ux = ay*bz - az*by
      uy = az*bx - ax*bz
      uz = ax*by - ay*bx
      
      vx = by*cz - bz*cy
      vy = bz*cx - bx*cz
      vz = bx*cy - by*cx
      
      dx1 = uy*vz - uz*vy
      dy1 = uz*vx - ux*vz
      dz1 = ux*vy - uy*vx
      
      dx1 = dot(dx1,dy1,dz1,bx,by,bz)

      ang = dacos(co)
      if(dx1.lt.0) ang = -ang

      if ( .not. doderiv ) return


c     ------ derivative calculation

      df = 1.d0
      yy = dsin(ang)
      if ( abs(yy).gt.1.d-06 ) then
         df = -df/yy
      else
         if((ang.gt.0.d0 .and. ang.lt.(pi/2.d0)) .or. (ang.lt.0.d0 .and.
     $        ang.gt.(-pi/2.d0)) ) then 
            df = df*1d06
         else
            df = -df*1d06
         end if	
      end if
		 


c     calculate derivatives by chain rule
      co1 = 0.5d0*co*den

      a0x = -bc*bx + bb*cx
      a0y = -bc*by + bb*cy
      a0z = -bc*bz + bb*cz
      
      b0x = ab*cx + bc*ax -2.d0*ac*bx
      b0y = ab*cy + bc*ay -2.d0*ac*by
      b0z = ab*cz + bc*az -2.d0*ac*bz
      
      c0x = ab*bx - bb*ax
      c0y = ab*by - bb*ay
      c0z = ab*bz - bb*az
      
      uu2 = 2.d0*uu

      a1x = uu2*(-cc*bx + bc*cx)
      a1y = uu2*(-cc*by + bc*cy)
      a1z = uu2*(-cc*bz + bc*cz)
      
      b1x = uu2*(bb*cx - bc*bx)
      b1y = uu2*(bb*cy - bc*by)
      b1z = uu2*(bb*cz - bc*bz)

      vv2 = 2.d0*vv
      
      a2x = -vv2*(bb*ax - ab*bx)
      a2y = -vv2*(bb*ay - ab*by)
      a2z = -vv2*(bb*az - ab*bz)
      
      b2x = vv2*(aa*bx - ab*ax)
      b2y = vv2*(aa*by - ab*ay)
      b2z = vv2*(aa*bz - ab*az)

      den = den * df

      dang(1,1) = (a0x - a2x*co1)*den
      dang(2,1) = (a0y - a2y*co1)*den
      dang(3,1) = (a0z - a2z*co1)*den

      dang(1,2) = (-a0x - b0x - (a1x - a2x - b2x)*co1)*den
      dang(2,2) = (-a0y - b0y - (a1y - a2y - b2y)*co1)*den
      dang(3,2) = (-a0z - b0z - (a1z - a2z - b2z)*co1)*den

      dang(1,3) = (b0x - c0x - (-a1x - b1x + b2x)*co1)*den
      dang(2,3) = (b0y - c0y - (-a1y - b1y + b2y)*co1)*den
      dang(3,3) = (b0z - c0z - (-a1z - b1z + b2z)*co1)*den

      dang(1,4) = (c0x - b1x*co1)*den 
      dang(2,4) = (c0y - b1y*co1)*den 
      dang(3,4) = (c0z - b1z*co1)*den


      return
      end
