      subroutine cnstrn()
        implicit none
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c
      integer icnst,level
      integer i,j,k,l
      double precision ener,dpoten(3,4)
      character*6 name
      integer namel
c      
c     e_cnst = total constrain energy (ENERGY.BLOCK)
c
c     initialize e_cnst
      e_cnst = 0.d0
      name = 'cnstrn'
      namel = 6
c     initialize loop over constrained torsions in ichunk chunks

      do 100 icnst = 1,ncnst
c
c........call the torsion atoms
         
         i = icnst1(icnst)
         j = icnst2(icnst)
         k = icnst3(icnst)
         l = icnst4(icnst)
c
         if (debug) then
            write(*,*)'       '
            write(*,*)'icnst=',icnst
         end if

c
         if (debug) then
            write(stdo,*)' i j k l',i,j,k,l
            write(stdo,*)' kcns cnseq ',kcns(icnst),cnseq(icnst)
         end if
c
         if (l.gt.0)then
c...........calculating torsional constrain
            call cnsttors(i,j,k,l,kcns(icnst),cnseq(icnst),ener,dpoten)
c         
            e_cnst = e_cnst + ener
c
c...........calculate forces
            dpot(1,i)  = dpot(1,i) + dpoten(1,1)
            dpot(2,i)  = dpot(2,i) + dpoten(2,1)
            dpot(3,i)  = dpot(3,i) + dpoten(3,1)
            
            dpot(1,j)  = dpot(1,j) + dpoten(1,2)
            dpot(2,j)  = dpot(2,j) + dpoten(2,2)
            dpot(3,j)  = dpot(3,j) + dpoten(3,2)
            
            dpot(1,k)  = dpot(1,k) + dpoten(1,3)
            dpot(2,k)  = dpot(2,k) + dpoten(2,3)
            dpot(3,k)  = dpot(3,k) + dpoten(3,3)
            
            dpot(1,l)  = dpot(1,l) + dpoten(1,4)
            dpot(2,l)  = dpot(2,l) + dpoten(2,4)
            dpot(3,l)  = dpot(3,l) + dpoten(3,4)
c
         else if (k.gt.0) then
c
c...........calculating torsional constrain
            call cnstangl(i,j,k,kcns(icnst),cnseq(icnst),ener,dpoten)
c         
            e_cnst = e_cnst + ener
c
c...........calculate forces
            dpot(1,i)  = dpot(1,i) + dpoten(1,1)
            dpot(2,i)  = dpot(2,i) + dpoten(2,1)
            dpot(3,i)  = dpot(3,i) + dpoten(3,1)
            
            dpot(1,j)  = dpot(1,j) + dpoten(1,2)
            dpot(2,j)  = dpot(2,j) + dpoten(2,2)
            dpot(3,j)  = dpot(3,j) + dpoten(3,2)
            
            dpot(1,k)  = dpot(1,k) + dpoten(1,3)
            dpot(2,k)  = dpot(2,k) + dpoten(2,3)
            dpot(3,k)  = dpot(3,k) + dpoten(3,3)
c
         else if (j.gt.0) then
c
c...........calculating distance  constrain
            call cnstdist(i,j,kcns(icnst),cnseq(icnst),ener,dpoten)
c         
            e_cnst = e_cnst + ener
c
c...........calculate forces
            dpot(1,i)  = dpot(1,i) + dpoten(1,1)
            dpot(2,i)  = dpot(2,i) + dpoten(2,1)
            dpot(3,i)  = dpot(3,i) + dpoten(3,1)
            
            dpot(1,j)  = dpot(1,j) + dpoten(1,2)
            dpot(2,j)  = dpot(2,j) + dpoten(2,2)
            dpot(3,j)  = dpot(3,j) + dpoten(3,2)
c 
         else if (i .gt. 0) then 

           call cnstpos(i,kcns(icnst),cnseq(icnst),ener,dpoten)
c
            e_cnst = e_cnst + ener
c
c...........calculate forces
            dpot(MOD(i-1,3)+1,INT((i+2)/3)) = 
     &       dpot(MOD(i-1,3)+1,INT((i+2)/3)) + dpoten(1,1)
        
         else
c
            level = 1
            call alert(name,namel,'something very wrong with constrain',
     &           35,level)
           
         end if
c                 
 100  continue
      
      
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine cnsttors_old(i,j,k,l,kcns,cnseq,ener,dpoten)
c
c     calculate dihedral constrains energies and forces according
c     to T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 
c     with local variant (no delta funcion involved)

c   If phi0=0.
c       V = 2*kimp (1 - cos(phi))
c       dV/dcos(phi) = -2*kimp
c   If abs(phi0)=pi
c       V = 2*kimp (1 + cos(phi))
c       dV/dcos(phi) = 2*kimp
c   If |sin(phi)|.gt.1e-06 use the angle directly
c       V = k ( phi -phi0 ) ^2
c       dV/dcos(phi) = -2*k/sin(phi)
c   if |sin(phi)|.le.1e-06 but far from minimum,
c   set sin(phi) to 1e-06*sign(sin(phi)) arbitrarily
c
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
c
      integer i,j,k,l
      double precision kcns,cnseq,ener,dpoten(3,4)
c
      double precision delta
      double precision co ,phi,den,co1
      double precision uu,vv,uv
      double precision ax,bx,cx, ay,by,cy, az,bz,cz
      double precision ux,uy,uz, vx,vy,vz 
      double precision dx1,dy1,dz1  
      double precision a0x,b0x,c0x, a0y,b0y,c0y, a0z,b0z,c0z
      double precision a1x,b1x, a1y,b1y, a1z,b1z
      double precision a2x,b2x, a2y,b2y, a2z,b2z
      double precision dd1x,dd2x,dd3x,dd4x
      double precision dd1y,dd2y,dd3y,dd4y
      double precision dd1z,dd2z,dd3z,dd4z
      double precision df,yy
      double precision aa,bb,cc,ab,bc,ac
      double precision d1,d2,d3,e1,e2,e3
      double precision pi
      double precision uu2,vv2
c
c     ener = constrain torsion energy (non-acumulated)
c     co = cos(phi)
c
c     define function for dot product calculation
c     dx,ex = dummy variables
c
      double precision dot
      dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3
c      
      pi = 4.d0*datan(1.d0)
c     initialize ener
c
      ener = 0.d0
c     calculate distances
c     
      ax = coor(1,j) - coor(1,i)
      ay = coor(2,j) - coor(2,i)
      az = coor(3,j) - coor(3,i)
c
      bx = coor(1,k) - coor(1,j)
      by = coor(2,k) - coor(2,j)
      bz = coor(3,k) - coor(3,j)
c                  
      cx = coor(1,l) - coor(1,k)
      cy = coor(2,l) - coor(2,k)
      cz = coor(3,l) - coor(3,k)
c                  
c     calculate dot products
      ab = dot(ax,ay,az,bx,by,bz)
      bc = dot(bx,by,bz,cx,cy,cz)
      ac = dot(ax,ay,az,cx,cy,cz)
      aa = dot(ax,ay,az,ax,ay,az)
      bb = dot(bx,by,bz,bx,by,bz)
      cc = dot(cx,cy,cz,cx,cy,cz)
         
                        
c     calculate cos(phi)        
      uu = (aa * bb) - (ab*ab)
      vv = (bb * cc ) - (bc * bc)
      uv = (ab * bc) - (ac * bb)
      den= 1.d0/dsqrt(uu*vv) 
      co = uv * den
         
      if (co.gt.1.d0) co = 1.d0 
      if (co.lt.-1.d0) co = -1.d0 
c
      if (debug) then
         write(stdo,*) 'ab bc ac aa bb cc '
         write(stdo,*) ab,bc,ac,aa,bb,cc
         write(stdo,*) 'uu vv uv den co '
         write(stdo,*) uu,vv,uv,den,co
      end if
c         
c There are 3 options here:
c (a) impeq=0   V = 2*k*(1-cos(phi))
      if (dabs(cnseq).lt.1.d-6) then
         df = -2.d0*kcns
         ener  = -df*(1.d0-co)
c     
c (b) |impeq|=pi  V = 2*k*(1+cos(phi))
      else if (dabs(dabs(cnseq)-pi).lt.1.d-6) then
         df = 2.d0*kcns
         ener  = df*(1.d0+co) 
c     
c (c) impeq=something else V = k(phi-phi0)^2
      else
         
c     now calculate sin(phi) because cos(phi) is symmetric, so
c     we can decide between +-phi.
         
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
         phi = dacos(co)
         
         if(dx1.lt.0) phi = -phi 
         
         
c calculate the deviation from equilibrium
         
         delta = phi-cnseq
         
         if (delta.gt.pi) delta = delta - 2*pi
         if (delta.lt.-pi) delta = delta + 2*pi
         
         df = kcns* delta
         ener = df*delta
         yy = dsin(phi)
         
         if (debug) then
            write(*,*)'phi=',phi/pi180
            write(*,*)'cnseq=',cnseq/pi180
            write(*,*)'delta=',delta/pi180
            write(*,*)'dsin(phi)=',yy
            write(*,*)'energy=',ener
         end if
         
         if ( abs(yy).gt.1.d-06) then
            
            df = -2.d0*df/yy
            
         else   
               
            df = -2.d0*df*1.d6*sign(1.d0,yy)
            
         end if
                 
c
      end if
c
c     calculate derivatives by chain rule
      co1 = 0.5d0 * co * den
         
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
      vv2 = 2.d0*vv
      
      a1x = uu2*(-cc*bx + bc*cx)
      a1y = uu2*(-cc*by + bc*cy)
      a1z = uu2*(-cc*bz + bc*cz)
         
      b1x = uu2*(bb*cx - bc*bx)
      b1y = uu2*(bb*cy - bc*by)
      b1z = uu2*(bb*cz - bc*bz)
         
      a2x = -vv2*(bb*ax - ab*bx)
      a2y = -vv2*(bb*ay - ab*by)
      a2z = -vv2*(bb*az - ab*bz)
         
      b2x = vv2*(aa*bx - ab*ax)
      b2y = vv2*(aa*by - ab*ay)
      b2z = vv2*(aa*bz - ab*az)
         
         
      dd1x = (a0x - a2x*co1)*den
      dd1y = (a0y - a2y*co1)*den
      dd1z = (a0z - a2z*co1)*den
         
      dd2x = (-a0x - b0x - (a1x - a2x - b2x)*co1)*den
      dd2y = (-a0y - b0y - (a1y - a2y - b2y)*co1)*den
      dd2z = (-a0z - b0z - (a1z - a2z - b2z)*co1)*den
         
      dd3x = (b0x - c0x - (-a1x - b1x + b2x)*co1)*den
      dd3y = (b0y - c0y - (-a1y - b1y + b2y)*co1)*den
      dd3z = (b0z - c0z - (-a1z - b1z + b2z)*co1)*den
         
      dd4x = (c0x - b1x*co1)*den 
      dd4y = (c0y - b1y*co1)*den 
      dd4z = (c0z - b1z*co1)*den 
         
c     calculate forces
      dpoten(1,1)  =  df*dd1x
      dpoten(2,1)  =  df*dd1y
      dpoten(3,1)  =  df*dd1z
      
      dpoten(1,2)  =  df*dd2x
      dpoten(2,2)  =  df*dd2y
      dpoten(3,2)  =  df*dd2z
         
      dpoten(1,3)  =  df*dd3x
      dpoten(2,3)  =  df*dd3y
      dpoten(3,3)  =  df*dd3z
         
      dpoten(1,4)  =  df*dd4x
      dpoten(2,4)  =  df*dd4y
      dpoten(3,4)  =  df*dd4z
c         
cdeb
c      if (dabs(dabs(cnseq)-pi).lt.1.d-6) then
c         write(*,*)'phi old, eq',acos(co)/pi*180,cnseq/pi*180, '(pi)'
c      else
c         write(*,*)'phi old, eq',phi/pi*180,cnseq/pi*180
c      end if
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine cnsttors(i,j,k,l,kcns,cnseq,ener,dpoten)
c
c     calculate dihedral constrains energies and forces according
c     to A. Blondel & M. Karplus, J.Comput.Chem, vol 17, 1132 (1996). 
c     It overcomes the singularities of the previous T. Schlick
c     method.
c
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
c
      integer i,j,k,l
      double precision kcns,cnseq,ener,dpoten(3,4)
c
      double precision delta
      double precision co,phi,den
      double precision ax,bx, ay,by, az,bz
      double precision aa,bb,ab
      double precision dx1,dy1,dz1  
      double precision fx,gx,hx, fy,gy,hy, fz,gz,hz
      double precision dd1,dd2,dd3,dd4
      double precision dd1x,dd2x,dd3x,dd4x
      double precision dd1y,dd2y,dd3y,dd4y
      double precision dd1z,dd2z,dd3z,dd4z
      double precision df
      double precision gg,fg,gh,gabs
      double precision pi
      double precision d1,d2,d3,e1,e2,e3
c
c     ener = constrain torsion energy (non-acumulated)
c     co = cos(phi)
c
c     define function for dot product calculation
c     dx,ex = dummy variables
c
      double precision dot
      dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3
c      
      pi = 4.d0*datan(1.d0)
c     initialize ener
c
      ener = 0.d0
c     calculate distances
c     
      fx = coor(1,j) - coor(1,i)
      fy = coor(2,j) - coor(2,i)
      fz = coor(3,j) - coor(3,i)
c
      gx = coor(1,k) - coor(1,j)
      gy = coor(2,k) - coor(2,j)
      gz = coor(3,k) - coor(3,j)
c                  
      hx = coor(1,l) - coor(1,k)
      hy = coor(2,l) - coor(2,k)
      hz = coor(3,l) - coor(3,k)
c                  
c     calculate dot products
      fg = dot(fx,fy,fz,gx,gy,gz)
      gh = dot(gx,gy,gz,hx,hy,hz)
      gg = dot(gx,gy,gz,gx,gy,gz)
      gabs=dsqrt(gg)         
                        
c     calculate cos(phi)
      ax = fy*gz - fz*gy
      ay = fz*gx - fx*gz
      az = fx*gy - fy*gx
         
      bx = gy*hz - gz*hy
      by = gz*hx - gx*hz
      bz = gx*hy - gy*hx
c
      ab = dot(ax,ay,az,bx,by,bz)
      aa = dot(ax,ay,az,ax,ay,az)
      bb = dot(bx,by,bz,bx,by,bz)
c
      den = 1.d0/dsqrt(aa*bb)
      co = ab * den
c           
      if (co.gt.1.d0) co = 1.d0 
      if (co.lt.-1.d0) co = -1.d0 
c
      phi = dacos(co)
c
      dx1 = ay*bz - az*by
      dy1 = az*bx - ax*bz
      dz1 = ax*by - ay*bx
c            
      dx1 = dot(dx1,dy1,dz1,gx,gy,gz)
      if(dx1.lt.0) phi = -phi
c
      if (debug) then
         write(stdo,*) 'fg gh gg '
         write(stdo,*) fg,gh,gg
         write(stdo,*) 'den co '
         write(stdo,*) den,co
      end if
c         
c calculate the deviation from equilibrium
         
      delta = phi-cnseq
         
      if (delta.gt.pi) delta = delta - 2*pi
      if (delta.lt.-pi) delta = delta + 2*pi
         
      df = kcns*delta
      ener = df*delta
      df=2.d0*df
c               
      if (debug) then
         write(*,*)'phi=',phi/pi180
         write(*,*)'cnseq=',cnseq/pi180
         write(*,*)'delta=',delta/pi180
         write(*,*)'energy=',ener
      end if
c         
c
      dd1=df*gabs/aa 
      dd2=df*fg/(aa*gabs)
      dd3=df*gh/(bb*gabs)
      dd4=df*gabs/bb
c
      dd1x=-dd1*ax
      dd1y=-dd1*ay
      dd1z=-dd1*az
c
      dd2x=dd2*ax
      dd2y=dd2*ay
      dd2z=dd2*az
c
      dd3x=-dd3*bx
      dd3y=-dd3*by
      dd3z=-dd3*bz
c
      dd4x=dd4*bx
      dd4y=dd4*by
      dd4z=dd4*bz
c
c calculate forces
      dpoten(1,1)  =  dd1x
      dpoten(2,1)  =  dd1y
      dpoten(3,1)  =  dd1z
      
      dpoten(1,2)  =  dd2x-dd1x-dd3x
      dpoten(2,2)  =  dd2y-dd1y-dd3y
      dpoten(3,2)  =  dd2z-dd1z-dd3z
         
      dpoten(1,3)  =  dd3x-dd2x-dd4x
      dpoten(2,3)  =  dd3y-dd2y-dd4y
      dpoten(3,3)  =  dd3z-dd2z-dd4z
         
      dpoten(1,4)  =  dd4x
      dpoten(2,4)  =  dd4y
      dpoten(3,4)  =  dd4z
c
cdeb
c      write(*,*)'phi new, eq',phi/pi*180,cnseq/pi*180
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine cnstangl(i,j,k,kcns,cnseq,ener,dpoten)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
c
c calculate angle energies and forces
c
c
      integer i,j,k
      double precision kcns,cnseq,ener,dpoten(3,4)
c
      double precision dxi,dyi,dzi,dxj,dyj,dzj,ri2,rj2
      double precision ri,rj,rir,rjr,dxir,dyir,dzir,dxjr,dyjr,dzjr
      double precision cst,at,da,df,tmpx1,tmpy1,tmpz1
      double precision tmpx2,tmpy2,tmpz2

c
c     ener is the constrain angle energy (non-acumulated)
      ener=0.d0
c
c     calculate distances
      dxi = coor(1,i) - coor(1,j)
      dyi = coor(2,i) - coor(2,j)
      dzi = coor(3,i) - coor(3,j)
      dxj = coor(1,k) - coor(1,j)
      dyj = coor(2,k) - coor(2,j)
      dzj = coor(3,k) - coor(3,j)
C
      ri2=dxi*dxi+dyi*dyi+dzi*dzi
      rj2=dxj*dxj+dyj*dyj+dzj*dzj
      ri=dsqrt(ri2)
      rj=dsqrt(rj2)

      rir=1.d0/ri
      rjr=1.d0/rj

      dxir=dxi*rir
      dyir=dyi*rir
      dzir=dzi*rir
      dxjr=dxj*rjr
      dyjr=dyj*rjr
      dzjr=dzj*rjr

      cst=dxir*dxjr+dyir*dyjr+dzir*dzjr
      if (cst.gt.1.d0) cst = 1.d0
      if (cst.lt.-1.d0) cst = -1.d0

      at=dacos(cst)
      da=at-cnseq
      df=kcns*da
      ener=df*da
C@
      df = df+df
      at = dsin(at)
      if (at.lt.1.d-6 .and. at.gt.0.d0) then
         at = 1.d-6
      else if (at.lt.0.d0 .and. at.gt.-1.d-6) then
         at = -1.d-6
      end if
      df=-df/at
      
      tmpx1 = rir*(dxjr-cst*dxir)
      tmpy1 = rir*(dyjr-cst*dyir)
      tmpz1 = rir*(dzjr-cst*dzir)
      
      tmpx2 = rjr*(dxir-cst*dxjr)
      tmpy2 = rjr*(dyir-cst*dyjr)
      tmpz2 = rjr*(dzir-cst*dzjr)

      dpoten(1,1) = df*tmpx1
      dpoten(2,1) = df*tmpy1
      dpoten(3,1) = df*tmpz1

      dpoten(1,2) = -df*(tmpx1+tmpx2)
      dpoten(2,2) = -df*(tmpy1+tmpy2)
      dpoten(3,2) = -df*(tmpz1+tmpz2)

      dpoten(1,3) = df*tmpx2
      dpoten(2,3) = df*tmpy2
      dpoten(3,3) = df*tmpz2

      return

      end
c
c----------------------------------------------------------------------
c
      subroutine cnstdist(i,j,kcns,cnseq,ener,dpoten)
c
c     calculate distance constrains energies and forces.
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c
      integer i,j
      double precision kcns,cnseq,ener,dpoten(3,4)
c
      double precision ax,ay,az,aa,dij
      double precision delta,df
c
c     ener = constrain distance energy (non-acumulated)
c
c     calculate distances
c     
      ax = coor(1,j) - coor(1,i)
      ay = coor(2,j) - coor(2,i)
      az = coor(3,j) - coor(3,i)
c
c     calculate dot products
      aa = ax*ax+ay*ay+az*az
      dij= dsqrt(aa)
c                        
      if (debug) then
         write(stdo,*) 'distance', i,j
         write(stdo,*) dij
         write(stdo,*) 'constrain', '  force const.'
         write(stdo,*) cnseq,kcns
      end if
c
c calculate the deviation from equilibrium   
      delta = dij-cnseq 
      df = kcns*delta
      ener = df*delta
      df= -2.d0*df/dij
c         
c     calculate forces
      dpoten(1,1)  =  df*ax
      dpoten(2,1)  =  df*ay
      dpoten(3,1)  =  df*az
      
      dpoten(1,2)  =  -dpoten(1,1)
      dpoten(2,2)  =  -dpoten(2,1)
      dpoten(3,2)  =  -dpoten(3,1)
      if (debug) then
         write(*,*)'energy=',ener
         write(*,*)'forces=',(dpoten(i,1),i=1,3)
      end if
c         
      return
      end


c
c----------------------------------------------------------------------
c
      subroutine cnstpos(i,kcns,cnseq,ener,dpoten)
c
c     calculate posiiton constrains energies and forces.
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c
      integer i,iii,lll
      double precision kcns,cnseq,ener,dpoten(3,4)
c
      double precision ax,ay,az,aa,dij
      double precision delta,df
c
c     ener = constrain distance energy (non-acumulated)
c
      iii = INT((i+2)/3)
      lll = MOD(i-1,3)+1
C      write(6,*),i,iii,lll

c calculate the deviation from equilibrium
      delta = coor(lll,iii) - cnseq
      df = kcns*delta
      ener = df*delta
      df= 2.d0*df
c
c     calculate forces
      dpoten(1,1)  =  df
      
      if (debug) then
         write(*,*)'energy=',ener,coor(lll,iii),cnseq
         write(*,*)'forces=',(dpoten(i,1),i=1,3)
      end if
c
      return
      end
