      subroutine eimphi2_pcs()
c
c Calculate second derivative of improper torsion energy.
c
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c
c
c calculate second derivaives of the 
c improper dihedral energies and forces according
c to T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 
c with local variant (no delta funcion involved)

c If |sin(phi)|.gt.1e-03 use the angle directly
c    V = kimp ( phi -impeq ) ^2
c    dV/dcos(phi) = -2*k/sin(phi)

c If |sin(phi)|.lt.1e-03   and sin(impeq)=0.
c    V = 2*kimp (1 - cos(phi))
c    dV/dcos(phi) = -2*kimp

c If |sin(phi)|.lt.1e-03    and sin(impeq)=pi
c    V = 2*kimp (1 + cos(phi))
c    dV/dcos(phi) = 2*kimp

c We dont take care of very large displacements in which impeq !=0
c but phi is. In this case a warning is issued and no calculation is
c done
c
      integer namel,level
      character*7 name
      double precision phi,delta,pi,scalar1,scalar2
      double precision ux,uy,uz,vx,vy,vz,dx1,dy1,dz1

      integer iphi,iiphi,k,l,m,n,kk,ll,mm,nn
      double precision ax,ay,az,bx,by,bz,cx,cy,cz
      double precision ab,bc,ac,aa,bb,cc,uu,vv,uv
      double precision den,den2,den3,co

      double precision dxkuv,dykuv,dzkuv
      double precision dxluv,dyluv,dzluv
      double precision dxmuv,dymuv,dzmuv
      double precision dxnuv,dynuv,dznuv

      double precision dxkuu,dykuu,dzkuu
      double precision dxluu,dyluu,dzluu
      double precision dxmuu,dymuu,dzmuu

      double precision dxlvv,dylvv,dzlvv
      double precision dxmvv,dymvv,dzmvv
      double precision dxnvv,dynvv,dznvv

      double precision uvlxlx,uvlyly,uvlzlz,uvlxly,uvlxlz,uvlylz
      double precision uvmxmx,uvmymy,uvmzmz,uvmxmy,uvmxmz,uvmymz
      double precision uvkxlx,uvkyly,uvkzlz,uvkxly,uvkxlz,uvkylx,
     1     uvkylz,uvkzlx,uvkzly
      double precision uvkxmx,uvkymy,uvkzmz,uvkxmy,uvkxmz,uvkymx,
     1     uvkymz,uvkzmx,uvkzmy
      double precision uvkxnx,uvkyny,uvkznz,uvkxny,uvkxnz,uvkynx,
     1     uvkynz,uvkznx,uvkzny
      double precision uvlxmx,uvlymy,uvlzmz,uvlxmy,uvlxmz,uvlymx,
     1     uvlymz,uvlzmx,uvlzmy
      double precision uvlxnx,uvlyny,uvlznz,uvlxny,uvlxnz,uvlynx,
     1     uvlynz,uvlznx,uvlzny
      double precision uvmxnx,uvmyny,uvmznz,uvmxny,uvmxnz,uvmynx,
     1     uvmynz,uvmznx,uvmzny


      double precision uukxkx,uukyky,uukzkz,uukxky,uukxkz,uukykz
      double precision uukxlx,uukyly,uukzlz,uukxly,uukxlz,uukylx,
     1     uukylz,uukzlx,uukzly
      double precision uukxmx,uukymy,uukzmz,uukxmy,uukxmz,uukymx,
     1     uukymz,uukzmx,uukzmy
      double precision uulxlx,uulyly,uulzlz,uulxly,uulxlz,uulylz
      double precision uulxmx,uulymy,uulzmz,uulxmy,uulxmz,uulymx,
     1     uulymz,uulzmx,uulzmy
      double precision uumxmx,uumymy,uumzmz,uumxmy,uumxmz,uumymz

      double precision vvlxlx,vvlyly,vvlzlz,vvlxly,vvlxlz,vvlylz
      double precision vvlxmx,vvlymy,vvlzmz,vvlxmy,vvlxmz,vvlymx,
     1     vvlymz,vvlzmx,vvlzmy
      double precision vvlxnx,vvlyny,vvlznz,vvlxny,vvlxnz,vvlynx,
     1     vvlynz,vvlznx,vvlzny
      double precision vvmxmx,vvmymy,vvmzmz,vvmxmy,vvmxmz,vvmymz
      double precision vvmxnx,vvmyny,vvmznz,vvmxny,vvmxnz,vvmynx,
     1     vvmynz,vvmznx,vvmzny
      double precision vvnxnx,vvnyny,vvnznz,vvnxny,vvnxnz,vvnynz


      double precision cokx,coky,cokz,colx,coly,colz,comx,comy,comz
      double precision conx,cony,conz
      double precision cokxkx,cokyky,cokzkz,cokxky,cokxkz,cokykz
      double precision colxlx,colyly,colzlz,colxly,colxlz,colylz
      double precision comxmx,comymy,comzmz,comxmy,comxmz,comymz
      double precision conxnx,conyny,conznz,conxny,conxnz,conynz

      double precision cokxlx,cokxly,cokxlz,cokylx,cokyly,cokylz
      double precision cokzlx,cokzly,cokzlz
      double precision cokxmx,cokxmy,cokxmz,cokymx,cokymy,cokymz
      double precision cokzmx,cokzmy,cokzmz
      double precision cokxnx,cokxny,cokxnz,cokynx,cokyny,cokynz
      double precision cokznx,cokzny,cokznz
      double precision colxmx,colxmy,colxmz,colymx,colymy,colymz
      double precision colzmx,colzmy,colzmz
      double precision colxnx,colxny,colxnz,colynx,colyny,colynz
      double precision colznx,colzny,colznz
      double precision comxnx,comxny,comxnz,comynx,comyny,comynz
      double precision comznx,comzny,comznz

      pi = 4.d0*datan(1.d0)
      namel = 7
      name  = 'eimphi2'

      do 100 iphi = ifirst,ilast
c
c........call the torsion atoms
         k = iimp1(iphi)
         l = iimp2(iphi)
         m = iimp3(iphi)
         n = iimp4(iphi)


c........calculate distances
         ax = coor(1,l) - coor(1,k)
         ay = coor(2,l) - coor(2,k)
         az = coor(3,l) - coor(3,k)

         bx = coor(1,m) - coor(1,l)
         by = coor(2,m) - coor(2,l)
         bz = coor(3,m) - coor(3,l)

         cx = coor(1,n) - coor(1,m)
         cy = coor(2,n) - coor(2,m)
         cz = coor(3,n) - coor(3,m)

c........calculate dot products
         ab = ax*bx + ay*by + az*bz
         bc = bx*cx + by*cy + bz*cz
         ac = ax*cx + ay*cy + az*cz
         aa = ax*ax + ay*ay + az*az
         bb = bx*bx + by*by + bz*bz
         cc = cx*cx + cy*cy + cz*cz


c........calculate cos(phi)      
         uu  = (aa * bb) - (ab*ab)
         vv  = (bb * cc ) - (bc * bc)
         uv  = (ab * bc) - (ac * bb)
         den  = 1.d0/dsqrt(uu*vv) 
         den2 = den*den
         den3 = den2*den
         co   = uv * den
         if (co.gt.1.d0) co = 1.d0
         if (co.lt.-1.d0) co = -1.d0
c
C
C Alfredo:
c       ileana-- one single option, e = k(1+cos(n*phi -gamma))
c       ileana-- since always n=2 and gamma = pi in the new Amber ALL.PROP, e =
C       k - kcos(2phi)
C
                        if(.not.arith)then
c                          write(*,*)'use the old eimphi if not arith'
c there are three possibilities here:
c (a) The equilibrium position for the improper torsion is 0 (zero)
c     In this case the energy function is  V=2*kimp*(1-cos(phi))
c (b) The equilibrium position for the improper torsion is 3.14.. (pi)
c     In this case the energy function is  V=2*kimp*(1+cos(phi))
c (c) The equilibrium position is far from zero.
c      Energy: V = kimp*(phi-impeq)^2
c
c In cases (a) and (b) there is no need to calculate phi and/or sin(phi)
c which make the cvalculations more straightforward. Unfortunately no
c such simple solution is available for angles far from zero.
c

c Case (a)
         if (dabs(impeq(iphi)).lt.1.d-6) then
            scalar1 = 0.d0
            scalar2 = -2.d0*kimp(iphi)
c Case (b)
         else if (dabs(dabs(impeq(iphi))-pi).lt.1.d-6) then
            scalar1 = 0.d0
            scalar2 = 2.d0*kimp(iphi)
c Case (c)
         else
            phi = dacos(co)

c      now calculate sin(phi) because cos(phi) is symmetric, so
c      we can decide between +-phi.

            ux = ay*bz - az*by
            uy = az*bx - ax*bz
            uz = ax*by - ay*bx

            vx = by*cz - bz*cy
            vy = bz*cx - bx*cz
            vz = bx*cy - by*cx

            dx1 = uy*vz - uz*vy
            dy1 = uz*vx - ux*vz
            dz1 = ux*vy - uy*vx

            dx1 = dx1*bx+dy1*by+dz1*bz
            if(dx1.lt.0) phi = -phi

            sinphi  = sin(phi)
            if (sinphi*sinphi.lt.1.d-12) then
               level = 1
               call alert(name,namel,' sinphi=0 ',10,level)
            end if
            sinphi = 1.d0/sinphi
            delta = phi-impeq(iphi)
            scalar1 = 2.d0*kimp(iphi)*sinphi*sinphi*(
     1                1.d0-delta*co*sinphi)
            scalar2 = -2.d0*kimp(iphi)*sinphi*delta
         end if
C arith:
                           else

            phi = dacos(co)

c      now calculate sin(phi) because cos(phi) is symmetric, so
c      we can decide between +-phi.

            ux = ay*bz - az*by
            uy = az*bx - ax*bz
            uz = ax*by - ay*bx
 
            vx = by*cz - bz*cy
            vy = bz*cx - bx*cz
            vz = bx*cy - by*cx
 
            dx1 = uy*vz - uz*vy
            dy1 = uz*vx - ux*vz
            dz1 = ux*vy - uy*vx
 
            dx1 = dx1*bx+dy1*by+dz1*bz
            if(dx1.lt.0) phi = -phi 

            scalar1 = -4.0d0*kimp(iphi)
            scalar2 = -4.0d0*kimp(iphi)*co
         end if
C End of Alfredo modifications

c........calculate derivatives of scalar products 
c uv
         dxkuv = -bx*bc + cx*bb
         dykuv = -by*bc + cy*bb
         dzkuv = -bz*bc + cz*bb

         dxluv = (bx-ax)*bc-cx*(ab+bb)+2.d0*bx*ac
         dyluv = (by-ay)*bc-cy*(ab+bb)+2.d0*by*ac
         dzluv = (bz-az)*bc-cz*(ab+bb)+2.d0*bz*ac

c@
         dxmuv = ax*(bc+bb)+(cx-bx)*ab-2.d0*bx*ac
         dymuv = ay*(bc+bb)+(cy-by)*ab-2.d0*by*ac
         dzmuv = az*(bc+bb)+(cz-bz)*ab-2.d0*bz*ac

         dxnuv = bx*ab-ax*bb
         dynuv = by*ab-ay*bb
         dznuv = bz*ab-az*bb

c uu (note derivative with respect to n is zero and not written)
         dxkuu = 2.d0*(bx*ab-bb*ax)
         dykuu = 2.d0*(by*ab-bb*ay)
         dzkuu = 2.d0*(bz*ab-bb*az)

         dxluu = 2.d0*(ax*bb-bx*aa-(bx-ax)*ab)
         dyluu = 2.d0*(ay*bb-by*aa-(by-ay)*ab)
         dzluu = 2.d0*(az*bb-bz*aa-(bz-az)*ab)

         dxmuu = 2.d0*(bx*aa-ax*ab)
         dymuu = 2.d0*(by*aa-ay*ab)
         dzmuu = 2.d0*(bz*aa-az*ab)

c vv (note derivative with respect to k is zero and not written)
         dxlvv = 2.d0*(cx*bc-bx*cc)
         dylvv = 2.d0*(cy*bc-by*cc)
         dzlvv = 2.d0*(cz*bc-bz*cc)

         dxmvv = 2.d0*(bx*cc-cx*bb-(cx-bx)*bc)
         dymvv = 2.d0*(by*cc-cy*bb-(cy-by)*bc)
         dzmvv = 2.d0*(bz*cc-cz*bb-(cz-bz)*bc)
         
         dxnvv = 2.d0*(cx*bb-bx*bc)
         dynvv = 2.d0*(cy*bb-by*bc)
         dznvv = 2.d0*(cz*bb-bz*bc)

c Second derivatives of scalar products
c uv (kk,nn=0)
         uvlxlx = 2.d0*(-bc+ax*cx-ac+bx*cx)
         uvlyly = 2.d0*(-bc+ay*cy-ac+by*cy)
         uvlzlz = 2.d0*(-bc+az*cz-ac+bz*cz)
         uvlxly = ax*cy+ay*cx+cx*by+cy*bx
         uvlxlz = ax*cz+az*cx+cx*bz+cz*bx
         uvlylz = ay*cz+az*cy+cy*bz+cz*by

         uvmxmx = 2.d0*(-ab+ax*cx-ac+ax*bx)
         uvmymy = 2.d0*(-ab+ay*cy-ac+ay*by)
         uvmzmz = 2.d0*(-ab+az*cz-ac+az*bz)
         uvmxmy = ax*cy+ay*cx+ax*by+ay*bx
         uvmxmz = ax*cz+az*cx+ax*bz+az*bx
         uvmymz = ay*cz+az*cy+ay*bz+az*by

         uvkxlx = bc-bx*cx
         uvkyly = bc-by*cy
         uvkzlz = bc-bz*cz
         uvkxly = bx*cy-2.d0*cx*by
         uvkxlz = bx*cz-2.d0*cx*bz
         uvkylx = by*cx-2.d0*cy*bx
         uvkylz = by*cz-2.d0*cy*bz
         uvkzlx = bz*cx-2.d0*cz*bx
         uvkzly = bz*cy-2.d0*cz*by

         uvkxmx = -bc+bx*cx+bx*bx-bb
         uvkymy = -bc+by*cy+by*by-bb
         uvkzmz = -bc+bz*cz+bz*bz-bb
         uvkxmy = -bx*cy+bx*by+2.d0*cx*by
         uvkxmz = -bx*cz+bx*bz+2.d0*cx*bz
         uvkymx = -by*cx+by*bx+2.d0*cy*bx
         uvkymz = -by*cz+by*bz+2.d0*cy*bz
         uvkzmx = -bz*cx+bz*bx+2.d0*cz*bx
         uvkzmy = -bz*cy+bz*by+2.d0*cz*by

         uvkxnx = -bx*bx + bb
         uvkyny = -by*by + bb
         uvkznz = -bz*bz + bb
         uvkxny = -bx*by
         uvkxnz = -bx*bz
         uvkynx = uvkxny
         uvkynz = -by*bz
         uvkznx = uvkxnz
         uvkzny = uvkynz

         uvlxmx = bc+(bx-ax)*(cx-bx)+ab-
     1        ax*cx+bb+2.d0*(ac-cx*bx-ax*bx)
         uvlymy = bc+(by-ay)*(cy-by)+ab-
     1        ay*cy+bb+2.d0*(ac-cy*by-ay*by)
         uvlzmz = bc+(bz-az)*(cz-bz)+ab-
     1        az*cz+bb+2.d0*(ac-cz*bz-az*bz)
         uvlxmy =(bx-ax)*(cy-by)-ay*cx-2.d0*(cx*by+ay*bx)
         uvlxmz =(bx-ax)*(cz-bz)-az*cx-2.d0*(cx*bz+az*bx)
         uvlymx =(by-ay)*(cx-bx)-ax*cy-2.d0*(cy*bx+ax*by)
         uvlymz =(by-ay)*(cz-bz)-az*cy-2.d0*(cy*bz+az*by)
         uvlzmx =(bz-az)*(cx-bx)-ax*cz-2.d0*(cz*bx+ax*bz)
         uvlzmy =(bz-az)*(cy-by)-ay*cz-2.d0*(cz*by+ay*bz)

         uvlxnx = (bx+ax)*bx-ab-bb
         uvlyny = (by+ay)*by-ab-bb
         uvlznz = (bz+az)*bz-ab-bb
         uvlxny = (bx-ax)*by+2.d0*ay*bx
         uvlxnz = (bx-ax)*bz+2.d0*az*bx
         uvlynx = (by-ay)*bx+2.d0*ax*by
         uvlynz = (by-ay)*bz+2.d0*az*by
         uvlznx = (bz-az)*bx+2.d0*ax*bz
         uvlzny = (bz-az)*by+2.d0*ay*bz

         uvmxnx = -ax*bx+ab
         uvmyny = -ay*by+ab
         uvmznz = -az*bz+ab
         uvmxny = ax*by - 2.d0*bx*ay
         uvmxnz = ax*bz - 2.d0*bx*az
         uvmynx = ay*bx - 2.d0*by*ax
         uvmynz = ay*bz - 2.d0*by*az
         uvmznx = az*bx - 2.d0*bz*ax
         uvmzny = az*by - 2.d0*bz*ay

c uu (all n components are zero)
         uukxkx = 2.d0*(bb - bx*bx)
         uukyky = 2.d0*(bb - by*by)
         uukzkz = 2.d0*(bb - bz*bz)
         uukxky = -2.d0*bx*by
         uukxkz = -2.d0*bx*bz
         uukykz = -2.d0*by*bz

         uukxlx = 2.d0*(-bb-ab+bx*(bx+ax))
         uukyly = 2.d0*(-bb-ab+by*(by+ay))
         uukzlz = 2.d0*(-bb-ab+bz*(bz+az))
         uukxly = 2.d0*(2.d0*by*ax+bx*(by-ay))
         uukxlz = 2.d0*(2.d0*bz*ax+bx*(bz-az))
         uukylx = 2.d0*(2.d0*bx*ay+by*(bx-ax))
         uukylz = 2.d0*(2.d0*bz*ay+by*(bz-az))
         uukzlx = 2.d0*(2.d0*bx*az+bz*(bx-ax))
         uukzly = 2.d0*(2.d0*by*az+bz*(by-ay))

         uukxmx = 2.d0*(ab-ax*bx)
         uukymy = 2.d0*(ab-ay*by)
         uukzmz = 2.d0*(ab-az*bz)
         uukxmy = 2.d0*(bx*ay-2.d0*ax*by)
         uukxmz = 2.d0*(bx*az-2.d0*ax*bz)
         uukymx = 2.d0*(by*ax-2.d0*ay*bx)
         uukymz = 2.d0*(by*az-2.d0*ay*bz)
         uukzmx = 2.d0*(bz*ax-2.d0*az*bx)
         uukzmy = 2.d0*(bz*ay-2.d0*az*by)

         uulxlx = 2.d0*(aa+bb+2.d0*ab-4.d0*bx*ax
     1        -(ax-bx)**2)
         uulyly = 2.d0*(aa+bb+2.d0*ab-4.d0*by*ay
     1        -(ay-by)**2)
         uulzlz = 2.d0*(aa+bb+2.d0*ab-4.d0*bz*az
     1        -(az-bz)**2)
         uulxly = -2.d0*(bx-ax)*(by-ay)
     1        -4.d0*(ax*by+ay*bx)
         uulxlz = -2.d0*(bx-ax)*(bz-az)
     1        -4.d0*(ax*bz+az*bx)
         uulylz = -2.d0*(by-ay)*(bz-az)
     1        -4.d0*(ay*bz+az*by)

         uulxmx = 2.d0*(-aa-ab+(ax+bx)*ax)
         uulymy = 2.d0*(-aa-ab+(ay+by)*ay)
         uulzmz = 2.d0*(-aa-ab+(az+bz)*az)
         uulxmy = 4.d0*ax*by-2.d0*(bx-ax)*ay
         uulxmz = 4.d0*ax*bz-2.d0*(bx-ax)*az
         uulymx = 4.d0*ay*bx-2.d0*(by-ay)*ax
         uulymz = 4.d0*ay*bz-2.d0*(by-ay)*az
         uulzmx = 4.d0*az*bx-2.d0*(bz-az)*ax
         uulzmy = 4.d0*az*by-2.d0*(bz-az)*ay

         uumxmx = 2.d0*(aa-ax*ax)
         uumymy = 2.d0*(aa-ay*ay)
         uumzmz = 2.d0*(aa-az*az)
         uumxmy = -2.d0*ax*ay
         uumxmz = -2.d0*ax*az
         uumymz = -2.d0*ay*az

c vv (all k elements are zero)

         vvlxlx = 2.d0*(cc-cx*cx)
         vvlyly = 2.d0*(cc-cy*cy)
         vvlzlz = 2.d0*(cc-cz*cz)
         vvlxly = -2.d0*cx*cy
         vvlxlz = -2.d0*cx*cz
         vvlylz = -2.d0*cy*cz

         vvlxmx = 2.d0*(cx*(cx+bx)-cc-bc)
         vvlymy = 2.d0*((by+cy)*cy-cc-bc)
         vvlzmz = 2.d0*((bz+cz)*cz-cc-bc)
         vvlxmy = 4.d0*bx*cy+2.d0*cx*(cy-by)
         vvlxmz = 4.d0*bx*cz+2.d0*cx*(cz-bz)
         vvlymx = 4.d0*by*cx+2.d0*cy*(cx-bx)
         vvlymz = 4.d0*by*cz+2.d0*cy*(cz-bz)
         vvlzmx = 4.d0*bz*cx+2.d0*cz*(cx-bx)
         vvlzmy = 4.d0*bz*cy+2.d0*cz*(cy-by)

         vvlxnx = -2.d0*(-bc+bx*cx)
         vvlyny = -2.d0*(-bc+by*cy)
         vvlznz = -2.d0*(-bc+bz*cz)
         vvlxny = 2.d0*(by*cx-2.d0*bx*cy)
         vvlxnz = 2.d0*(bz*cx-2.d0*bx*cz)
         vvlynx = 2.d0*(bx*cy-2.d0*by*cx)
         vvlynz = 2.d0*(bz*cy-2.d0*by*cz)
         vvlznx = 2.d0*(bx*cz-2.d0*bz*cx)
         vvlzny = 2.d0*(by*cz-2.d0*bz*cy)

         vvmxmx = 2.d0*(cc+bb+2.d0*bc-(cx+bx)**2)
         vvmymy = 2.d0*(cc+bb+2.d0*bc-(cy+by)**2)
         vvmzmz = 2.d0*(cc+bb+2.d0*bc-(cz+bz)**2)
         vvmxmy = -2.d0*(cx+bx)*(cy+by)
         vvmxmz = -2.d0*(cx+bx)*(cz+bz)
         vvmymz = -2.d0*(cy+by)*(cz+bz)

         vvmxnx = -2.d0*(bb+bc-bx*(cx+bx))
         vvmyny = -2.d0*(bb+bc-by*(cy+by))
         vvmznz = -2.d0*(bb+bc-bz*(cz+bz))
         vvmxny = 4.d0*bx*cy-2.d0*by*(cx-bx)
         vvmxnz = 4.d0*bx*cz-2.d0*bz*(cx-bx)
         vvmynx = 4.d0*by*cx-2.d0*bx*(cy-by)
         vvmynz = 4.d0*by*cz-2.d0*bz*(cy-by)
         vvmznx = 4.d0*bz*cx-2.d0*bx*(cz-bz)
         vvmzny = 4.d0*bz*cy-2.d0*by*(cz-bz)

         vvnxnx = 2.d0*(bb-bx*bx)
         vvnyny = 2.d0*(bb-by*by)
         vvnznz = 2.d0*(bb-bz*bz)
         vvnxny = -2.d0*bx*by
         vvnxnz = -2.d0*bx*bz
         vvnynz = -2.d0*by*bz
c done with derivatives of scalar products

c Ahead with first cosine derivatives
         cokx=den*(dxkuv-0.5*uv*den2*(dxkuu*vv))
         coky=den*(dykuv-0.5*uv*den2*(dykuu*vv))
         cokz=den*(dzkuv-0.5*uv*den2*(dzkuu*vv))
         colx=den*(dxluv-0.5*uv*den2*(dxluu*vv+dxlvv*uu))
         coly=den*(dyluv-0.5*uv*den2*(dyluu*vv+dylvv*uu))
         colz=den*(dzluv-0.5*uv*den2*(dzluu*vv+dzlvv*uu))
         comx=den*(dxmuv-0.5*uv*den2*(dxmuu*vv+dxmvv*uu))
         comy=den*(dymuv-0.5*uv*den2*(dymuu*vv+dymvv*uu))
         comz=den*(dzmuv-0.5*uv*den2*(dzmuu*vv+dzmvv*uu))
         conx=den*(dxnuv-0.5*uv*den2*(dxnvv*uu))
         cony=den*(dynuv-0.5*uv*den2*(dynvv*uu))
         conz=den*(dznuv-0.5*uv*den2*(dznvv*uu))

c Second derivatives of cosine
         cokxkx=den3*(-dxkuv*vv*dxkuu+
     1        0.75d0*den2*uv*vv*vv*dxkuu*dxkuu-
     2        0.5d0*uv*uukxkx*vv)
         cokyky=den3*(-dykuv*vv*dykuu+
     1        0.75d0*den2*uv*vv*vv*dykuu*dykuu-
     2        0.5d0*uv*uukyky*vv)
         cokzkz=den3*(-dzkuv*vv*dzkuu+
     1        0.75d0*den2*uv*vv*vv*dzkuu*dzkuu-
     2        0.5d0*uv*uukzkz*vv)
         cokxky=den3*(-0.5*(dxkuv*dykuu*vv+dykuv*dxkuu*vv)+
     1        0.75d0*den2*uv*(vv*dxkuu*vv*dykuu)-
     2        0.5d0*uv*uukxky*vv)
         cokxkz=den3*(-0.5*(dxkuv*dzkuu*vv+dzkuv*dxkuu*vv)+
     1        0.75d0*den2*uv*(vv*dxkuu*vv*dzkuu)-
     2        0.5d0*uv*uukxkz*vv)
         cokykz=den3*(-0.5*(dykuv*dzkuu*vv+dzkuv*dykuu*vv)+
     1        0.75d0*den2*uv*(vv*dykuu*vv*dzkuu)-
     2        0.5d0*uv*uukykz*vv)

         conxnx=den3*(-dxnuv*uu*dxnvv+
     1        0.75d0*den2*uv*uu*uu*dxnvv*dxnvv-
     2        0.5d0*uv*uu*vvnxnx)
         conyny=den3*(-dynuv*uu*dynvv+
     1        0.75d0*den2*uv*uu*uu*dynvv*dynvv-
     2        0.5d0*uv*uu*vvnyny)
         conznz=den3*(-dznuv*uu*dznvv+
     1        0.75d0*den2*uv*uu*uu*dznvv*dznvv-
     2        0.5d0*uv*uu*vvnznz)
         conxny=den3*(-0.5d0*(dxnuv*dynvv*uu+dynuv*dxnvv*uu)+
     1        0.75d0*den2*uv*uu*uu*dxnvv*dynvv-
     2        0.5d0*uv*vvnxny*uu)
         conxnz=den3*(-0.5d0*(dxnuv*dznvv*uu+dznuv*dxnvv*uu)+
     1        0.75d0*den2*uv*uu*uu*dxnvv*dznvv-
     2        0.5d0*uv*vvnxnz*uu)
         conynz=den3*(-0.5d0*(dynuv*dznvv*uu+dznuv*dynvv*uu)+
     1        0.75d0*den2*uv*uu*uu*dynvv*dznvv-
     2        0.5d0*uv*vvnynz*uu)

         colxlx=den*uvlxlx+den3*(-dxluv*(dxluu*vv+dxlvv*uu)+
     1        0.75d0*den2*uv*(dxluu*vv+uu*dxlvv)**2-0.5d0*uv*
     2        (uulxlx*vv+uu*vvlxlx+2.d0*dxluu*dxlvv))
         colyly=den*uvlyly+den3*(-dyluv*(dyluu*vv+dylvv*uu)+
     1        0.75d0*den2*uv*(dyluu*vv+uu*dylvv)**2-0.5d0*uv*
     2        (uulyly*vv+uu*vvlyly+2.d0*dyluu*dylvv))
         colzlz=den*uvlzlz+den3*(-dzluv*(dzluu*vv+dzlvv*uu)+
     1        0.75d0*den2*uv*(dzluu*vv+uu*dzlvv)**2-0.5d0*uv*
     2        (uulzlz*vv+uu*vvlzlz+2.d0*dzluu*dzlvv))
         colxly=den*uvlxly+den3*(-0.5d0*dxluv*(dyluu*vv+dylvv*uu)-
     1        0.5d0*dyluv*(dxluu*vv+dxlvv*uu)+0.75d0*uv*den2*(dxluu
     2        *vv+uu*dxlvv)*(dyluu*vv+uu*dylvv)-0.5d0*uv*
     3        (uulxly*vv+uu*vvlxly+dxluu*dylvv+dxlvv*dyluu))
         colxlz=den*uvlxlz+den3*(-0.5d0*dxluv*(dzluu*vv+dzlvv*uu)-
     1        0.5d0*dzluv*(dxluu*vv+dxlvv*uu)+0.75d0*uv*den2*(dxluu
     2        *vv+uu*dxlvv)*(dzluu*vv+uu*dzlvv)-0.5d0*uv*
     3        (uulxlz*vv+uu*vvlxlz+dxluu*dzlvv+dxlvv*dzluu))
         colylz=den*uvlylz+den3*(-0.5d0*dyluv*(dzluu*vv+dzlvv*uu)-
     1        0.5d0*dzluv*(dyluu*vv+dylvv*uu)+0.75d0*uv*den2*(dyluu
     2        *vv+uu*dylvv)*(dzluu*vv+uu*dzlvv)-0.5d0*uv*
     3        (uulylz*vv+uu*vvlylz+dyluu*dzlvv+dylvv*dzluu))

         comxmx=den*uvmxmx+den3*(-dxmuv*(dxmuu*vv+dxmvv*uu)+
     1        0.75d0*den2*uv*(dxmuu*vv+uu*dxmvv)**2-0.5d0*uv*
     2        (uumxmx*vv+uu*vvmxmx+2.d0*dxmuu*dxmvv))
         comymy=den*uvmymy+den3*(-dymuv*(dymuu*vv+dymvv*uu)+
     1        0.75d0*den2*uv*(dymuu*vv+uu*dymvv)**2-0.5d0*uv*
     2        (uumymy*vv+uu*vvmymy+2.d0*dymuu*dymvv))
         comzmz=den*uvmzmz+den3*(-dzmuv*(dzmuu*vv+dzmvv*uu)+
     1        0.75d0*den2*uv*(dzmuu*vv+uu*dzmvv)**2-0.5d0*uv*
     2        (uumzmz*vv+uu*vvmzmz+2.d0*dzmuu*dzmvv))
         comxmy=den*uvmxmy+den3*(-0.5d0*dxmuv*(dymuu*vv+dymvv*uu)-
     1        0.5d0*dymuv*(dxmuu*vv+dxmvv*uu)+0.75d0*uv*den2*(dxmuu
     2        *vv+uu*dxmvv)*(dymuu*vv+uu*dymvv)-0.5d0*uv*
     3        (uumxmy*vv+uu*vvmxmy+dxmuu*dymvv+dxmvv*dymuu))
         comxmz=den*uvmxmz+den3*(-0.5d0*dxmuv*(dzmuu*vv+dzmvv*uu)-
     1        0.5d0*dzmuv*(dxmuu*vv+dxmvv*uu)+0.75d0*uv*den2*(dxmuu
     2        *vv+uu*dxmvv)*(dzmuu*vv+uu*dzmvv)-0.5d0*uv*
     3        (uumxmz*vv+uu*vvmxmz+dxmuu*dzmvv+dxmvv*dzmuu))
         comymz=den*uvmymz+den3*(-0.5d0*dymuv*(dzmuu*vv+dzmvv*uu)-
     1        0.5d0*dzmuv*(dymuu*vv+dymvv*uu)+0.75d0*uv*den2*(dymuu
     2        *vv+uu*dymvv)*(dzmuu*vv+uu*dzmvv)-0.5d0*uv*
     3        (uumymz*vv+uu*vvmymz+dymuu*dzmvv+dymvv*dzmuu))

c off diagonal elements
         cokxlx=uvkxlx*den+den3*(-0.5d0*(dxkuv*(dxluu*vv+dxlvv*uu)+
     1        dxluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dxluu*vv+
     2        uu*dxlvv))-0.5d0*uv*(uukxlx*vv+dxkuu*dxlvv))
         cokyly=uvkyly*den+den3*(-0.5d0*(dykuv*(dyluu*vv+dylvv*uu)+
     1        dyluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dyluu*vv+
     2        uu*dylvv))-0.5d0*uv*(uukyly*vv+dykuu*dylvv))
         cokzlz=uvkzlz*den+den3*(-0.5d0*(dzkuv*(dzluu*vv+dzlvv*uu)+
     1        dzluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzluu*vv+
     2        uu*dzlvv))-0.5d0*uv*(uukzlz*vv+dzkuu*dzlvv))
         cokxly=uvkxly*den+den3*(-0.5d0*(dxkuv*(dyluu*vv+dylvv*uu)+
     1        dyluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dyluu*vv+
     2        uu*dylvv))-0.5d0*uv*(uukxly*vv+dxkuu*dylvv))
         cokxlz=uvkxlz*den+den3*(-0.5d0*(dxkuv*(dzluu*vv+dzlvv*uu)+
     1        dzluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dzluu*vv+
     2        uu*dzlvv))-0.5d0*uv*(uukxlz*vv+dxkuu*dzlvv))
         cokylx=uvkylx*den+den3*(-0.5d0*(dykuv*(dxluu*vv+dxlvv*uu)+
     1        dxluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dxluu*vv+
     2        uu*dxlvv))-0.5d0*uv*(uukylx*vv+dykuu*dxlvv))
         cokylz=uvkylz*den+den3*(-0.5d0*(dykuv*(dzluu*vv+dzlvv*uu)+
     1        dzluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dzluu*vv+
     2        uu*dzlvv))-0.5d0*uv*(uukylz*vv+dykuu*dzlvv))
         cokzlx=uvkzlx*den+den3*(-0.5d0*(dzkuv*(dxluu*vv+dxlvv*uu)+
     1        dxluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dxluu*vv+
     2        uu*dxlvv))-0.5d0*uv*(uukzlx*vv+dzkuu*dxlvv))
         cokzly=uvkzly*den+den3*(-0.5d0*(dzkuv*(dyluu*vv+dylvv*uu)+
     1        dyluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dyluu*vv+
     2        uu*dylvv))-0.5d0*uv*(uukzly*vv+dzkuu*dylvv))
         cokzlz=uvkzlz*den+den3*(-0.5d0*(dzkuv*(dzluu*vv+dzlvv*uu)+
     1        dzluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzluu*vv+
     2        uu*dzlvv))-0.5d0*uv*(uukzlz*vv+dzkuu*dzlvv))

         
         cokxmx=uvkxmx*den+den3*(-0.5d0*(dxkuv*(dxmuu*vv+dxmvv*uu)+
     1        dxmuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dxmuu*vv+
     2        uu*dxmvv))-0.5d0*uv*(uukxmx*vv+dxkuu*dxmvv))
         cokymy=uvkymy*den+den3*(-0.5d0*(dykuv*(dymuu*vv+dymvv*uu)+
     1        dymuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dymuu*vv+
     2        uu*dymvv))-0.5d0*uv*(uukymy*vv+dykuu*dymvv))
         cokzmz=uvkzmz*den+den3*(-0.5d0*(dzkuv*(dzmuu*vv+dzmvv*uu)+
     1        dzmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzmuu*vv+
     2        uu*dzmvv))-0.5d0*uv*(uukzmz*vv+dzkuu*dzmvv))
         cokxmy=uvkxmy*den+den3*(-0.5d0*(dxkuv*(dymuu*vv+dymvv*uu)+
     1        dymuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dymuu*vv+
     2        uu*dymvv))-0.5d0*uv*(uukxmy*vv+dxkuu*dymvv))
         cokxmz=uvkxmz*den+den3*(-0.5d0*(dxkuv*(dzmuu*vv+dzmvv*uu)+
     1        dzmuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dzmuu*vv+
     2        uu*dzmvv))-0.5d0*uv*(uukxmz*vv+dxkuu*dzmvv))
         cokymx=uvkymx*den+den3*(-0.5d0*(dykuv*(dxmuu*vv+dxmvv*uu)+
     1        dxmuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dxmuu*vv+
     2        uu*dxmvv))-0.5d0*uv*(uukymx*vv+dykuu*dxmvv))
         cokymz=uvkymz*den+den3*(-0.5d0*(dykuv*(dzmuu*vv+dzmvv*uu)+
     1        dzmuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dzmuu*vv+
     2        uu*dzmvv))-0.5d0*uv*(uukymz*vv+dykuu*dzmvv))
         cokzmx=uvkzmx*den+den3*(-0.5d0*(dzkuv*(dxmuu*vv+dxmvv*uu)+
     1        dxmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dxmuu*vv+
     2        uu*dxmvv))-0.5d0*uv*(uukzmx*vv+dzkuu*dxmvv))
         cokzmy=uvkzmy*den+den3*(-0.5d0*(dzkuv*(dymuu*vv+dymvv*uu)+
     1        dymuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dymuu*vv+
     2        uu*dymvv))-0.5d0*uv*(uukzmy*vv+dzkuu*dymvv))
         cokzmz=uvkzmz*den+den3*(-0.5d0*(dzkuv*(dzmuu*vv+dzmvv*uu)+
     1        dzmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzmuu*vv+
     2        uu*dzmvv))-0.5d0*uv*(uukzmz*vv+dzkuu*dzmvv))
         
         
         cokxnx=uvkxnx*den+den3*(-0.5d0*(dxkuv*(dxnvv*uu)+
     1        dxnuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2        uu*dxnvv))-0.5d0*uv*(dxkuu*dxnvv))
         cokyny=uvkyny*den+den3*(-0.5d0*(dykuv*(dynvv*uu)+
     1        dynuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2        uu*dynvv))-0.5d0*uv*(dykuu*dynvv))
         cokznz=uvkznz*den+den3*(-0.5d0*(dzkuv*(dznvv*uu)+
     1        dznuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2        uu*dznvv))-0.5d0*uv*(dzkuu*dznvv))
         cokxny=uvkxny*den+den3*(-0.5d0*(dxkuv*(dynvv*uu)+
     1        dynuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2        uu*dynvv))-0.5d0*uv*(dxkuu*dynvv))
         cokxnz=uvkxnz*den+den3*(-0.5d0*(dxkuv*(dznvv*uu)+
     1        dznuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2        uu*dznvv))-0.5d0*uv*(dxkuu*dznvv))
         cokynx=uvkynx*den+den3*(-0.5d0*(dykuv*(dxnvv*uu)+
     1        dxnuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2        uu*dxnvv))-0.5d0*uv*(dykuu*dxnvv))
         cokynz=uvkynz*den+den3*(-0.5d0*(dykuv*(dznvv*uu)+
     1        dznuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2        uu*dznvv))-0.5d0*uv*(dykuu*dznvv))
         cokznx=uvkznx*den+den3*(-0.5d0*(dzkuv*(dxnvv*uu)+
     1        dxnuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2        uu*dxnvv))-0.5d0*uv*(dzkuu*dxnvv))
         cokzny=uvkzny*den+den3*(-0.5d0*(dzkuv*(dynvv*uu)+
     1        dynuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2        uu*dynvv))-0.5d0*uv*(dzkuu*dynvv))

         colxmx=uvlxmx*den+den3*(-0.5d0*(dxluv*(dxmuu*vv+dxmvv*uu)
     1        +dxmuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2        dxluu*vv+uu*dxlvv)*(dxmuu*vv+uu*dxmvv)-
     3        0.5d0*uv*(uulxmx*vv+uu*vvlxmx+dxluu*dxmvv
     4        +dxmuu*dxlvv))
         colymy=uvlymy*den+den3*(-0.5d0*(dyluv*(dymuu*vv+dymvv*uu)
     1        +dymuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2        dyluu*vv+uu*dylvv)*(dymuu*vv+uu*dymvv)-
     3        0.5d0*uv*(uulymy*vv+uu*vvlymy+dyluu*dymvv
     4        +dymuu*dylvv))
         colzmz=uvlzmz*den+den3*(-0.5d0*(dzluv*(dzmuu*vv+dzmvv*uu)
     1        +dzmuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2        dzluu*vv+uu*dzlvv)*(dzmuu*vv+uu*dzmvv)-
     3        0.5d0*uv*(uulzmz*vv+uu*vvlzmz+dzluu*dzmvv
     4        +dzmuu*dzlvv))
         colxmy=uvlxmy*den+den3*(-0.5d0*(dxluv*(dymuu*vv+dymvv*uu)
     1        +dymuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2        dxluu*vv+uu*dxlvv)*(dymuu*vv+uu*dymvv)-
     3        0.5d0*uv*(uulxmy*vv+uu*vvlxmy+dxluu*dymvv
     4        +dymuu*dxlvv))
         colxmz=uvlxmz*den+den3*(-0.5d0*(dxluv*(dzmuu*vv+dzmvv*uu)
     1        +dzmuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2        dxluu*vv+uu*dxlvv)*(dzmuu*vv+uu*dzmvv)-
     3        0.5d0*uv*(uulxmz*vv+uu*vvlxmz+dxluu*dzmvv
     4        +dzmuu*dxlvv))
         colymx=uvlymx*den+den3*(-0.5d0*(dyluv*(dxmuu*vv+dxmvv*uu)
     1        +dxmuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2        dyluu*vv+uu*dylvv)*(dxmuu*vv+uu*dxmvv)-
     3        0.5d0*uv*(uulymx*vv+uu*vvlymx+dyluu*dxmvv
     4        +dxmuu*dylvv))
         colymz=uvlymz*den+den3*(-0.5d0*(dyluv*(dzmuu*vv+dzmvv*uu)
     1        +dzmuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2        dyluu*vv+uu*dylvv)*(dzmuu*vv+uu*dzmvv)-
     3        0.5d0*uv*(uulymz*vv+uu*vvlymz+dyluu*dzmvv
     4        +dzmuu*dylvv))
         colzmx=uvlzmx*den+den3*(-0.5d0*(dzluv*(dxmuu*vv+dxmvv*uu)
     1        +dxmuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2        dzluu*vv+uu*dzlvv)*(dxmuu*vv+uu*dxmvv)-
     3        0.5d0*uv*(uulzmx*vv+uu*vvlzmx+dzluu*dxmvv
     4        +dxmuu*dzlvv))
         colzmy=uvlzmy*den+den3*(-0.5d0*(dzluv*(dymuu*vv+dymvv*uu)
     1        +dymuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2        dzluu*vv+uu*dzlvv)*(dymuu*vv+uu*dymvv)-
     3        0.5d0*uv*(uulzmy*vv+uu*vvlzmy+dzluu*dymvv
     4        +dymuu*dzlvv))
         
         colxnx=uvlxnx*den+den3*(-0.5d0*(dxluv*(dxnvv*uu)+dxnuv*
     1        (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2        *dxlvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlxnx+
     3        dxluu*dxnvv))
         colyny=uvlyny*den+den3*(-0.5d0*(dyluv*(dynvv*uu)+dynuv*
     1        (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2        *dylvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlyny+
     3        dyluu*dynvv))
         colznz=uvlznz*den+den3*(-0.5d0*(dzluv*(dznvv*uu)+dznuv*
     1        (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2        *dzlvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlznz+
     3        dzluu*dznvv))
         colxny=uvlxny*den+den3*(-0.5d0*(dxluv*(dynvv*uu)+dynuv*
     1        (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2        *dxlvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlxny+
     3        dxluu*dynvv))
         colxnz=uvlxnz*den+den3*(-0.5d0*(dxluv*(dznvv*uu)+dznuv*
     1        (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2        *dxlvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlxnz+
     3        dxluu*dznvv))
         colynx=uvlynx*den+den3*(-0.5d0*(dyluv*(dxnvv*uu)+dxnuv*
     1        (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2        *dylvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlynx+
     3        dyluu*dxnvv))
         colynz=uvlynz*den+den3*(-0.5d0*(dyluv*(dznvv*uu)+dznuv*
     1        (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2        *dylvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlynz+
     3        dyluu*dznvv))
         colznx=uvlznx*den+den3*(-0.5d0*(dzluv*(dxnvv*uu)+dxnuv*
     1        (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2        *dzlvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlznx+
     3        dzluu*dxnvv))
         colzny=uvlzny*den+den3*(-0.5d0*(dzluv*(dynvv*uu)+dynuv*
     1        (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2        *dzlvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlzny+
     3        dzluu*dynvv))
         
         
         comxnx=uvmxnx*den+den3*(-0.5d0*(dxmuv*(dxnvv*uu)+dxnuv*
     1        (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2        *dxmvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmxnx+
     3        dxmuu*dxnvv))
         comyny=uvmyny*den+den3*(-0.5d0*(dymuv*(dynvv*uu)+dynuv*
     1        (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2        *dymvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmyny+
     3        dymuu*dynvv))
         comznz=uvmznz*den+den3*(-0.5d0*(dzmuv*(dznvv*uu)+dznuv*
     1        (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2        *dzmvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmznz+
     3        dzmuu*dznvv))
         comxny=uvmxny*den+den3*(-0.5d0*(dxmuv*(dynvv*uu)+dynuv*
     1        (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2        *dxmvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmxny+
     3        dxmuu*dynvv))
         comxnz=uvmxnz*den+den3*(-0.5d0*(dxmuv*(dznvv*uu)+dznuv*
     1        (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2        *dxmvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmxnz+
     3        dxmuu*dznvv))
         comynx=uvmynx*den+den3*(-0.5d0*(dymuv*(dxnvv*uu)+dxnuv*
     1        (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2        *dymvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmynx+
     3        dymuu*dxnvv))
         comynz=uvmynz*den+den3*(-0.5d0*(dymuv*(dznvv*uu)+dznuv*
     1        (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2        *dymvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmynz+
     3        dymuu*dznvv))
         comznx=uvmznx*den+den3*(-0.5d0*(dzmuv*(dxnvv*uu)+dxnuv*
     1        (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2        *dzmvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmznx+
     3        dzmuu*dxnvv))
         comzny=uvmzny*den+den3*(-0.5d0*(dzmuv*(dynvv*uu)+dynuv*
     1        (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2        *dzmvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmzny+
     3        dzmuu*dynvv))
         
c     Done with off diagonal elements :-)
         
C     DIAGONAL ELEMENTS OF THE SECOND DERIVATIVE MATRIX
         kk = 6*(k-1)+1
         ll = 6*(l-1)+1
         mm = 6*(m-1)+1
         nn = 6*(n-1)+1
c@
         diag(kk) = diag(kk) + scalar1*cokx*cokx + scalar2*cokxkx
c@
         diag(kk+1) = diag(kk+1) + scalar1*cokx*coky + scalar2*cokxky
         diag(kk+2) = diag(kk+2) + scalar1*cokx*cokz + scalar2*cokxkz
         diag(kk+3) = diag(kk+3) + scalar1*coky*coky + scalar2*cokyky
         diag(kk+4) = diag(kk+4) + scalar1*coky*cokz + scalar2*cokykz
         diag(kk+5) = diag(kk+5) + scalar1*cokz*cokz + scalar2*cokzkz

         diag(ll) = diag(ll) + scalar1*colx*colx + scalar2*colxlx
         diag(ll+1) = diag(ll+1) + scalar1*colx*coly + scalar2*colxly
         diag(ll+2) = diag(ll+2) + scalar1*colx*colz + scalar2*colxlz
         diag(ll+3) = diag(ll+3) + scalar1*coly*coly + scalar2*colyly
         diag(ll+4) = diag(ll+4) + scalar1*coly*colz + scalar2*colylz
         diag(ll+5) = diag(ll+5) + scalar1*colz*colz + scalar2*colzlz

         diag(mm) = diag(mm) + scalar1*comx*comx + scalar2*comxmx
         diag(mm+1) = diag(mm+1) + scalar1*comx*comy + scalar2*comxmy
         diag(mm+2) = diag(mm+2) + scalar1*comx*comz + scalar2*comxmz
         diag(mm+3) = diag(mm+3) + scalar1*comy*comy + scalar2*comymy
         diag(mm+4) = diag(mm+4) + scalar1*comy*comz + scalar2*comymz
         diag(mm+5) = diag(mm+5) + scalar1*comz*comz + scalar2*comzmz

         diag(nn) = diag(nn) + scalar1*conx*conx + scalar2*conxnx
         diag(nn+1) = diag(nn+1) + scalar1*conx*cony + scalar2*conxny
         diag(nn+2) = diag(nn+2) + scalar1*conx*conz + scalar2*conxnz
         diag(nn+3) = diag(nn+3) + scalar1*cony*cony + scalar2*conyny
         diag(nn+4) = diag(nn+4) + scalar1*cony*conz + scalar2*conynz
         diag(nn+5) = diag(nn+5) + scalar1*conz*conz + scalar2*conznz

C OFF DIAGONAL SECOND DERIVATIVE MATRIX ELEMENTS

         iiphi = 54*(iphi-ifirst)+1

c kl
         offdiag(iiphi) = scalar1*cokx*colx + scalar2*cokxlx
         offdiag(iiphi+1) = scalar1*coky*colx + scalar2*cokylx
         offdiag(iiphi+2) = scalar1*cokz*colx + scalar2*cokzlx
         offdiag(iiphi+3) = scalar1*cokx*coly + scalar2*cokxly
         offdiag(iiphi+4) = scalar1*coky*coly + scalar2*cokyly
         offdiag(iiphi+5) = scalar1*cokz*coly + scalar2*cokzly
         offdiag(iiphi+6) = scalar1*cokx*colz + scalar2*cokxlz
         offdiag(iiphi+7) = scalar1*coky*colz + scalar2*cokylz
         offdiag(iiphi+8) = scalar1*cokz*colz + scalar2*cokzlz
c km
         offdiag(iiphi+9) = scalar1*cokx*comx + scalar2*cokxmx
         offdiag(iiphi+10) = scalar1*coky*comx + scalar2*cokymx
         offdiag(iiphi+11) = scalar1*cokz*comx + scalar2*cokzmx
         offdiag(iiphi+12) = scalar1*cokx*comy + scalar2*cokxmy
         offdiag(iiphi+13) = scalar1*coky*comy + scalar2*cokymy
         offdiag(iiphi+14) = scalar1*cokz*comy + scalar2*cokzmy
         offdiag(iiphi+15) = scalar1*cokx*comz + scalar2*cokxmz
         offdiag(iiphi+16) = scalar1*coky*comz + scalar2*cokymz
         offdiag(iiphi+17) = scalar1*cokz*comz + scalar2*cokzmz

c kn
         offdiag(iiphi+18) = scalar1*cokx*conx + scalar2*cokxnx
         offdiag(iiphi+19) = scalar1*coky*conx + scalar2*cokynx
         offdiag(iiphi+20) = scalar1*cokz*conx + scalar2*cokznx
         offdiag(iiphi+21) = scalar1*cokx*cony + scalar2*cokxny
         offdiag(iiphi+22) = scalar1*coky*cony + scalar2*cokyny
         offdiag(iiphi+23) = scalar1*cokz*cony + scalar2*cokzny
         offdiag(iiphi+24) = scalar1*cokx*conz + scalar2*cokxnz
         offdiag(iiphi+25) = scalar1*coky*conz + scalar2*cokynz
         offdiag(iiphi+26) = scalar1*cokz*conz + scalar2*cokznz

c lm
         offdiag(iiphi+27) = scalar1*colx*comx + scalar2*colxmx
         offdiag(iiphi+28) = scalar1*coly*comx + scalar2*colymx
         offdiag(iiphi+29) = scalar1*colz*comx + scalar2*colzmx
         offdiag(iiphi+30) = scalar1*colx*comy + scalar2*colxmy
         offdiag(iiphi+31) = scalar1*coly*comy + scalar2*colymy
         offdiag(iiphi+32) = scalar1*colz*comy + scalar2*colzmy
         offdiag(iiphi+33) = scalar1*colx*comz + scalar2*colxmz
         offdiag(iiphi+34) = scalar1*coly*comz + scalar2*colymz
         offdiag(iiphi+35) = scalar1*colz*comz + scalar2*colzmz
c ln
         offdiag(iiphi+36) = scalar1*colx*conx + scalar2*colxnx
         offdiag(iiphi+37) = scalar1*coly*conx + scalar2*colynx
         offdiag(iiphi+38) = scalar1*colz*conx + scalar2*colznx
         offdiag(iiphi+39) = scalar1*colx*cony + scalar2*colxny
         offdiag(iiphi+40) = scalar1*coly*cony + scalar2*colyny
         offdiag(iiphi+41) = scalar1*colz*cony + scalar2*colzny
         offdiag(iiphi+42) = scalar1*colx*conz + scalar2*colxnz
         offdiag(iiphi+43) = scalar1*coly*conz + scalar2*colynz
         offdiag(iiphi+44) = scalar1*colz*conz + scalar2*colznz
c mn
         offdiag(iiphi+45) = scalar1*comx*conx + scalar2*comxnx
         offdiag(iiphi+46) = scalar1*comy*conx + scalar2*comynx
         offdiag(iiphi+47) = scalar1*comz*conx + scalar2*comznx
         offdiag(iiphi+48) = scalar1*comx*cony + scalar2*comxny
         offdiag(iiphi+49) = scalar1*comy*cony + scalar2*comyny
         offdiag(iiphi+50) = scalar1*comz*cony + scalar2*comzny
         offdiag(iiphi+51) = scalar1*comx*conz + scalar2*comxnz
         offdiag(iiphi+52) = scalar1*comy*conz + scalar2*comynz
         offdiag(iiphi+53) = scalar1*comz*conz + scalar2*comznz

 100  continue

c All done (?!)

      return
      end










