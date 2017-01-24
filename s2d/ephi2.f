        subroutine etors2()

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SCNDRV.BLOCK'

c       calculate dihedral energies and forces according
c       to tricks with scalar products
c       of T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 

c       V = sum [ k(n)(1 + cos(n*phi + delta)]   (1=< n <= 3)
c       delta is either 0 or pi, so the formula become
c       less complicated upon expansion.



        integer iphi,iiphi,k,l,m,n,kk,ll,mm,nn
        double precision ax,ay,az,bx,by,bz,cx,cy,cz
        double precision ab,bc,ac,aa,bb,cc,uu,vv,uv
        double precision den,den2,den3,co,dcos2,d2cos

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
     1          uvkylz,uvkzlx,uvkzly
        double precision uvkxmx,uvkymy,uvkzmz,uvkxmy,uvkxmz,uvkymx,
     1          uvkymz,uvkzmx,uvkzmy
        double precision uvkxnx,uvkyny,uvkznz,uvkxny,uvkxnz,uvkynx,
     1          uvkynz,uvkznx,uvkzny
        double precision uvlxmx,uvlymy,uvlzmz,uvlxmy,uvlxmz,uvlymx,
     1          uvlymz,uvlzmx,uvlzmy
        double precision uvlxnx,uvlyny,uvlznz,uvlxny,uvlxnz,uvlynx,
     1          uvlynz,uvlznx,uvlzny
        double precision uvmxnx,uvmyny,uvmznz,uvmxny,uvmxnz,uvmynx,
     1          uvmynz,uvmznx,uvmzny


        double precision uukxkx,uukyky,uukzkz,uukxky,uukxkz,uukykz
        double precision uukxlx,uukyly,uukzlz,uukxly,uukxlz,uukylx,
     1          uukylz,uukzlx,uukzly
        double precision uukxmx,uukymy,uukzmz,uukxmy,uukxmz,uukymx,
     1          uukymz,uukzmx,uukzmy
        double precision uulxlx,uulyly,uulzlz,uulxly,uulxlz,uulylz
        double precision uulxmx,uulymy,uulzmz,uulxmy,uulxmz,uulymx,
     1          uulymz,uulzmx,uulzmy
        double precision uumxmx,uumymy,uumzmz,uumxmy,uumxmz,uumymz

        double precision vvlxlx,vvlyly,vvlzlz,vvlxly,vvlxlz,vvlylz
        double precision vvlxmx,vvlymy,vvlzmz,vvlxmy,vvlxmz,vvlymx,
     1          vvlymz,vvlzmx,vvlzmy
        double precision vvlxnx,vvlyny,vvlznz,vvlxny,vvlxnz,vvlynx,
     1          vvlynz,vvlznx,vvlzny
        double precision vvmxmx,vvmymy,vvmzmz,vvmxmy,vvmxmz,vvmymz
        double precision vvmxnx,vvmyny,vvmznz,vvmxny,vvmxnz,vvmynx,
     1          vvmynz,vvmznx,vvmzny
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

        do 100 iphi = 1,ntors


c       call the torsion atoms
                        k = itor1(iphi)
                        l = itor2(iphi)
                        m = itor3(iphi)
                        n = itor4(iphi)


c       calculate distances
                        ax = coor(1,l) - coor(1,k)
                        ay = coor(2,l) - coor(2,k)
                        az = coor(3,l) - coor(3,k)

                        bx = coor(1,m) - coor(1,l)
                        by = coor(2,m) - coor(2,l)
                        bz = coor(3,m) - coor(3,l)

                        cx = coor(1,n) - coor(1,m)
                        cy = coor(2,n) - coor(2,m)
                        cz = coor(3,n) - coor(3,m)

c       calculate dot products
                        ab = ax*bx + ay*by + az*bz
                        bc = bx*cx + by*cy + bz*cz
                        ac = ax*cx + ay*cy + az*cz
                        aa = ax*ax + ay*ay + az*az
                        bb = bx*bx + by*by + bz*bz
                        cc = cx*cx + cy*cy + cz*cz


c       calculate cos(phi)      
                        uu  = (aa * bb) - (ab*ab)
                        vv  = (bb * cc ) - (bc * bc)
                        uv  = (ab * bc) - (ac * bb)
                        den  = 1.d0/dsqrt(uu*vv) 
                        den2 = den*den
                        den3 = den2*den
                        co   = uv * den
c@
C@                      write(*,*)' uu vv uv den co'
C@                      write(*,*)uu,vv,uv,den,co

c       calculate multplicative constants
        if(.not.arith)then
                        d2cos = phase3(iphi) * 12.d0*ktors3(iphi)*co*co
     1                        + phase2(iphi) * 4.d0 * ktors2(iphi)*co
     2                        + phase1(iphi) * ktors1(iphi) 
     2                        - phase3(iphi) * 3.d0*ktors3(iphi)
                        dcos2 = phase3(iphi) * 24.d0*ktors3(iphi)*co+
     1                          phase2(iphi) * 4.d0*ktors2(iphi)
         else
           if(period(iphi).eq.1) then
            d2cos = ktors2(iphi)*phase1(iphi)
            dcos2 = 0.0d0
           else if(period(iphi).eq.2) then
            d2cos = 4.0d0*ktors2(iphi)*phase1(iphi)*co
            dcos2 = 4.0d0*ktors2(iphi)*phase1(iphi)
           else if(period(iphi).eq.3) then
            d2cos = ktors2(iphi)*phase1(iphi)*(12.0d0*co*co-3.0d0)
            dcos2 = ktors2(iphi)*phase1(iphi)*(24.0d0*co)
           else if(period(iphi).eq.4) then
            d2cos = 16.0d0*ktors2(iphi)*phase1(iphi)*co*
     &              (2.0d0*co*co-1.0d0)
            dcos2 = 16.0d0*ktors2(iphi)*phase1(iphi)*
     &              (6.0d0*co*co-1.0d0)
           endif
          endif


c@
C@                      write(*,*)' dcos2 d2cos ',dcos2,d2cos
c       calculate derivatives of scalar products 
c uv
                        dxkuv = -bx*bc + cx*bb
                        dykuv = -by*bc + cy*bb
                        dzkuv = -bz*bc + cz*bb

                        dxluv = (bx-ax)*bc-cx*(ab+bb)+2.d0*bx*ac
                        dyluv = (by-ay)*bc-cy*(ab+bb)+2.d0*by*ac
                        dzluv = (bz-az)*bc-cz*(ab+bb)+2.d0*bz*ac

c@
C@                      write(*,*)' dxluv dyluv dzluv '
C@                      write(*,*)dxluv,dyluv,dzluv
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
     1                          ax*cx+bb+2.d0*(ac-cx*bx-ax*bx)
                        uvlymy = bc+(by-ay)*(cy-by)+ab-
     1                          ay*cy+bb+2.d0*(ac-cy*by-ay*by)
                        uvlzmz = bc+(bz-az)*(cz-bz)+ab-
     1                          az*cz+bb+2.d0*(ac-cz*bz-az*bz)
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
     1                          -(ax-bx)**2)
                        uulyly = 2.d0*(aa+bb+2.d0*ab-4.d0*by*ay
     1                          -(ay-by)**2)
                        uulzlz = 2.d0*(aa+bb+2.d0*ab-4.d0*bz*az
     1                          -(az-bz)**2)
                        uulxly = -2.d0*(bx-ax)*(by-ay)
     1                          -4.d0*(ax*by+ay*bx)
                        uulxlz = -2.d0*(bx-ax)*(bz-az)
     1                          -4.d0*(ax*bz+az*bx)
                        uulylz = -2.d0*(by-ay)*(bz-az)
     1                          -4.d0*(ay*bz+az*by)

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
     1                  0.75d0*den2*uv*vv*vv*dxkuu*dxkuu-
     2                  0.5d0*uv*uukxkx*vv)
                cokyky=den3*(-dykuv*vv*dykuu+
     1                  0.75d0*den2*uv*vv*vv*dykuu*dykuu-
     2                  0.5d0*uv*uukyky*vv)
C@              write(*,*)' coky cokyky = ',coky,cokyky
C@              write(*,*)' dykuu dykvv uukyky ',dykuu,dykuv,uukyky
                cokzkz=den3*(-dzkuv*vv*dzkuu+
     1                  0.75d0*den2*uv*vv*vv*dzkuu*dzkuu-
     2                  0.5d0*uv*uukzkz*vv)
                cokxky=den3*(-0.5*(dxkuv*dykuu*vv+dykuv*dxkuu*vv)+
     1                  0.75d0*den2*uv*(vv*dxkuu*vv*dykuu)-
     2                  0.5d0*uv*uukxky*vv)
                cokxkz=den3*(-0.5*(dxkuv*dzkuu*vv+dzkuv*dxkuu*vv)+
     1                  0.75d0*den2*uv*(vv*dxkuu*vv*dzkuu)-
     2                  0.5d0*uv*uukxkz*vv)
                cokykz=den3*(-0.5*(dykuv*dzkuu*vv+dzkuv*dykuu*vv)+
     1                  0.75d0*den2*uv*(vv*dykuu*vv*dzkuu)-
     2                  0.5d0*uv*uukykz*vv)

                conxnx=den3*(-dxnuv*uu*dxnvv+
     1                  0.75d0*den2*uv*uu*uu*dxnvv*dxnvv-
     2                  0.5d0*uv*uu*vvnxnx)
                conyny=den3*(-dynuv*uu*dynvv+
     1                  0.75d0*den2*uv*uu*uu*dynvv*dynvv-
     2                  0.5d0*uv*uu*vvnyny)
                conznz=den3*(-dznuv*uu*dznvv+
     1                  0.75d0*den2*uv*uu*uu*dznvv*dznvv-
     2                  0.5d0*uv*uu*vvnznz)
                conxny=den3*(-0.5d0*(dxnuv*dynvv*uu+dynuv*dxnvv*uu)+
     1                  0.75d0*den2*uv*uu*uu*dxnvv*dynvv-
     2                  0.5d0*uv*vvnxny*uu)
                conxnz=den3*(-0.5d0*(dxnuv*dznvv*uu+dznuv*dxnvv*uu)+
     1                  0.75d0*den2*uv*uu*uu*dxnvv*dznvv-
     2                  0.5d0*uv*vvnxnz*uu)
                conynz=den3*(-0.5d0*(dynuv*dznvv*uu+dznuv*dynvv*uu)+
     1                  0.75d0*den2*uv*uu*uu*dynvv*dznvv-
     2                  0.5d0*uv*vvnynz*uu)

                colxlx=den*uvlxlx+den3*(-dxluv*(dxluu*vv+dxlvv*uu)+
     1                  0.75d0*den2*uv*(dxluu*vv+uu*dxlvv)**2-0.5d0*uv*
     2                  (uulxlx*vv+uu*vvlxlx+2.d0*dxluu*dxlvv))
C@              write(*,*)'colxlx ',colxlx
C@              write(*,*)' uvlxlx dxluv dxluu dxlvv uulxlx vvlxlx '
C@              write(*,*)uvlxlx,dxluv,dxluu,dxlvv,uulxlx,vvlxlx
                colyly=den*uvlyly+den3*(-dyluv*(dyluu*vv+dylvv*uu)+
     1                  0.75d0*den2*uv*(dyluu*vv+uu*dylvv)**2-0.5d0*uv*
     2                  (uulyly*vv+uu*vvlyly+2.d0*dyluu*dylvv))
                colzlz=den*uvlzlz+den3*(-dzluv*(dzluu*vv+dzlvv*uu)+
     1                  0.75d0*den2*uv*(dzluu*vv+uu*dzlvv)**2-0.5d0*uv*
     2                  (uulzlz*vv+uu*vvlzlz+2.d0*dzluu*dzlvv))
        colxly=den*uvlxly+den3*(-0.5d0*dxluv*(dyluu*vv+dylvv*uu)-
     1          0.5d0*dyluv*(dxluu*vv+dxlvv*uu)+0.75d0*uv*den2*(dxluu
     2                  *vv+uu*dxlvv)*(dyluu*vv+uu*dylvv)-0.5d0*uv*
     3                  (uulxly*vv+uu*vvlxly+dxluu*dylvv+dxlvv*dyluu))
        colxlz=den*uvlxlz+den3*(-0.5d0*dxluv*(dzluu*vv+dzlvv*uu)-
     1          0.5d0*dzluv*(dxluu*vv+dxlvv*uu)+0.75d0*uv*den2*(dxluu
     2                  *vv+uu*dxlvv)*(dzluu*vv+uu*dzlvv)-0.5d0*uv*
     3                  (uulxlz*vv+uu*vvlxlz+dxluu*dzlvv+dxlvv*dzluu))
C@      write(*,*)' uvlxlz dxluv dzluu dzlvv ',uvlxlz,dxluv,dzluu,dzlvv
C@      write(*,*)' colxlz colx colz ',colxlz,colx,colz
        colylz=den*uvlylz+den3*(-0.5d0*dyluv*(dzluu*vv+dzlvv*uu)-
     1          0.5d0*dzluv*(dyluu*vv+dylvv*uu)+0.75d0*uv*den2*(dyluu
     2                  *vv+uu*dylvv)*(dzluu*vv+uu*dzlvv)-0.5d0*uv*
     3                  (uulylz*vv+uu*vvlylz+dyluu*dzlvv+dylvv*dzluu))

                comxmx=den*uvmxmx+den3*(-dxmuv*(dxmuu*vv+dxmvv*uu)+
     1                  0.75d0*den2*uv*(dxmuu*vv+uu*dxmvv)**2-0.5d0*uv*
     2                  (uumxmx*vv+uu*vvmxmx+2.d0*dxmuu*dxmvv))
C@              write(*,*)'comxmx ',comxmx
C@              write(*,*)' uvmxmx dxmuv dxmuu dxmvv uumxmx vvmxmx '
C@              write(*,*)uvmxmx,dxmuv,dxmuu,dxmvv,uumxmx,vvmxmx
                comymy=den*uvmymy+den3*(-dymuv*(dymuu*vv+dymvv*uu)+
     1                  0.75d0*den2*uv*(dymuu*vv+uu*dymvv)**2-0.5d0*uv*
     2                  (uumymy*vv+uu*vvmymy+2.d0*dymuu*dymvv))
                comzmz=den*uvmzmz+den3*(-dzmuv*(dzmuu*vv+dzmvv*uu)+
     1                  0.75d0*den2*uv*(dzmuu*vv+uu*dzmvv)**2-0.5d0*uv*
     2                  (uumzmz*vv+uu*vvmzmz+2.d0*dzmuu*dzmvv))
        comxmy=den*uvmxmy+den3*(-0.5d0*dxmuv*(dymuu*vv+dymvv*uu)-
     1          0.5d0*dymuv*(dxmuu*vv+dxmvv*uu)+0.75d0*uv*den2*(dxmuu
     2                  *vv+uu*dxmvv)*(dymuu*vv+uu*dymvv)-0.5d0*uv*
     3                  (uumxmy*vv+uu*vvmxmy+dxmuu*dymvv+dxmvv*dymuu))
        comxmz=den*uvmxmz+den3*(-0.5d0*dxmuv*(dzmuu*vv+dzmvv*uu)-
     1          0.5d0*dzmuv*(dxmuu*vv+dxmvv*uu)+0.75d0*uv*den2*(dxmuu
     2                  *vv+uu*dxmvv)*(dzmuu*vv+uu*dzmvv)-0.5d0*uv*
     3                  (uumxmz*vv+uu*vvmxmz+dxmuu*dzmvv+dxmvv*dzmuu))
        comymz=den*uvmymz+den3*(-0.5d0*dymuv*(dzmuu*vv+dzmvv*uu)-
     1          0.5d0*dzmuv*(dymuu*vv+dymvv*uu)+0.75d0*uv*den2*(dymuu
     2                  *vv+uu*dymvv)*(dzmuu*vv+uu*dzmvv)-0.5d0*uv*
     3                  (uumymz*vv+uu*vvmymz+dymuu*dzmvv+dymvv*dzmuu))

c off diaggonal elements
        cokxlx=uvkxlx*den+den3*(-0.5d0*(dxkuv*(dxluu*vv+dxlvv*uu)+
     1          dxluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dxluu*vv+
     2                  uu*dxlvv))-0.5d0*uv*(uukxlx*vv+dxkuu*dxlvv))
        cokyly=uvkyly*den+den3*(-0.5d0*(dykuv*(dyluu*vv+dylvv*uu)+
     1          dyluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dyluu*vv+
     2                  uu*dylvv))-0.5d0*uv*(uukyly*vv+dykuu*dylvv))
        cokzlz=uvkzlz*den+den3*(-0.5d0*(dzkuv*(dzluu*vv+dzlvv*uu)+
     1          dzluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzluu*vv+
     2                  uu*dzlvv))-0.5d0*uv*(uukzlz*vv+dzkuu*dzlvv))
        cokxly=uvkxly*den+den3*(-0.5d0*(dxkuv*(dyluu*vv+dylvv*uu)+
     1          dyluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dyluu*vv+
     2                  uu*dylvv))-0.5d0*uv*(uukxly*vv+dxkuu*dylvv))
        cokxlz=uvkxlz*den+den3*(-0.5d0*(dxkuv*(dzluu*vv+dzlvv*uu)+
     1          dzluv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dzluu*vv+
     2                  uu*dzlvv))-0.5d0*uv*(uukxlz*vv+dxkuu*dzlvv))
        cokylx=uvkylx*den+den3*(-0.5d0*(dykuv*(dxluu*vv+dxlvv*uu)+
     1          dxluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dxluu*vv+
     2                  uu*dxlvv))-0.5d0*uv*(uukylx*vv+dykuu*dxlvv))
        cokylz=uvkylz*den+den3*(-0.5d0*(dykuv*(dzluu*vv+dzlvv*uu)+
     1          dzluv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dzluu*vv+
     2                  uu*dzlvv))-0.5d0*uv*(uukylz*vv+dykuu*dzlvv))
        cokzlx=uvkzlx*den+den3*(-0.5d0*(dzkuv*(dxluu*vv+dxlvv*uu)+
     1          dxluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dxluu*vv+
     2                  uu*dxlvv))-0.5d0*uv*(uukzlx*vv+dzkuu*dxlvv))
        cokzly=uvkzly*den+den3*(-0.5d0*(dzkuv*(dyluu*vv+dylvv*uu)+
     1          dyluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dyluu*vv+
     2                  uu*dylvv))-0.5d0*uv*(uukzly*vv+dzkuu*dylvv))
        cokzlz=uvkzlz*den+den3*(-0.5d0*(dzkuv*(dzluu*vv+dzlvv*uu)+
     1          dzluv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzluu*vv+
     2                  uu*dzlvv))-0.5d0*uv*(uukzlz*vv+dzkuu*dzlvv))

                
        cokxmx=uvkxmx*den+den3*(-0.5d0*(dxkuv*(dxmuu*vv+dxmvv*uu)+
     1          dxmuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dxmuu*vv+
     2                  uu*dxmvv))-0.5d0*uv*(uukxmx*vv+dxkuu*dxmvv))
        cokymy=uvkymy*den+den3*(-0.5d0*(dykuv*(dymuu*vv+dymvv*uu)+
     1          dymuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dymuu*vv+
     2                  uu*dymvv))-0.5d0*uv*(uukymy*vv+dykuu*dymvv))
        cokzmz=uvkzmz*den+den3*(-0.5d0*(dzkuv*(dzmuu*vv+dzmvv*uu)+
     1          dzmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzmuu*vv+
     2                  uu*dzmvv))-0.5d0*uv*(uukzmz*vv+dzkuu*dzmvv))
        cokxmy=uvkxmy*den+den3*(-0.5d0*(dxkuv*(dymuu*vv+dymvv*uu)+
     1          dymuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dymuu*vv+
     2                  uu*dymvv))-0.5d0*uv*(uukxmy*vv+dxkuu*dymvv))
        cokxmz=uvkxmz*den+den3*(-0.5d0*(dxkuv*(dzmuu*vv+dzmvv*uu)+
     1          dzmuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(dzmuu*vv+
     2                  uu*dzmvv))-0.5d0*uv*(uukxmz*vv+dxkuu*dzmvv))
        cokymx=uvkymx*den+den3*(-0.5d0*(dykuv*(dxmuu*vv+dxmvv*uu)+
     1          dxmuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dxmuu*vv+
     2                  uu*dxmvv))-0.5d0*uv*(uukymx*vv+dykuu*dxmvv))
        cokymz=uvkymz*den+den3*(-0.5d0*(dykuv*(dzmuu*vv+dzmvv*uu)+
     1          dzmuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(dzmuu*vv+
     2                  uu*dzmvv))-0.5d0*uv*(uukymz*vv+dykuu*dzmvv))
        cokzmx=uvkzmx*den+den3*(-0.5d0*(dzkuv*(dxmuu*vv+dxmvv*uu)+
     1          dxmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dxmuu*vv+
     2                  uu*dxmvv))-0.5d0*uv*(uukzmx*vv+dzkuu*dxmvv))
        cokzmy=uvkzmy*den+den3*(-0.5d0*(dzkuv*(dymuu*vv+dymvv*uu)+
     1          dymuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dymuu*vv+
     2                  uu*dymvv))-0.5d0*uv*(uukzmy*vv+dzkuu*dymvv))
        cokzmz=uvkzmz*den+den3*(-0.5d0*(dzkuv*(dzmuu*vv+dzmvv*uu)+
     1          dzmuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(dzmuu*vv+
     2                  uu*dzmvv))-0.5d0*uv*(uukzmz*vv+dzkuu*dzmvv))
                
                
        cokxnx=uvkxnx*den+den3*(-0.5d0*(dxkuv*(dxnvv*uu)+
     1          dxnuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2                  uu*dxnvv))-0.5d0*uv*(dxkuu*dxnvv))
        cokyny=uvkyny*den+den3*(-0.5d0*(dykuv*(dynvv*uu)+
     1          dynuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2                  uu*dynvv))-0.5d0*uv*(dykuu*dynvv))
        cokznz=uvkznz*den+den3*(-0.5d0*(dzkuv*(dznvv*uu)+
     1          dznuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2                  uu*dznvv))-0.5d0*uv*(dzkuu*dznvv))
        cokxny=uvkxny*den+den3*(-0.5d0*(dxkuv*(dynvv*uu)+
     1          dynuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2                  uu*dynvv))-0.5d0*uv*(dxkuu*dynvv))
        cokxnz=uvkxnz*den+den3*(-0.5d0*(dxkuv*(dznvv*uu)+
     1          dznuv*dxkuu*vv)+0.75d0*uv*den2*(dxkuu*vv*(
     2                  uu*dznvv))-0.5d0*uv*(dxkuu*dznvv))
        cokynx=uvkynx*den+den3*(-0.5d0*(dykuv*(dxnvv*uu)+
     1          dxnuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2                  uu*dxnvv))-0.5d0*uv*(dykuu*dxnvv))
        cokynz=uvkynz*den+den3*(-0.5d0*(dykuv*(dznvv*uu)+
     1          dznuv*dykuu*vv)+0.75d0*uv*den2*(dykuu*vv*(
     2                  uu*dznvv))-0.5d0*uv*(dykuu*dznvv))
        cokznx=uvkznx*den+den3*(-0.5d0*(dzkuv*(dxnvv*uu)+
     1          dxnuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2                  uu*dxnvv))-0.5d0*uv*(dzkuu*dxnvv))
        cokzny=uvkzny*den+den3*(-0.5d0*(dzkuv*(dynvv*uu)+
     1          dynuv*dzkuu*vv)+0.75d0*uv*den2*(dzkuu*vv*(
     2                  uu*dynvv))-0.5d0*uv*(dzkuu*dynvv))

        colxmx=uvlxmx*den+den3*(-0.5d0*(dxluv*(dxmuu*vv+dxmvv*uu)
     1          +dxmuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2                  dxluu*vv+uu*dxlvv)*(dxmuu*vv+uu*dxmvv)-
     3                  0.5d0*uv*(uulxmx*vv+uu*vvlxmx+dxluu*dxmvv
     4                  +dxmuu*dxlvv))
        colymy=uvlymy*den+den3*(-0.5d0*(dyluv*(dymuu*vv+dymvv*uu)
     1          +dymuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2                  dyluu*vv+uu*dylvv)*(dymuu*vv+uu*dymvv)-
     3                  0.5d0*uv*(uulymy*vv+uu*vvlymy+dyluu*dymvv
     4                  +dymuu*dylvv))
        colzmz=uvlzmz*den+den3*(-0.5d0*(dzluv*(dzmuu*vv+dzmvv*uu)
     1          +dzmuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2                  dzluu*vv+uu*dzlvv)*(dzmuu*vv+uu*dzmvv)-
     3                  0.5d0*uv*(uulzmz*vv+uu*vvlzmz+dzluu*dzmvv
     4                  +dzmuu*dzlvv))
        colxmy=uvlxmy*den+den3*(-0.5d0*(dxluv*(dymuu*vv+dymvv*uu)
     1          +dymuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2                  dxluu*vv+uu*dxlvv)*(dymuu*vv+uu*dymvv)-
     3                  0.5d0*uv*(uulxmy*vv+uu*vvlxmy+dxluu*dymvv
     4                  +dymuu*dxlvv))
        colxmz=uvlxmz*den+den3*(-0.5d0*(dxluv*(dzmuu*vv+dzmvv*uu)
     1          +dzmuv*(dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(
     2                  dxluu*vv+uu*dxlvv)*(dzmuu*vv+uu*dzmvv)-
     3                  0.5d0*uv*(uulxmz*vv+uu*vvlxmz+dxluu*dzmvv
     4                  +dzmuu*dxlvv))
        colymx=uvlymx*den+den3*(-0.5d0*(dyluv*(dxmuu*vv+dxmvv*uu)
     1          +dxmuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2                  dyluu*vv+uu*dylvv)*(dxmuu*vv+uu*dxmvv)-
     3                  0.5d0*uv*(uulymx*vv+uu*vvlymx+dyluu*dxmvv
     4                  +dxmuu*dylvv))
        colymz=uvlymz*den+den3*(-0.5d0*(dyluv*(dzmuu*vv+dzmvv*uu)
     1          +dzmuv*(dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(
     2                  dyluu*vv+uu*dylvv)*(dzmuu*vv+uu*dzmvv)-
     3                  0.5d0*uv*(uulymz*vv+uu*vvlymz+dyluu*dzmvv
     4                  +dzmuu*dylvv))
        colzmx=uvlzmx*den+den3*(-0.5d0*(dzluv*(dxmuu*vv+dxmvv*uu)
     1          +dxmuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2                  dzluu*vv+uu*dzlvv)*(dxmuu*vv+uu*dxmvv)-
     3                  0.5d0*uv*(uulzmx*vv+uu*vvlzmx+dzluu*dxmvv
     4                  +dxmuu*dzlvv))
        colzmy=uvlzmy*den+den3*(-0.5d0*(dzluv*(dymuu*vv+dymvv*uu)
     1          +dymuv*(dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(
     2                  dzluu*vv+uu*dzlvv)*(dymuu*vv+uu*dymvv)-
     3                  0.5d0*uv*(uulzmy*vv+uu*vvlzmy+dzluu*dymvv
     4                  +dymuu*dzlvv))

        colxnx=uvlxnx*den+den3*(-0.5d0*(dxluv*(dxnvv*uu)+dxnuv*
     1          (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2                  *dxlvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlxnx+
     3                  dxluu*dxnvv))
        colyny=uvlyny*den+den3*(-0.5d0*(dyluv*(dynvv*uu)+dynuv*
     1          (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2                  *dylvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlyny+
     3                  dyluu*dynvv))
        colznz=uvlznz*den+den3*(-0.5d0*(dzluv*(dznvv*uu)+dznuv*
     1          (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2                  *dzlvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlznz+
     3                  dzluu*dznvv))
        colxny=uvlxny*den+den3*(-0.5d0*(dxluv*(dynvv*uu)+dynuv*
     1          (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2                  *dxlvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlxny+
     3                  dxluu*dynvv))
        colxnz=uvlxnz*den+den3*(-0.5d0*(dxluv*(dznvv*uu)+dznuv*
     1          (dxluu*vv+dxlvv*uu))+0.75d0*den2*uv*(dxluu*vv+uu
     2                  *dxlvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlxnz+
     3                  dxluu*dznvv))
        colynx=uvlynx*den+den3*(-0.5d0*(dyluv*(dxnvv*uu)+dxnuv*
     1          (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2                  *dylvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlynx+
     3                  dyluu*dxnvv))
        colynz=uvlynz*den+den3*(-0.5d0*(dyluv*(dznvv*uu)+dznuv*
     1          (dyluu*vv+dylvv*uu))+0.75d0*den2*uv*(dyluu*vv+uu
     2                  *dylvv)*(uu*dznvv)-0.5d0*uv*(uu*vvlynz+
     3                  dyluu*dznvv))
        colznx=uvlznx*den+den3*(-0.5d0*(dzluv*(dxnvv*uu)+dxnuv*
     1          (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2                  *dzlvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvlznx+
     3                  dzluu*dxnvv))
        colzny=uvlzny*den+den3*(-0.5d0*(dzluv*(dynvv*uu)+dynuv*
     1          (dzluu*vv+dzlvv*uu))+0.75d0*den2*uv*(dzluu*vv+uu
     2                  *dzlvv)*(uu*dynvv)-0.5d0*uv*(uu*vvlzny+
     3                  dzluu*dynvv))


        comxnx=uvmxnx*den+den3*(-0.5d0*(dxmuv*(dxnvv*uu)+dxnuv*
     1          (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2                  *dxmvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmxnx+
     3                  dxmuu*dxnvv))
        comyny=uvmyny*den+den3*(-0.5d0*(dymuv*(dynvv*uu)+dynuv*
     1          (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2                  *dymvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmyny+
     3                  dymuu*dynvv))
        comznz=uvmznz*den+den3*(-0.5d0*(dzmuv*(dznvv*uu)+dznuv*
     1          (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2                  *dzmvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmznz+
     3                  dzmuu*dznvv))
        comxny=uvmxny*den+den3*(-0.5d0*(dxmuv*(dynvv*uu)+dynuv*
     1          (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2                  *dxmvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmxny+
     3                  dxmuu*dynvv))
        comxnz=uvmxnz*den+den3*(-0.5d0*(dxmuv*(dznvv*uu)+dznuv*
     1          (dxmuu*vv+dxmvv*uu))+0.75d0*den2*uv*(dxmuu*vv+uu
     2                  *dxmvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmxnz+
     3                  dxmuu*dznvv))
        comynx=uvmynx*den+den3*(-0.5d0*(dymuv*(dxnvv*uu)+dxnuv*
     1          (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2                  *dymvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmynx+
     3                  dymuu*dxnvv))
        comynz=uvmynz*den+den3*(-0.5d0*(dymuv*(dznvv*uu)+dznuv*
     1          (dymuu*vv+dymvv*uu))+0.75d0*den2*uv*(dymuu*vv+uu
     2                  *dymvv)*(uu*dznvv)-0.5d0*uv*(uu*vvmynz+
     3                  dymuu*dznvv))
        comznx=uvmznx*den+den3*(-0.5d0*(dzmuv*(dxnvv*uu)+dxnuv*
     1          (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2                  *dzmvv)*(uu*dxnvv)-0.5d0*uv*(uu*vvmznx+
     3                  dzmuu*dxnvv))
        comzny=uvmzny*den+den3*(-0.5d0*(dzmuv*(dynvv*uu)+dynuv*
     1          (dzmuu*vv+dzmvv*uu))+0.75d0*den2*uv*(dzmuu*vv+uu
     2                  *dzmvv)*(uu*dynvv)-0.5d0*uv*(uu*vvmzny+
     3                  dzmuu*dynvv))

c Done with off diaggonal elements :-)

C DIAGGONAL ELEMENTS OF THE SECOND DERIVATIVE MATRIX
                kk = 6*(k-1)+1
                ll = 6*(l-1)+1
                mm = 6*(m-1)+1
                nn = 6*(n-1)+1
c@
C@              write(*,*)' kk cokx cokxkx ',kk,cokx,cokxkx
                diag(kk) = diag(kk) + dcos2*cokx*cokx + d2cos*cokxkx
c@
C@              write(*,*)'#1 diag(1) =',diag(1)
                diag(kk+1) = diag(kk+1) + dcos2*cokx*coky + d2cos*cokxky
                diag(kk+2) = diag(kk+2) + dcos2*cokx*cokz + d2cos*cokxkz
                diag(kk+3) = diag(kk+3) + dcos2*coky*coky + d2cos*cokyky
                diag(kk+4) = diag(kk+4) + dcos2*coky*cokz + d2cos*cokykz
                diag(kk+5) = diag(kk+5) + dcos2*cokz*cokz + d2cos*cokzkz

                diag(ll) = diag(ll) + dcos2*colx*colx + d2cos*colxlx
                diag(ll+1) = diag(ll+1) + dcos2*colx*coly + d2cos*colxly
                diag(ll+2) = diag(ll+2) + dcos2*colx*colz + d2cos*colxlz
                diag(ll+3) = diag(ll+3) + dcos2*coly*coly + d2cos*colyly
                diag(ll+4) = diag(ll+4) + dcos2*coly*colz + d2cos*colylz
                diag(ll+5) = diag(ll+5) + dcos2*colz*colz + d2cos*colzlz

                diag(mm) = diag(mm) + dcos2*comx*comx + d2cos*comxmx
                diag(mm+1) = diag(mm+1) + dcos2*comx*comy + d2cos*comxmy
                diag(mm+2) = diag(mm+2) + dcos2*comx*comz + d2cos*comxmz
                diag(mm+3) = diag(mm+3) + dcos2*comy*comy + d2cos*comymy
                diag(mm+4) = diag(mm+4) + dcos2*comy*comz + d2cos*comymz
                diag(mm+5) = diag(mm+5) + dcos2*comz*comz + d2cos*comzmz

                diag(nn) = diag(nn) + dcos2*conx*conx + d2cos*conxnx
                diag(nn+1) = diag(nn+1) + dcos2*conx*cony + d2cos*conxny
                diag(nn+2) = diag(nn+2) + dcos2*conx*conz + d2cos*conxnz
                diag(nn+3) = diag(nn+3) + dcos2*cony*cony + d2cos*conyny
                diag(nn+4) = diag(nn+4) + dcos2*cony*conz + d2cos*conynz
                diag(nn+5) = diag(nn+5) + dcos2*conz*conz + d2cos*conznz

C OFF DIAGONAL SECOND DERIVATIVE MATRIX ELEMENTS

                iiphi = 54*(iphi-1)+1

c kl
                d2phi(iiphi) = dcos2*cokx*colx + d2cos*cokxlx
                d2phi(iiphi+1) = dcos2*coky*colx + d2cos*cokylx
                d2phi(iiphi+2) = dcos2*cokz*colx + d2cos*cokzlx
                d2phi(iiphi+3) = dcos2*cokx*coly + d2cos*cokxly
                d2phi(iiphi+4) = dcos2*coky*coly + d2cos*cokyly
                d2phi(iiphi+5) = dcos2*cokz*coly + d2cos*cokzly
                d2phi(iiphi+6) = dcos2*cokx*colz + d2cos*cokxlz
                d2phi(iiphi+7) = dcos2*coky*colz + d2cos*cokylz
                d2phi(iiphi+8) = dcos2*cokz*colz + d2cos*cokzlz
c km
                d2phi(iiphi+9) = dcos2*cokx*comx + d2cos*cokxmx
                d2phi(iiphi+10) = dcos2*coky*comx + d2cos*cokymx
                d2phi(iiphi+11) = dcos2*cokz*comx + d2cos*cokzmx
                d2phi(iiphi+12) = dcos2*cokx*comy + d2cos*cokxmy
                d2phi(iiphi+13) = dcos2*coky*comy + d2cos*cokymy
                d2phi(iiphi+14) = dcos2*cokz*comy + d2cos*cokzmy
                d2phi(iiphi+15) = dcos2*cokx*comz + d2cos*cokxmz
                d2phi(iiphi+16) = dcos2*coky*comz + d2cos*cokymz
                d2phi(iiphi+17) = dcos2*cokz*comz + d2cos*cokzmz

c kn
                d2phi(iiphi+18) = dcos2*cokx*conx + d2cos*cokxnx
                d2phi(iiphi+19) = dcos2*coky*conx + d2cos*cokynx
                d2phi(iiphi+20) = dcos2*cokz*conx + d2cos*cokznx
                d2phi(iiphi+21) = dcos2*cokx*cony + d2cos*cokxny
                d2phi(iiphi+22) = dcos2*coky*cony + d2cos*cokyny
                d2phi(iiphi+23) = dcos2*cokz*cony + d2cos*cokzny
                d2phi(iiphi+24) = dcos2*cokx*conz + d2cos*cokxnz
                d2phi(iiphi+25) = dcos2*coky*conz + d2cos*cokynz
                d2phi(iiphi+26) = dcos2*cokz*conz + d2cos*cokznz

c lm
                d2phi(iiphi+27) = dcos2*colx*comx + d2cos*colxmx
                d2phi(iiphi+28) = dcos2*coly*comx + d2cos*colymx
                d2phi(iiphi+29) = dcos2*colz*comx + d2cos*colzmx
                d2phi(iiphi+30) = dcos2*colx*comy + d2cos*colxmy
                d2phi(iiphi+31) = dcos2*coly*comy + d2cos*colymy
                d2phi(iiphi+32) = dcos2*colz*comy + d2cos*colzmy
                d2phi(iiphi+33) = dcos2*colx*comz + d2cos*colxmz
                d2phi(iiphi+34) = dcos2*coly*comz + d2cos*colymz
                d2phi(iiphi+35) = dcos2*colz*comz + d2cos*colzmz
c ln
                d2phi(iiphi+36) = dcos2*colx*conx + d2cos*colxnx
                d2phi(iiphi+37) = dcos2*coly*conx + d2cos*colynx
                d2phi(iiphi+38) = dcos2*colz*conx + d2cos*colznx
                d2phi(iiphi+39) = dcos2*colx*cony + d2cos*colxny
                d2phi(iiphi+40) = dcos2*coly*cony + d2cos*colyny
                d2phi(iiphi+41) = dcos2*colz*cony + d2cos*colzny
                d2phi(iiphi+42) = dcos2*colx*conz + d2cos*colxnz
                d2phi(iiphi+43) = dcos2*coly*conz + d2cos*colynz
                d2phi(iiphi+44) = dcos2*colz*conz + d2cos*colznz
c mn
                d2phi(iiphi+45) = dcos2*comx*conx + d2cos*comxnx
                d2phi(iiphi+46) = dcos2*comy*conx + d2cos*comynx
                d2phi(iiphi+47) = dcos2*comz*conx + d2cos*comznx
                d2phi(iiphi+48) = dcos2*comx*cony + d2cos*comxny
                d2phi(iiphi+49) = dcos2*comy*cony + d2cos*comyny
                d2phi(iiphi+50) = dcos2*comz*cony + d2cos*comzny
                d2phi(iiphi+51) = dcos2*comx*conz + d2cos*comxnz
                d2phi(iiphi+52) = dcos2*comy*conz + d2cos*comynz
                d2phi(iiphi+53) = dcos2*comz*conz + d2cos*comznz

100     continue

c All done (?!)

        return
        end
