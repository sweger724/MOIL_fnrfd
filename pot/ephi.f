        subroutine etors()
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
c mauro add for muta
      include 'COMMON/MUTA.BLOCK'

c       calculate dihedral energies and forces according
c       to T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 
c       with local variant (no delta funcion involved)

c       V = sum [ k(n)(1 + cos(n*phi + delta)]   (1=< n <= 3)
c       delta is either 0 or pi, so the formula become
c       less complicated upon expansion.


c       Note that VECTOR.BLOCK contains temporary vectors used for
c       vectorization. Limit of 15000 proper dihedrals.

        double precision e
        double precision co ,den,co1
        double precision uu,vv,uv
        double precision ax,bx,cx
        double precision ay,by,cy
        double precision az,bz,cz
        double precision a0x,b0x,c0x
        double precision a0y,b0y,c0y
        double precision a0z,b0z,c0z
        double precision a1x,b1x
        double precision a1y,b1y
        double precision a1z,b1z
        double precision a2x,b2x
        double precision a2y,b2y
        double precision a2z,b2z
        double precision dd1x,dd2x,dd3x,dd4x
        double precision dd1y,dd2y,dd3y,dd4y
        double precision dd1z,dd2z,dd3z,dd4z
        double precision df
        double precision aa,bb,cc,ab,bc,ac
        double precision d1,d2,d3,e1,e2,e3
        double precision ktot
        double precision uu2,vv2, phi
        integer iphi
        integer i,j,k,l

c       e_tors = total proper torsion energy (ENERGY.BLOCK)
c       e = proper torsion energy (non-acumulated)
c       co = cos(phi)
c       ktot = k(1) + k(2) + k(3)


c       define function for dot product calculation
c       dx,ex = dummy variables

        double precision dot
        dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3


c       initialize etors
        e_tors = 0.d0

c       initialize loop over torsions in ichunk chunks
        do 100 iphi = 1,ntors

        

c       call the torsion atoms
                        i = itor1(iphi)
                        j = itor2(iphi)
                        k = itor3(iphi)
                        l = itor4(iphi)
c                       write(*,*)'i j k l  is',i,j,k,l
c                       
c                       write(*,*)'ktors is', ktors2(iphi)
c                       write(*,*)'phase1 ',phase1(iphi)

c       calculate distances
                        ax = coor(1,j) - coor(1,i)
                        ay = coor(2,j) - coor(2,i)
                        az = coor(3,j) - coor(3,i)

                        bx = coor(1,k) - coor(1,j)
                        by = coor(2,k) - coor(2,j)
                        bz = coor(3,k) - coor(3,j)

                        cx = coor(1,l) - coor(1,k)
                        cy = coor(2,l) - coor(2,k)
                        cz = coor(3,l) - coor(3,k)

c       calculate dot products
                        ab = dot(ax,ay,az,bx,by,bz)
                        bc = dot(bx,by,bz,cx,cy,cz)
                        ac = dot(ax,ay,az,cx,cy,cz)
                        aa = dot(ax,ay,az,ax,ay,az)
                        bb = dot(bx,by,bz,bx,by,bz)
                        cc = dot(cx,cy,cz,cx,cy,cz)


c       calculate cos(phi)      
                        uu = (aa * bb) - (ab*ab)
                        vv = (bb * cc ) - (bc * bc)
                        uv = (ab * bc) - (ac * bb)
                        den= 1.d0/dsqrt(uu*vv) 
                        co = uv * den

c       calculate energies

c       ileana

                        if(arith) then
                        
                           

              if(period(iphi).eq.1) then
                 e = ktors2(iphi) +phase2(iphi)* ktors2(iphi)*co
              endif

              if(period(iphi).eq.2) then

               e = ktors2(iphi) +phase2(iphi)*ktors2(iphi)
     &         *(2.0d0*co*co -1.0d0)
               endif

               if(period(iphi).eq.3)then

               e = ktors2(iphi) +phase2(iphi)*
     &         ktors2(iphi)*co*(4.0d0*co*co-3.0d0)
               endif

               if(period(iphi).eq.4)then

               e = ktors2(iphi)+phase2(iphi)*ktors2(iphi)*
     &         (2.0d0*(2.0d0*co*co-1.0d0)*
     1         (2.0d0*co*co-1.0d0)-1.0d0)
               endif   



               else
               
                 e1 = ktors1(iphi) * co
                 e2 = ktors2(iphi) * ( 2.d0*co*co - 1.d0 )
                 e3 = ktors3(iphi) * co * ( 4.d0*co*co -3.d0 )
               
                 ktot = ktors1(iphi)+ktors2(iphi)+ktors3(iphi)
                 e = ktot + phase1(iphi)*e1 + phase2(iphi)*e2 
     1                    + phase3(iphi)*e3 
                        
               endif

c       ileana -- just print individual terms

C                       write(6,*)'itors is ',iphi
C                       write(*,*)'i j k l  is',i,j,k,l
C                      write(*,*)'N: ',ptnm(i),ptnm(j),ptnm(k),ptnm(l)
C                      write(*,*)'res index',poimon(i),poimon(j)
C     & ,poimon(k),poimon(l)
C                     write(*,*)'res name ',moname(poimon(i)),' ',
C     & moname(poimon(j)),' ',moname(poimon(k)),' ',moname(poimon(l))
C                       write(*,*)'ktors is ',ktors1(iphi),ktors2(iphi)
C     & ,ktors3(iphi)
C                       write(*,*)'cosine is',co
C                       write(6,'(a,f8.4)')'etors is ',e,e_tors
C                       write(*,*)'Phi is',phi*pi180 
C                       write(*,*)'~~~~~~~~~~~~~~~~~~~~~'
c mauro - <U1-U2>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1
     &                 .or. mutaid(k).eq.1 .or. mutaid(l).eq.1) then
                           tmp_e_lambda=tmp_e_lambda+e/lambda
                        endif
                        if (mutaid(i).eq.2 .or. mutaid(j).eq.2
     &                 .or. mutaid(k).eq.2 .or. mutaid(l).eq.2) then
                           tmp_e_lambda=tmp_e_lambda-e/(1.d0-lambda)
                        endif
c end mauro - <U1-U2>l

                        e_tors = e_tors + e
                       phi = dacos(co)


c       calculate derivatives by chain rule
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

c       calculate forces
c       ileana-- slight change in df due to k1-k2-k3 

                        if(.not.arith)then
                          df = phase1(iphi)*ktors1(iphi)
     1                   + phase2(iphi)*ktors2(iphi)*4.d0*co
     2                   + phase3(iphi)*ktors3(iphi)*(12.d0*co*co-3.d0)
                        
                        else

                        if(period(iphi).eq.1) then
                        df = ktors2(iphi) * phase2(iphi)
                        endif
                        if(period(iphi).eq.2) then
                        df = 4.0d0*ktors2(iphi) * phase2(iphi) *co
                        endif 
                        if(period(iphi).eq.3) then
                       df=ktors2(iphi)*phase2(iphi)*(12.0d0*co*co-3.0d0)
                        endif 
                        if(period(iphi).eq.4) then
                        df=16.0d0*ktors2(iphi)*phase2(iphi)*co*
     &                    (2.0d0*co*co-1.0d0)
                        endif 

                        
                        endif


                        dpot(1,i) = dpot(1,i) + df*dd1x
                        dpot(2,i) = dpot(2,i) + df*dd1y
                        dpot(3,i) = dpot(3,i) + df*dd1z

                        dpot(1,j) = dpot(1,j) + df*dd2x
                        dpot(2,j) = dpot(2,j) + df*dd2y
                        dpot(3,j) = dpot(3,j) + df*dd2z

                        dpot(1,k) = dpot(1,k) + df*dd3x
                        dpot(2,k) = dpot(2,k) + df*dd3y
                        dpot(3,k) = dpot(3,k) + df*dd3z

                        dpot(1,l) = dpot(1,l) + df*dd4x
                        dpot(2,l) = dpot(2,l) + df*dd4y
                        dpot(3,l) = dpot(3,l) + df*dd4z


100     continue

        return
        end













