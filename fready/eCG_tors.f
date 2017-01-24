        subroutine eCG_tors()
        implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c       calculate dihedral energies and forces

        double precision e, signn, at
        double precision co, si, den, co1, tmp
        double precision v_ab_x, v_ab_y, v_ab_z
        double precision uu,vv,uv
        double precision a(3),b(3),c(3)
        double precision a0(3),b0(3),c0(3)
        double precision a1(3),b1(3),a2(3),b2(3)
        double precision dd1(3),dd2(3),dd3(3),dd4(3)
        double precision df_co, df_si
        double precision aa,bb,cc,ab,bc,ac
        double precision d1,d2,d3,e1,e2,e3
        double precision uu2,vv2
        integer iphi
        integer i,j,k,l,m

C sinus force variables
        double precision mx,my,mz,nx,ny,nz,mnx,mny,mnz
        double precision s_m,s_n,s_mn
        double precision dmx(12),dmy(12),dmz(12),dnx(12),dny(12),dnz(12)
        double precision dmnx(12),dmny(12),dmnz(12)
        double precision dsm(12),dsn(12),dsmn(12),dsi(12)

c       e_tors = total proper torsion energy (ENERGY.BLOCK)
c       e = proper torsion energy (non-acumulated)
c       co = cos(phi)
c       si = sin(phi)

c       define function for dot product calculation
c       dx,ex = dummy variables

        integer T,mm
        double precision dot
        dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3


c       initialize etors
        e_tors = 0.d0
        mm = 0

c       initialize loop over torsions
        do iphi = 1,ntors

C This is part of our cooperatity model
C   only for CA-CA-CA-CA torsions
        if (torType(iphi).le.400) then
          mm = mm + 1
        endif

c       name the torsion atoms
                        i = itor1(iphi)
                        j = itor2(iphi)
                        k = itor3(iphi)
                        l = itor4(iphi) 

        if (CGstr(j).eq.0.0 .or. CGstr(k).eq.0.0
     &  .or. ( .not.Fix_2nd_structure) ) then

c       calculate distances
                     do m = 1, 3
                        a(m) = coor(m,j) - coor(m,i)
                        b(m) = coor(m,k) - coor(m,j)
                        c(m) = coor(m,l) - coor(m,k)
                     end do

c calculate vector products m = a x b , n = b x c

                        mx = a(2)*b(3) - a(3)*b(2)
                        my = a(3)*b(1) - a(1)*b(3)
                        mz = a(1)*b(2) - a(2)*b(1)

                        nx = b(2)*c(3) - b(3)*c(2)
                        ny = b(3)*c(1) - b(1)*c(3)
                        nz = b(1)*c(2) - b(2)*c(1)

c calculate vector product mn = m x n

                        mnx = my*nz - mz*ny
                        mny = mz*nx - mx*nz
                        mnz = mx*ny - my*nx

c calculate sizes of vectors m,n,mn

                        s_m = dsqrt(mx**2 + my**2 + mz**2)
                        s_n = dsqrt(nx**2 + ny**2 + nz**2)
                        s_mn = dsqrt(mnx**2 + mny**2 + mnz**2)

C check colinear angles
        if (s_m .lt. 1.d-5) s_m = 1.d-5
        if (s_n .lt. 1.d-5) s_n = 1.d-5
        if (s_mn .lt.1.d-6) s_mn = 1.d-6

c calculate sinus

                        si = s_mn / ( s_m * s_n )
        if (s_m .eq. 1.d-5 .or. s_n .eq. 1.d-5) si = 0.d0

c       calculate dot products
                        ab = dot(a(1),a(2),a(3),b(1),b(2),b(3))
                        bc = dot(b(1),b(2),b(3),c(1),c(2),c(3))
                        ac = dot(a(1),a(2),a(3),c(1),c(2),c(3))
                        aa = dot(a(1),a(2),a(3),a(1),a(2),a(3))
                        bb = dot(b(1),b(2),b(3),b(1),b(2),b(3))
                        cc = dot(c(1),c(2),c(3),c(1),c(2),c(3))

c       calculate cos(phi)      
                        uu = (aa * bb) - (ab*ab)
                        vv = (bb * cc ) - (bc * bc)
                        uv = (ab * bc) - (ac * bb)

C check colinear angles
        if (uu .lt. 1.d-6) uu = 1.d-6
        if (vv .lt. 1.d-6) vv = 1.d-6

                        den= 1.d0/dsqrt(uu*vv) 
                        co = uv * den

                        
c       calculate sin(phi)
C       v_ab_? stands for particular coordinate of the vector product (a x b)
                        
                        v_ab_x =  a(2)*b(3) - a(3)*b(2)
                        v_ab_y = -a(1)*b(3) + a(3)*b(1)
                        v_ab_z =  a(1)*b(2) - a(2)*b(1)
                        
C  the equation of the plane defined by vectors a,b is
C  dot(v_ab)(x-j) = 0 ... if this for x:=l is > 0 we keep the sin(phi) 
C  as calculated before, otherwise we assign sin(phi) := -sin(phi)

        signn = dot(v_ab_x,v_ab_y,v_ab_z,coor(1,l)-coor(1,j)
     1              ,coor(2,l)-coor(2,j),coor(3,l)-coor(3,j))
        if (signn .lt. 0.d0)  then
           signn = -1.d0
        else 
           signn = 1.d0
        endif

        si = signn * si
        
        at=dacos(co)
        if (signn .lt. 0.d0)     at = -at
        
c       calculate energy
        T = torType(iphi) 
        
C          write(6,*)"Dihedral angle: ",at * pi180, T   
          call eCG_Torsion(si,co,T,E,df_co,df_si)
        
        e_tors = e_tors + E


C      Add cooperativity derivatives (see eCG_coop.f)
C        if (T.le.400) then
C          df_co = df_co + ddco(mm)
C          df_si = df_si + ddsi(mm)
C        endif


C                       write(*,*)'itors is   ',iphi, i,j,k,l
C                       write(*,*)'i is       ',coor(1,i),coor(2,i),coor(3,i)
C                       write(*,*)'j is       ',coor(1,j),coor(2,j),coor(3,j)
C                       write(*,*)'k is       ',coor(1,k),coor(2,k),coor(3,k)
C                       write(*,*)'l is       ',coor(1,l),coor(2,l),coor(3,l)
C                       write(*,*)'cosine is  ',co, uu, vv,bb,cc,bc
C                       write(*,*)'sine is    ',si
C                       write(*,*)'etors is   ',e
C       write(*,*)'A,B: ',CGtor(iphi,1),CGtor(iphi,2)
C                       write(*,*)'~~~~~~~~~~~~~~~~~~~~~'

c       calculate derivatives by chain rule
                        co1 = 0.5d0*co*den

                        do m = 1,3
                          a0(m) = -bc*b(m) + bb*c(m)
                          b0(m) =  ab*c(m) + bc*a(m) -2.d0*ac*b(m)
                          c0(m) =  ab*b(m) - bb*a(m)
                        end do
        
                        uu2 = 2.d0*uu
                        vv2 = 2.d0*vv

                        do m = 1,3
                          a1(m) = uu2*(-cc*b(m) + bc*c(m))
                          b1(m) = uu2*( bb*c(m) - bc*b(m))
                          a2(m) = -vv2*(bb*a(m) - ab*b(m))
                          b2(m) = vv2*(aa*b(m) - ab*a(m))
                        end do

                        do m = 1,3
                          dd1(m) = (a0(m) - a2(m)*co1)*den

                          dd2(m) = (-a0(m) - b0(m) - 
     &                             (a1(m) - a2(m) - b2(m))*co1)*den

                          dd3(m) = (b0(m) - c0(m) - 
     &                             (-a1(m) - b1(m) + b2(m))*co1)*den

                          dd4(m) = (c0(m) - b1(m)*co1)*den 
                        end do

c ******************************************************************
c Calculate derivatives of  sin(phi) with respect to 12 positions of 
c coresponding points defyning the angle phi.
C ******************************************************************
        dmx(1) = 0.d0
        dmx(2) = -b(3)
        dmx(3) =  b(2)

        dmx(4) = 0.d0
        dmx(5) =  b(3) + a(3)
        dmx(6) = -b(2) - a(2)

        dmx(7) = 0.d0
        dmx(8) = -a(3)
        dmx(9) =  a(2)

        dmx(10) = 0.d0
        dmx(11) = 0.d0
        dmx(12) = 0.d0



        dmy(1) =  b(3)
        dmy(2) = 0.d0
        dmy(3) = -b(1)

        dmy(4) = -a(3) - b(3)
        dmy(5) = 0.d0
        dmy(6) =  b(1) + a(1)

        dmy(7) =  a(3)
        dmy(8) = 0.d0
        dmy(9) = -a(1)

        dmy(10) = 0.d0
        dmy(11) = 0.d0
        dmy(12) = 0.d0

        
        dmz(1) = -b(2)
        dmz(2) =  b(1)
        dmz(3) = 0.d0

        dmz(4) =  b(2) + a(2)
        dmz(5) = -b(1) - a(1)
        dmz(6) = 0.d0

        dmz(7) = -a(2)
        dmz(8) =  a(1)
        dmz(9) = 0.d0

        dmz(10) = 0.d0
        dmz(11) = 0.d0
        dmz(12) = 0.d0


        dnx(1) = 0.d0
        dnx(2) = 0.d0
        dnx(3) = 0.d0

        dnx(4) = 0.d0
        dnx(5) = -c(3)
        dnx(6) =  c(2)

        dnx(7) = 0.d0
        dnx(8) =  c(3) + b(3)
        dnx(9) = -c(2) - b(2)

        dnx(10) = 0.d0
        dnx(11) = -b(3)
        dnx(12) =  b(2)


        dny(1) = 0.d0
        dny(2) = 0.d0
        dny(3) = 0.d0

        dny(4) =  c(3)
        dny(5) = 0.d0
        dny(6) = -c(1)

        dny(7) = -b(3) - c(3)
        dny(8) = 0.d0
        dny(9) =  c(1) + b(1)

        dny(10) =  b(3)
        dny(11) = 0.d0
        dny(12) = -b(1)


        dnz(1) = 0.d0
        dnz(2) = 0.d0
        dnz(3) = 0.d0

        dnz(4) = -c(2)
        dnz(5) =  c(1)
        dnz(6) = 0.d0

        dnz(7) =  c(2) + b(2)
        dnz(8) = -c(1) - b(1)
        dnz(9) = 0.d0

        dnz(10) = -b(2)
        dnz(11) =  b(1)
        dnz(12) = 0.d0


        do t =1,12
          dmnx(t) = dmy(t)*nz + my*dnz(t) - ( dmz(t)*ny + mz*dny(t) ) 
          dmny(t) = dmz(t)*nx + mz*dnx(t) - ( dmx(t)*nz + mx*dnz(t) )
          dmnz(t) = dmx(t)*ny + mx*dny(t) - ( dmy(t)*nx + my*dnx(t) )
          dsm(t)  = ( mx * dmx(t) + my * dmy(t) + mz * dmz(t) )/s_m
          dsn(t)  = ( nx * dnx(t) + ny * dny(t) + nz * dnz(t) )/s_n
          dsmn(t) = ( mnx*dmnx(t) + mny*dmny(t) + mnz*dmnz(t) )/s_mn
          dsi(t)  = signn * ( dsmn(t)/(s_m*s_n) 
     &            - s_mn*(dsm(t)*s_n + dsn(t)*s_m)/((s_m*s_n)**2) )
        enddo

c       calculate forces
                        
           do m = 1,3
             dpot(m,i) = dpot(m,i) + df_co*dd1(m) + df_si*dsi(m)
             dpot(m,j) = dpot(m,j) + df_co*dd2(m) + df_si*dsi(m+3)
             dpot(m,k) = dpot(m,k) + df_co*dd3(m) + df_si*dsi(m+6)
             dpot(m,l) = dpot(m,l) + df_co*dd4(m) + df_si*dsi(m+9)
           end do

          endif   

        end do

        return
        end
