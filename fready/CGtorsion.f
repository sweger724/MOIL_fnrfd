        function CGtorsion(i,j,k,l,coor2)
        implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/FREADY.BLOCK'
      include 'COMMON/CONVERT.BLOCK'

c       calculate dihedral energies and forces

        double precision e, signn, at
        double precision co, den, co1, tmp
        double precision v_ab_x, v_ab_y, v_ab_z
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
        double precision df_co, df_si
        double precision aa,bb,cc,ab,bc,ac
        double precision d1,d2,d3,e1,e2,e3
        double precision uu2,vv2
        integer iphi, mm
        integer i,j,k,l
        double precision CGtorsion, coor2(3,maxpt)
        integer terminal(maxmono)


c       e_tors = total proper torsion energy (ENERGY.BLOCK)
c       e = proper torsion energy (non-acumulated)
c       co = cos(phi)
c       si = sin(phi)

c       define function for dot product calculation
c       dx,ex = dummy variables

        double precision dot
        dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3


c       calculate distances
                        ax = coor2(1,j) - coor2(1,i)
                        ay = coor2(2,j) - coor2(2,i)
                        az = coor2(3,j) - coor2(3,i)

                        bx = coor2(1,k) - coor2(1,j)
                        by = coor2(2,k) - coor2(2,j)
                        bz = coor2(3,k) - coor2(3,j)

                        cx = coor2(1,l) - coor2(1,k)
                        cy = coor2(2,l) - coor2(2,k)
                        cz = coor2(3,l) - coor2(3,k)


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

C check colinear angles
        if (uu .lt. 1.d-6) uu = 1.d-6
        if (vv .lt. 1.d-6) vv = 1.d-6

                        den= 1.d0/dsqrt(uu*vv) 
                        co = uv * den

                        
C       v_ab_? stands for particular coordinate of the vector product (a x b)
                        
                        v_ab_x =  ay*bz - az*by
                        v_ab_y = -ax*bz + az*bx
                        v_ab_z =  ax*by - ay*bx
                        
C  the equation of the plane defined by vectors a,b is
C  dot(v_ab)(x-j) = 0 ... if this for x:=l is > 0 we keep the sin(phi) 
C  as calculated before, otherwise we assign sin(phi) := -sin(phi)

        signn = dot(v_ab_x,v_ab_y,v_ab_z,coor2(1,l)-coor2(1,j)
     1           ,coor2(2,l)-coor2(2,j),coor2(3,l)-coor2(3,j))
        if (signn .lt. 0.d0)  then
           signn = -1.d0
        else 
           signn = 1.d0
        endif

        CGtorsion=dacos(co) 
        if (signn .lt. 0.d0)     CGtorsion = -CGtorsion
        !write(6,*)"CCC",CGtorsion,i,j,k,l

        return 
        end
