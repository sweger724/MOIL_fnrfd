        subroutine eimphi()
        implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
c mauro add for muta
      include 'COMMON/MUTA.BLOCK'

c       calculate improper dihedral energies and forces according
c       to T.Schlick, J.Comput.Chem, vol 10, 7 (1989) 
c       with local variant (no delta funcion involved)

c   If |sin(phi)|.gt.1e-03 use the angle directly
c       V = kimp ( phi -impeq ) ^2
c       dV/dcos(phi) = -2*k/sin(phi)

c   If |sin(phi)|.lt.1e-03 and (phi-cnseq(iimp)).lt.0.1
c   and cnseq(phi)=0.
c       V = 2*kimp (1 - cos(phi))
c       dV/dcos(phi) = -2*kimp

c   If |sin(phi)|.lt.1e-03 and (phi-cnseq(iimp)).lt.0.1
c   and cnseq(phi)=pi
c       V = 2*kimp (1 + cos(phi))
c       dV/dcos(phi) = 2*kimp

c   Now, if |sin(phi)|.lt.1e-03 but far from minimum,
c   set sin(phi) to 1e-03 arbitrarily

c       Note that VECTOR.BLOCK contains temporary vectors used for
c       vectorization. Limit of 15000 proper dihedrals.

        double precision e
        double precision co ,phi,den,co1
        double precision uu,vv,uv
        double precision ax,bx,cx
        double precision ay,by,cy
        double precision az,bz,cz
        double precision ux,uy,uz
        double precision vx,vy,vz 
        double precision dx1,dy1,dz1 
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
        double precision df,delta,pi,yy
        double precision aa,bb,cc,ab,bc,ac
        double precision d1,d2,d3,e1,e2,e3
        double precision uu2, vv2
        integer iimp
        integer i,j,k,l


c       e_imp = total improper torsion energy (ENERGY.BLOCK)
c       e = improper torsion energy (non-acumulated)
c       co = cos(phi)

c       define function for dot product calculation
c       dx,ex = dummy variables

        double precision dot
        dot(d1,d2,d3,e1,e2,e3) = d1*e1 + d2*e2 + d3*e3


        pi=4.d0*datan(1.d0)

c       initialize e_imp
        e_imp = 0.d0

        call init_var()

c       initialize loop over improper torsions in ichunk chunks
        do 100 iimp = 1,nimp

        
c       call the torsion atoms
c assume the following arrangements
c      K   J          J   K
c       \ /            \ /
c        I              I
c        |              |
c        L              L
c they will give exactly the same improper angle.
c       ileana -- modification in atom order i,j,k,l
           if(.not.arith) then
                        i = iimp1(iimp)
                        j = iimp2(iimp)
                        k = iimp3(iimp)
                        l = iimp4(iimp)

                        !k = iimp1(iimp)
                        !j = iimp2(iimp)
                        !l = iimp3(iimp)
                        !i = iimp4(iimp)


            else

                           k=iimp1(iimp)
                           j=iimp2(iimp)
                           l=iimp3(iimp)
                           i=iimp4(iimp)
             endif

c       calculate distances
c 
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

                        if (co.gt.1.d0) co = 1.d0 
                        if (co.lt.-1.d0) co = -1.d0 


c       ileana -- if not arith use the normal MOIL equation for improper torsion energy

                        if(.not.arith) then


c There are 3 options here:
c (a) impeq=0   V = 2*k*(1-cos(phi))
                if (dabs(impeq(iimp)).lt.1.d-6) then
                         df = -2.d0*kimp(iimp)
                         e  = -df*(1.d0-co)
                        
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

                         
                         !write(6,*)"eee",e,phi * 180/pi
                         !write(6,*)"fff",i,j,k,l
c (b) impeq=pi  V = 2*k*(1+cos(phi))
                else if (dabs(dabs(impeq(iimp))-pi).lt.1.d-6) then
                        df = 2.d0*kimp(iimp)
                        e  = df*(1.d0+co)

c (c) impeq=something else V = k(phi-impeq)^2
                else

c       now calculate sin(phi) because cos(phi) is symmetric, so
c       we can decide between +-phi.

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

c       calculate deviation form equilibrium
                        delta = phi-impeq(iimp)
c       check for |delta| not to be bigger that pi
                        if (delta.gt.pi) delta = delta - 2.d0*pi
                        if (delta.lt.-pi) delta = delta + 2.d0*pi



                        df = kimp(iimp)* delta
                        e = df*delta
                        yy = dsin(phi)

c       the 1e-03 value could be lowered for increased precision.
c       This insures ~1e-05% error for sin(phi)=0.001
c       Increased accuracy by setting boundary to 1.d-6 instead
c       of 1.d-3
                        if ( abs(yy).gt.1.d-06) then

                                df = -2.d0*df/yy

                        else
                  if((phi.gt.0.d0 .and. phi.lt.(pi/2.d0)) .or.
     *            (phi.lt.0.d0 .and. phi.gt.(-pi/2.d0)) ) then
                        df = df*1d06
                  else
                        df = -df*1d06
                  end if        
                        end if
                
        if (debug ) then
	write(*,*)' i j k l',i,j,k,l
        write(*,*)' phi impeq ',phi*180.d0/pi,impeq(iimp)*180.d0/pi
                write(*,*)' delta ',delta*180.d0/pi
                write(*,*)' df = ',df
                write(*,*)' e yy = ',e,yy
        end if

                end if


                else
c       -- that is, if arith
c       ileana-- one single option, e = k(1+cos(n*phi -gamma))
c       ileana-- since always n=2 and gamma = pi in the new Amber ALL.PROP, e = k - kcos(2phi)


                        e = kimp(iimp)-kimp(iimp)*(2.d0*co*co-1.d0)
                         df = -4.d0*kimp(iimp)*co
                   

c       ileana -- just printing individual terms

c               write(stdo,*)'iimp is ',iimp
c               write(stdo,*)'cosine phi is ',co
c               write(stdo,*)'k is ',kimp(iimp)

c       now calculate sin(phi) because cos(phi) is symmetric, so
c       we can decide between +-phi.

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


c               write(stdo,200)i,j,k,l,e,phi*180.d0/pi 
 200            format('impropers',4(i5,2x),f8.3, " phi",f8.3)
                         

                         endif
c       end "if arith"

c mauro - <U1-U2>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1
     &                 .or. mutaid(l).eq.1 .or. mutaid(k).eq.1) then
                           tmp_e_lambda=tmp_e_lambda+e/lambda
                        endif
                        if (mutaid(i).eq.2 .or. mutaid(j).eq.2
     &                 .or. mutaid(l).eq.2 .or. mutaid(k).eq.2) then
                           tmp_e_lambda=tmp_e_lambda-e/(1.0d0-lambda)
                        endif
c end mauro - <U1-U2>l


        e_imp = e_imp + e



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
                        dpot(1,i) = dpot(1,i) + df * dd1x
                        dpot(2,i) = dpot(2,i) + df * dd1y
                        dpot(3,i) = dpot(3,i) + df * dd1z

                        dpot(1,j) = dpot(1,j) + df * dd2x
                        dpot(2,j) = dpot(2,j) + df * dd2y
                        dpot(3,j) = dpot(3,j) + df * dd2z

                        dpot(1,k) = dpot(1,k) + df * dd3x
                        dpot(2,k) = dpot(2,k) + df * dd3y
                        dpot(3,k) = dpot(3,k) + df * dd3z

                        dpot(1,l) = dpot(1,l) + df * dd4x
                        dpot(2,l) = dpot(2,l) + df * dd4y
                        dpot(3,l) = dpot(3,l) + df * dd4z


100     continue

        return
        end

