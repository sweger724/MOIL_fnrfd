      subroutine eCG_Torsion(si,co,itor,E,dE_co,dE_si)
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c calculate  torsion energy and forces for i-th torsion.
c in dE_co return dE/d(cos) 
c in dE_si return dE/d(sin)

      double precision si,co,db,dE_co, dE_si,E, dG1, dG2
      integer itor

      double precision si2,si3,si4,si5,co2,co3,co4,co5

c  DATASTRUCTURE CGTor
c  column   1   2   3   4   5   6   7   8   9   10     
c  data    a1  a2  a3  a4  a5  b1  b2  b3  b4   b5     
c
c
c   E = a1 cos(x) + a2 cos(2x) + a3 cos(3x) + a4 cos(4x) + a5 cos(5x)
c     + b1 sin(x) + b2 sin(2x) + b3 sin(3x) + b4 cos(4x) + b5 cos(5x)
c
      co2 = co**2
      co3 = co2*co
      co4 = co3*co
      co5 = co4*co

      si2 = si**2
      si3 = si2*si
      si4 = si3*si
      si5 = si4*si

        E = CGTor(itor,1)*co + CGTor(itor,2)*(2.d0*co2 -1.d0) 
     &    + CGTor(itor,3) * co *(4.d0*co2 -3.d0)
     &    + CGTor(itor,6)*si + CGTor(itor,7)*2.d0*co*si
     &    + CGTor(itor,8)*si * (3.d0 - 4.d0*si2)

        E = E + CGTor(itor,4)*(8.d0*co4 - 8.d0*co2 + 1.d0)
     &        + CGTor(itor,5)*(16.d0*co5 - 20.d0*co3 + 5.d0*co)
     &        + CGTor(itor,9)*(4.d0*co*si*(1.d0-2.d0*si2))
     &        + CGTor(itor,10)*(16.d0*si5 - 20.d0*si3 + 5.d0*si)

       dE_co = CGTor(itor,1) + CGTor(itor,2)*4.d0*co
     &    + CGTor(itor,3) *(12.d0*co2 -3.d0) 
     &    + CGTor(itor,7)*2.d0*si 

       dE_co = dE_co + CGTor(itor,4)*(32.d0*co3-16.d0*co)
     &       + CGTor(itor,5)*(80.d0*co4-60.d0*co2+5.d0)
     &       + CGTor(itor,9)*(4.d0*si-8.d0*si3)

       dE_si = CGTor(itor,6) + CGTor(itor,7)*2.d0*co
     &       + CGTor(itor,8)*(3.d0 - 12.d0*si2)

       dE_si = dE_si + CGTor(itor,9)*co*(4.d0-24.d0*si2)
     &       + CGTor(itor,10)*(80.d0*si4-60.d0*si2+5.d0)
    
      if ( E .gt. E_CG_max ) then
        write(6,*)"Torsion: co,si,Ene,dE_co,dE_si,c1,s1,type:",
     &               co,si, E, dE_co, dE_si, 
     &               CGTor(itor,1),CGTor(itor,6),itor
C        write(6,*)"Cs: ",CGTor(itor,1),CGTor(itor,2),CGTor(itor,3),
C     &                   CGTor(itor,4),CGTor(itor,5)
C        write(6,*)"Ss: ",CGTor(itor,6),CGTor(itor,7),CGTor(itor,8),
C     &                   CGTor(itor,9),CGTor(itor,10)
      endif

      return
      end
