      subroutine ctors()

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'

C     Calculate the value (e) and gradient (de) of a penalty function
C     that biases the phi and psi torsions of alanine dipeptide to
C     satisfy the linear constraint: phi_k*phi + psi_k*psi = phipsi_ceq.
C     The constant str determines the strength of the penalty: e = k(a
C     *phi + b*psi - c)^2.
C     Alanine dipeptide coordinates are assumed to be the first 12
C     coordinates of the global variable coor.
C     
C     IMPORTANT: This code assumes that arith is off (false).
C     NOTE: This file is modeled after pot/eimphi.f. 9-May-2005, -TF.

      character*8 name
      integer namel, level
      double precision z, f, junk
      double precision phi, psi
      double precision dphi(3,4), dpsi(3,4)
      integer i, j, phii(4), psii(4)


      name = 'ctors'
      namel = 5

      if (arith) then
         level = 1
         call alert(name,namel,'ctors called with arith',30,level)
      end if
      
      if ( abs(keq) .gt. 3.142 ) then
         level = 1
         call alert(name,namel,'kequ is not in radians',30,level)
      end if


      phii(1) = 2
      phii(2) = 4
      phii(3) = 6
      phii(4) = 8

      psii(1) = 4
      psii(2) = 6
      psii(3) = 8
      psii(4) = 10

      call torsion ( phii(1), phii(2), phii(3), phii(4), phi, dphi,
     $     .true. )
      call torsion ( psii(1), psii(2), psii(3), psii(4), psi, dpsi,
     $     .true. )

      z = (kphi * phi) + (kpsi * psi) - keq

c     energy
      e_ctors = ctstr * z * z


c     part of the force (factor of 2 supplied by torsion, below)
      f = 2.0 * ctstr * z

c     force
      do 1 i = 1,4
         do 2 j = 1,3
            dpot(j,phii(i)) = dpot(j,phii(i)) + f*kphi*dphi(j,i)
            dpot(j,psii(i)) = dpot(j,psii(i)) + f*kpsi*dpsi(j,i)
 2       continue
 1    continue


      end
