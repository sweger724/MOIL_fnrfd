      subroutine alabarrier()

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'

C     Create a gaussian barrier in psi at psi = 180 degrees for
c     alanine dipeptide. 
C     IMPORTANT: This code assumes that arith is off (false).
C     NOTE: This file is modeled after ctors.f. 9-May-2005, -TF.

      character*20 name
      integer namel, level
      double precision psi,k,sig,del,dels,f
      double precision dpsi(3,4)
      integer i, j, psii(4)


      name = 'alabarrier'
      namel = 10

      if (arith) then
         level = 1
         call alert(name,namel,'alabarrier called with arith',30,level)
      end if
      
      psii(1) = 4
      psii(2) = 6
      psii(3) = 8
      psii(4) = 10

      call torsion(psii(1),psii(2),psii(3),psii(4),psi,dpsi,.true.)
c      write(*,*) 'psi: ', psi

c     barrier height
      k = 20.d0

c     standard deviation (radians)
      sig = 0.05
      sig = 2.d0 * sig * sig

c     write(*,*) 'alabar k, sig:', k, sig

c     this places the barrier at +/- pi
c      if ( psi .ge. 0 ) then
c         del = psi - 3.14159
c      else
c         del = psi + 3.14159
c      end if

c     this places the barrier at psi = -100 degrees
c      del = psi + 1.7452 
      del = psi - abarpsi

      dels = del / sig
c     write(*,*) 'alabar dels: ', dels
      
c     energy
      e_alabar = k * exp ( -del * dels )

c      write(*,*) 'alabar energy: ', e

c     part of the force
      f = -2.d0 * e_alabar * dels

c     force
      do 1 i = 1,4
         do 2 j = 1,3
            dpot(j,psii(i)) = dpot(j,psii(i)) + f * dpsi(j,i)
 2       continue
 1    continue


      end
