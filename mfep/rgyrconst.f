      subroutine rgyrconst()

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/RGYRCST.BLOCK'

C     Energy penalty enforcing a given radius of gyration (squared).
c     NOTE: This code assumes that the center of mass is zero
C     NOTE: This code assumes that arith is off (false).
C     NOTE: This file is modeled after ctors.f. 22-Oct-2005, -TF.

      character*20 name
      integer namel, level
      double precision rg2, del, f, m
      integer i, j, k

c     rgk is the penalty strength
c     rg20 is the desired radius of gyration SQUARED


      name = 'rgyrconst'
      namel = 10

      if (arith) then
         level = 1
         call alert(name,namel,'rgyrconst called with arith',30,level)
      end if

      
c     radius of gyration squared
      rg2 = 0.d0
      m = 0.d0
      do 100 i = 1,rgcnpick
         j = rgcnpointr(i)
         m = m + ptms(j)
         do 101 k = 1,3
            rg2 = rg2 + ptms(j) * coor(1,j) * coor(1,j)
 101     continue
 100  continue
      m = 1.d0 / m
      rg2 = rg2 * m

      del = rg2 - rg20
            

c     energy
      e_rgyr = rgk * del * del


c     force
      f = 4 * rgk * del * m
      do 1 i = 1,rgcnpick
         j = rgcnpointr(i)
         do 2 k = 1,3
            dpot(k,j) = dpot(k,j) + f * ptms(j) * coor(k,j)
 2       continue
 1    continue


      end
