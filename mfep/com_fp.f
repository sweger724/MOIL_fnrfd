      subroutine com_fp (coor, ptms, cm, nselec, pointr)

c     IN
c     coor: [d n]: coordinates
c     ptms: [1 n]: particle masses
c     nselec, pointr: particle selection

c     OUT
c     cm: [d 1]: center of mass

      implicit none
      integer i, j, k, nselec, pointr(*)
      double precision cm(3), ptms(*)
      double precision coor(3,*)

c     total mass
      double precision m


      m = 0.d0
      do 300 j = 1, nselec
         i = pointr(j)
         m = m + ptms(i)
 300  continue
      m = 1.d0 / m


      do 100 k = 1,3
         cm(k) = 0.d0
         do 200 j = 1,nselec
            i = pointr(j)
            cm(k) = cm(k) + (coor(k,i) * ptms(i))
c            cm(k) = cm(k) + (coor(k,i))
 200     continue
         cm(k) = cm(k) * m
c         cm(k) = cm(k) / nselec
 100  continue


      end
