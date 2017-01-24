      subroutine urey()
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
c muta
      include 'COMMON/MUTA.BLOCK'

      double precision rx,ry,rz,r2,s,r,db,df,e
      integer mm,i,j,k,ith

      do 100 ith=1,nangl

c        if ((pnm_a(ith).ne.2.and.pnm_a(ith).ne.3)) goto 100
c  add a spring between i and k
c
c      i           i
c     /           / \
c    j---k   =>  j---k
c     
c


                i=iangl1(ith)
                j=iangl2(ith)
                k=iangl3(ith)

                rx=coor(1,i)-coor(1,k)
                ry=coor(2,i)-coor(2,k)
                rz=coor(3,i)-coor(3,k)
                r2=rx*rx + ry*ry + rz*rz
                s=dsqrt(r2)
                r=2.d0/s
                db=s-angleq_urey(ith)
                df=kangl_urey(ith)*db
                e=df*db
c@@@
C               write(6,*)'bond',i,j,s,req(mm)
c@@@

c muta - store different part of interaction in different
c        elements of the array:
c        1 P-P;
c        2 P-M;
c        3 P-N;
c        4 M-M;
c        5 N-N;
c        6 M-N;
c        if no muta everething goes in 1.

                !pnm_eb(pnm_a(ith)) = pnm_eb(pnm_a(ith)) + e 

                e_theta = e_theta + e
                !write (777,*) ith,e
                df=df*r
                dpot(1,i) = dpot(1,i) + rx*df
                dpot(2,i) = dpot(2,i) + ry*df
                dpot(3,i) = dpot(3,i) + rz*df
                dpot(1,k) = dpot(1,k) - rx*df
                dpot(2,k) = dpot(2,k) - ry*df
                dpot(3,k) = dpot(3,k) - rz*df

c        if (i.eq.25.or.k.eq.25) write (250,*) 'brady',i,k,rx*df


c               by assumption no bonds between slow and fast

100     continue

      return
      end
