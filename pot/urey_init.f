      subroutine urey_init(ii)
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
c muta
      include 'COMMON/MUTA.BLOCK'

      double precision rx,ry,rz,r2,s,r,db,df,e
      double precision th,rji,rjk,ks,pi
      integer mm,i,j,k,ii,namel
      character*6 name

      namel=6
      name='gmanbo'

      i=iangl1(ii)
      j=iangl2(ii)
      k=iangl3(ii)

      pi = 4.d0*datan(1.d0)

c  add a spring between i and k
c
c      i           i
c     /           / \
c    j---k   =>  j---k
c     
c

      do 100 mm=1,nb
        write (777,*) ib1(mm),ib2(mm)

        if ((ib1(mm).eq.i.and.ib2(mm).eq.j).or.
     1      (ib1(mm).eq.j.and.ib2(mm).eq.i)) then
                rji = req(mm)
         else if ((ib1(mm).eq.k.and.ib2(mm).eq.j).or.
     1            (ib1(mm).eq.j.and.ib2(mm).eq.k)) then
                rjk = req(mm)
        endif
100     continue

       th = angleq(ii)
       ks = kangl(ii)

       angleq_urey(ii) = (rji**(2.d0)+rjk**(2.d0)-
     1              2.d0*rji*rjk*dcos(th))**(0.5d0)
       if (th.eq.0.d0 .or. th.eq.0.5d0*pi .or. th.eq.pi) then
        write (*,*) i,j,k,th
        call alert(name,namel,'non valid angle for Brady',25,1)
       endif
       kangl_urey(ii) = 
     1 ks*(1.d0-(2.d0*rji*rjk)/(rji**(2.d0)+rjk**(2.d0))*dcos(th))/
     2 (rji**(2.d0)*rjk**(2.d0)*dsin(th)**(2.d0)/
     3 (rji**(2.d0)+rjk**(2.d0)))

       kangl(ii) = 0.d0
       !write (777,'(3i5)') i,j,k
       !write (777,'(4f10.5)') rji,rjk,th,ks
       !write (777,'(2f10.5)') angleq(ii),kangl(ii)
       !write (777,*) '******************'

      return
      end
