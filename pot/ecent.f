      subroutine ecent()
        implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'

c calculate center of mass constraint for selected particles.


      double precision xfor,yfor,zfor,xcm,ycm,zcm,xdif,ydif,zdif
      double precision invcent
      integer i,j

      xcm=0.d0
      ycm=0.d0
      zcm=0.d0

      xdif=0.d0
      ydif=0.d0
      zdif=0.d0

      invcent=dble(1.d0/icenter)
      do 100 i=1,icenter
         j = center(i)
         if (xcenter) xcm=coor(1,j)+xcm
         if (ycenter) ycm=coor(2,j)+ycm
         if (zcenter) zcm=coor(3,j)+zcm
100   continue
      if (xcenter) xdif=xcm*invcent-xeq
      if (ycenter) ydif=ycm*invcent-yeq
      if (zcenter) zdif=zcm*invcent-zeq

c      write (*,*) 'xf yf zf',xcenter,ycenter,zcenter
c      write (*,*) 'xcm=',xcm,xdif,invcent
c      write (*,*) 'ycm=',ycm,ydif,invcent
c      write (*,*) 'zcm=',zcm,zdif,invcent

c      write (*,*) 'kc,xcm,ycm,zcm,xeq,yeq,zeq'
c      write (*,*) kcenter,xcm,ycm,zcm,xeq,yeq,zeq
      e_cent=kcenter*(xdif*xdif+ydif*ydif+zdif*zdif)
     
      invcent = 2.d0 * kcenter * invcent

c      write (*,*) 'kcenter and e_cent',kcenter,e_cent

      xfor=xdif*invcent
      yfor=ydif*invcent
      zfor=zdif*invcent

c      write (*,*) 'xfor=',xfor
c      write (*,*) 'yfor=',yfor
c      write (*,*) 'zfor=',zfor

c all particles in pick get same force increment
      do 300 i=1,icenter
                j = center(i)
                if (xcenter) dpot(1,j) = dpot(1,j) + xfor
                if (ycenter) dpot(2,j) = dpot(2,j) + yfor
                if (zcenter) dpot(3,j) = dpot(3,j) + zfor
300   continue  
      return
      end
