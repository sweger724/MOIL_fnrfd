      subroutine ecent()

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'

c calculate center of mass constraint for selected particles.


      double precision xfor,yfor,zfor,xcm,ycm,zcm,xdif,ydif,zdif
      double precision invcent
      integer i,j

      vircent=0.d0

      xcm=0.d0
      ycm=0.d0
      zcm=0.d0

      invcent=dble(1.d0/icenter)
      do 100 i=1,icenter
	 j = center(i)
	 xcm=coor(1,j)+xcm
	 ycm=coor(2,j)+ycm
	 zcm=coor(3,j)+zcm
100   continue
      xdif=xcm*invcent-xeq
      ydif=ycm*invcent-yeq
      zdif=zcm*invcent-zeq

c      write (*,*) 'kc,xcm,ycm,zcm,xeq,yeq,zeq'
c      write (*,*) kcenter,xcm,ycm,zcm,xeq,yeq,zeq

      e_cent=kcenter*(xdif*xdif+ydif*ydif+zdif*zdif)
     
      invcent = 2.d0 * kcenter * invcent

      xfor=xdif*invcent
      yfor=ydif*invcent
      zfor=zdif*invcent

c all particles in pick get same force increment
      do 300 i=1,icenter
		j = center(i)
		dpot(1,j) = dpot(1,j) + xfor
		dpot(2,j) = dpot(2,j) + yfor
		dpot(3,j) = dpot(3,j) + zfor
c
c accumulate virial
c here sign is as vircent = r*force = - r*dpot,
c so reasonable convention. 
c DLP has virteth = - r*force !!! 
c switched vircent to DLP convention in dyna.f
c for consistency with other terms

		vircent = vircent - coor(1,j)*xfor
     &                    - coor(2,j)*yfor
     &                    - coor(3,j)*zfor

300   continue	
      return
      end
