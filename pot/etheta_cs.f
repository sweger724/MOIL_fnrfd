      subroutine etheta()
c---This is etheta.cos written by  M. Kara-Ivanov
c-- it will help to get rid of singularities sintheta appr. 0
c---in the case if the angles are close to 180 deg !
c---Idea is to use E=k(cos - coseq)*(cos-coseq) instead of angles!

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'

c calculate angle energies and forces

c note that vector.block contains temporary vectors used for
c vectorization.

      double precision dxi,dyi,dzi,dxj,dyj,dzj,ri2,rj2
      double precision ri,rj,rir,rjr,dxir,dyir,dzir,dxjr,dyjr,dzjr
      double precision cst,da,df,e,tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2
      integer ith,i,j,k


c initialize e_theta
      e_theta=0.0

c initialize loop over angles in ichunk chunks
      do 100 ith=1,nangl


                i=iangl1(ith)
                j=iangl2(ith)
                k=iangl3(ith)

		dxi = coor(1,i) - coor(1,j)
		dyi = coor(2,i) - coor(2,j)
		dzi = coor(3,i) - coor(3,j)

		dxj = coor(1,k) - coor(1,j)
		dyj = coor(2,k) - coor(2,j)
		dzj = coor(3,k) - coor(3,j)

      		ri2=dxi*dxi+dyi*dyi+dzi*dzi
      		rj2=dxj*dxj+dyj*dyj+dzj*dzj
      		ri=dsqrt(ri2)
      		rj=dsqrt(rj2)

c		if (ri*rj.eq.0) then
c                   write(*,*)i,j,ri
c                   goto 200
c                end if

      		rir=1.d0/ri
      		rjr=1.d0/rj

      		dxir=dxi*rir
      		dyir=dyi*rir
      		dzir=dzi*rir
      		dxjr=dxj*rjr
      		dyjr=dyj*rjr
      		dzjr=dzj*rjr

      		cst=dxir*dxjr+dyir*dyjr+dzir*dzjr
		if (cst.gt.1.d0) cst = 1.d0
		if (cst.lt.-1.d0) cst = -1.d0

c      		at=dacos(cst)
c---Now: e=k(cst-coseq)*(cst-coseq)=df*da This is angle energy!
c---Derivative: de/dx=2k(cst-coseq)*dcst/dx
c
c
      		da=cst-dcos(angleq(ith))
      		df=kangl(ith)*da
      		e = df*da
      		e_theta = e_theta + e
      		df = df+df
c---Now df=2k(cst-coseq); The only left is the derivative dcst/dx
c		at = dsin(at)
c		if (at.lt.1.d-3) then
c		 at = 1.d-3
c		else if (at.lt.0 .and. at.gt.-1.d-3) then
c		 at = -1.d-3
c		end if
c      		df=-df/at

                tmpx1 = rir*(dxjr-cst*dxir)
                tmpy1 = rir*(dyjr-cst*dyir)
                tmpz1 = rir*(dzjr-cst*dzir)

                tmpx2 = rjr*(dxir-cst*dxjr)
                tmpy2 = rjr*(dyir-cst*dyjr)
                tmpz2 = rjr*(dzir-cst*dzjr)

                dpot(1,i) = dpot(1,i) + df*tmpx1
                dpot(2,i) = dpot(2,i) + df*tmpy1
                dpot(3,i) = dpot(3,i) + df*tmpz1

                dpot(1,j) = dpot(1,j) - df*(tmpx1+tmpx2)
                dpot(2,j) = dpot(2,j) - df*(tmpy1+tmpy2)
                dpot(3,j) = dpot(3,j) - df*(tmpz1+tmpz2)

                dpot(1,k) = dpot(1,k) + df*tmpx2
                dpot(2,k) = dpot(2,k) + df*tmpy2
                dpot(3,k) = dpot(3,k) + df*tmpz2


100   continue

      return	

      end











