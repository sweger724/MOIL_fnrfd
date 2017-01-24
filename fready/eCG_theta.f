      subroutine eCG_theta()
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/FREADY.BLOCK'
      include 'COMMON/CONVERT.BLOCK'

c calculate angle energies and forces

c note that vector.block contains temporary vectors used for
c vectorization.

      double precision dxi,dyi,dzi,dxj,dyj,dzj,ri2,rj2
      double precision ri,rj,rir,rjr,dxir,dyir,dzir,dxjr,dyjr,dzjr
      double precision cst,at,da,df,e,tmpx1,tmpy1,tmpz1
      double precision tmpx2,tmpy2,tmpz2
      integer ith,i,j,k, T


c initialize e_theta
      e_theta=0.0

c initialize loop over angles in ichunk chunks
      do 100 ith = 1 , nangl

                i=iangl1(ith)
                j=iangl2(ith)
                k=iangl3(ith)

                dxi = coor(1,i) - coor(1,j)
                dyi = coor(2,i) - coor(2,j)
                dzi = coor(3,i) - coor(3,j)
                dxj = coor(1,k) - coor(1,j)
                dyj = coor(2,k) - coor(2,j)
                dzj = coor(3,k) - coor(3,j)
C@
C               write(*,*)' i j k coor(i) coor(j) coor(k) '
C               write(*,*)i,j,k
C               write(*,*)(coor(T,i),T=1,3)
C               write(*,*)(coor(T,j),T=1,3)
C               write(*,*)(coor(T,k),T=1,3)
C               write(*,*)' angleq ',angleq(ith),kangl(ith)

                ri2=dxi*dxi+dyi*dyi+dzi*dzi
                rj2=dxj*dxj+dyj*dyj+dzj*dzj
                ri=dsqrt(ri2)
                rj=dsqrt(rj2)

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

                at=dacos(cst)
                T = anglType(ith)
C               write(6,*)"Angle style: ",T,pi180*at
               if (T .le. 60 ) then
                  call eCG_WellAnglePot(at*pi180,T,E,df)
                  if (E .gt. E_CG_max) then
                 write (6,*) "Residue: ",poimon(i),poimon(j),poimon(k),T
                  endif
               else
                 E =0.d0
                 df = 0.d0
               endif

               e_theta = e_theta + E
               df = df * pi180

C               write(*,*)' Angle energy: ',e, at
C               df = df+df
                at = dsin(at)
                if (at.lt.1.d-6 .and. at.gt.0.d0) then
                 at = 1.d-6
                else if (at.lt.0.d0 .and. at.gt.-1.d-6) then
                 at = -1.d-6
                end if
                df=-df/at

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

C      write(6,*)"Average angle is:", SumAt/nangl

      return    

      end
