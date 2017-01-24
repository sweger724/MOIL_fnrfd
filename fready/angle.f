      function  angle(i,j,k)
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'

c calculate an angle between particles i,j,k

      double precision dxi,dyi,dzi,dxj,dyj,dzj,ri2,rj2
      double precision angle
      double precision ri,rj,rir,rjr,dxir,dyir,dzir,dxjr,dyjr,dzjr
      double precision cst
      integer i,j,k

                dxi = coor2(1,i) - coor2(1,j)
                dyi = coor2(2,i) - coor2(2,j)
                dzi = coor2(3,i) - coor2(3,j)
                dxj = coor2(1,k) - coor2(1,j)
                dyj = coor2(2,k) - coor2(2,j)
                dzj = coor2(3,k) - coor2(3,j)

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

                angle=dacos(cst)

      return

      end
