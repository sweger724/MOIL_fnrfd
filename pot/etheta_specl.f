      subroutine etheta_specl(tmpnmb)
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONSPECL1.BLOCK'
      include 'COMMON/CONSPECL2.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SPECL.BLOCK'

c calculate angle energies and forces

c note that vector.block contains temporary vectors used for
c vectorization.

      double precision tempdx(maxspcl),tempdy(maxspcl)
      double precision tempdz(maxspcl)
      double precision dfx,dfy,dfz
      double precision dxi,dyi,dzi,dxj,dyj,dzj,ri2,rj2
      double precision ri,rj,rir,rjr,dxir,dyir,dzir,dxjr,dyjr,dzjr
      double precision cst,at,da,df,e
      double precision tmpx1,tmpx2,tmpy1,tmpy2,tmpz1,tmpz2
      integer ith,i,j,k,n,ii,jj,kk
      integer initial,tmpnmb
      integer factr(maxspcl)

c initialize e_theta
      se_theta1=0.0d0
      se_theta2=0.0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      n = tmpnmb
      snang1(0) = 0
      initial = snang1(n-1) + 1
c initialize loop over angles in ichunk chunks
      do 100 ith=initial,snang1(n)

                ii=siangl11(ith)
                jj=siangl12(ith)
                kk=siangl13(ith)
                i = newpoit(ii)
                j = newpoit(jj)
                k = newpoit(kk)


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
                da=at-sangleq1(ith)
                df=skangl1(ith)*da
                e = df*da
                se_theta1 = se_theta1 + e
                df = df+df
                at = dsin(at)
                if (at.lt.1.d-3)  then
                 at = 1.d-3
                else if (at.lt.0 .and. at.gt.-1.d-3) then
                 at = -1.d-3
                end if

                df=-df/at

                tmpx1 = df* rir * (dxjr - cst * dxir)
                tmpx2 = df* rjr * (dxir - cst * dxjr)
                tmpy1 = df* rir * (dyjr - cst * dyir)
                tmpy2 = df* rjr * (dyir - cst * dyjr)
                tmpz1 = df* rir * (dzjr - cst * dzir)
                tmpz2 = df* rjr * (dzir - cst * dzjr)

                dx1(ii) = dx1(ii) + tmpx1
                dx1(kk) = dx1(kk) + tmpx2
                dx1(jj) = dx1(jj) - (tmpx1 + tmpx2)
                dy1(ii) = dy1(ii) + tmpy1
                dy1(kk) = dy1(kk) + tmpy2
                dy1(jj) = dy1(jj) - (tmpy1 + tmpy2)
                dz1(ii) = dz1(ii) + tmpz1
                dz1(kk) = dz1(kk) + tmpz2
                dz1(jj) = dz1(jj) - (tmpz1 + tmpz2)
c =======================================================
                 if (i .eq. imb1(n)) then
                  factr(ii) = 1
                 else if (i .eq. imb2(n)) then
                  factr(ii) = -1
                 else
                  factr(ii) = 0
                 end if
                 if (j .eq. imb1(n)) then
                  factr(jj) = 1
                 else if (j .eq. imb2(n)) then
                  factr(jj) = -1
                 else
                  factr(jj) = 0
                 end if
                 if (k .eq. imb1(n)) then
                  factr(kk) = 1
                 else if (k .eq. imb2(n)) then
                  factr(kk) = -1
                 else
                  factr(kk) = 0
                 end if
 
                dfx = tempdfx(n) * e
                dfy = tempdfy(n) * e
                dfz = tempdfz(n) * e

                tempdx(ii) = f(n) * tmpx1   + factr(ii)*dfx
                tempdx(kk) = f(n) * tmpx2  + factr(kk)*dfx
                tempdx(jj) = f(n) * (-tmpx1 - tmpx2)
     $          + factr(jj) * dfx
                tempdy(ii) = f(n) * tmpy1   + factr(ii)*dfy
                tempdy(kk) = f(n) * tmpy2  + factr(kk)*dfy
                tempdy(jj) = f(n) * (-tmpy1 - tmpy2)
     $          + factr(jj) * dfy
                tempdz(ii) = f(n) * tmpz1   + factr(ii)*dfz
                tempdz(kk) = f(n) * tmpz2  + factr(kk)*dfz
                tempdz(jj) = f(n) * (-tmpz1 - tmpz2)
     $          + factr(jj) * dfz
                
                sdx(ii) = sdx(ii) + tempdx(ii)
                sdx(kk) = sdx(kk) + tempdx(kk)
                sdx(jj) = sdx(jj) + tempdx(jj)
                sdy(ii) = sdy(ii) + tempdy(ii)
                sdy(kk) = sdy(kk) + tempdy(kk)
                sdy(jj) = sdy(jj) + tempdy(jj)
                sdz(ii) = sdz(ii) + tempdz(ii)
                sdz(kk) = sdz(kk) + tempdz(kk)
                sdz(jj) = sdz(jj) + tempdz(jj)

100   continue
cccccccccccccccccccccccccccccccccccccccccccccc
      snang2(0) = 0
      initial = snang2(n-1) + 1
cccccccccccccccccccccccccccccccccccccccccccccccc
c initialize loop over angles in ichunk chunks
      do 1000 ith=initial,snang2(n)



                ii = siangl21(ith)
                jj = siangl22(ith)
                kk = siangl23(ith)
                i = newpoit(ii)
                j = newpoit(jj)
                k = newpoit(kk)


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
                da=at-sangleq2(ith)
                df=skangl2(ith)*da
                e = df*da
                se_theta2 = se_theta2 + e
                df = df+df
                at = dsin(at)
                if (at.lt.1.d-3) then
                 at = 1.d-3
                else if (at.lt.0 .and. at.gt.-1.d-3) then
                 at = -1.d-3
                end if

                df=-df/at

                tmpx1 = df * rir * (dxjr - cst * dxir)
                tmpx2 = df * rjr * (dxir - cst * dxjr)
                tmpy1 = df * rir * (dyjr - cst * dyir)
                tmpy2 = df * rjr * (dyir - cst * dyjr)
                tmpz1 = df * rir * (dzjr - cst * dzir)
                tmpz2 = df * rjr * (dzir - cst * dzjr)


c =======================================================
                 if (i .eq. imb1(n)) then
                  factr(ii)=1
                 else if (i .eq. imb2(n)) then
                  factr(ii)=-1
                 else
                  factr(ii)=0
                 end if
                 if (j .eq. imb1(n)) then
                  factr(jj)=1
                 else if (j .eq. imb2(n)) then
                  factr(jj)=-1
                 else
                  factr(jj)=0
                 end if
                 if (k .eq. imb1(n)) then
                  factr(kk)=1
                 else if (k .eq. imb2(n)) then
                  factr(kk)=-1
                 else
                  factr(kk)=0
                 end if

                 dfx = tempdfx(n) * e
                 dfy = tempdfy(n) * e
                 dfz = tempdfz(n) * e
                 tempdx(ii) = (1-f(n)) * tmpx1   - factr(ii)*dfx
                 tempdx(kk) = (1-f(n)) * tmpx2  - factr(kk)*dfx
                 tempdx(jj) = (1-f(n)) * (-tmpx1 - tmpx2)
     $           - factr(jj) * dfx
                 tempdy(ii) = (1-f(n))  * tmpy1   - factr(ii)*dfy
                 tempdy(kk) = (1-f(n))  * tmpy2  - factr(kk)*dfy
                 tempdy(jj) = (1-f(n))  * (-tmpy1 - tmpy2)
     $           - factr(jj) * dfy
                 tempdz(ii) = (1-f(n)) * tmpz1  - factr(ii)*dfz
                 tempdz(kk) = (1-f(n)) * tmpz2 - factr(kk)*dfz
                 tempdz(jj) = (1-f(n)) * (-tmpz1 - tmpz2)
     $           - factr(jj) *dfz
                
                sdx(ii) = sdx(ii) + tempdx(ii)
                sdx(kk) = sdx(kk) + tempdx(kk)
                sdx(jj) = sdx(jj) + tempdx(jj)
                sdy(ii) = sdy(ii) + tempdy(ii)
                sdy(kk) = sdy(kk) + tempdy(kk)
                sdy(jj) = sdy(jj) + tempdy(jj)
                sdz(ii) = sdz(ii) + tempdz(ii)
                sdz(kk) = sdz(kk) + tempdz(kk)
                sdz(jj) = sdz(jj) + tempdz(jj)

1000    continue

c      print*,'se_theta=',se_theta
      return    

      end
