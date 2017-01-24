      subroutine ebond_specl(tmpnmb)
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONSPECL1.BLOCK'
      include 'COMMON/CONSPECL2.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SPECL.BLOCK'

c calculate bond energies and forces.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision dfx(maxmorsb),dfy(maxmorsb),dfz(maxmorsb)
      double precision tempdx1(maxspcl),tempdy1(maxspcl)
      double precision tempdz1(maxspcl)
      double precision xtemp1(maxspcl),ytemp1(maxspcl),ztemp1(maxspcl)
      double precision tempdx2(maxspcl),tempdy2(maxspcl)
      double precision tempdz2(maxspcl)
      double precision xtemp2(maxspcl),ytemp2(maxspcl),ztemp2(maxspcl)

      double precision rx,ry,rz,r2,s,r,db,df,e
      integer mm,i,j,n,ii,jj,tmpnmb
      integer iv,rv,nn,initial
      integer factr(maxspcl)

c e_bond=total bond energy 
c rx,ry,rz = distance squared between particles in that axis
c r2 = distance squared between particles
c s = distance between particles
c db,df temporary variables
c e = bond energy (non-acumulated)
c also used : ichunk for vectorizacion length
c xtemp,ytemp,ztemp for vectorization purposes



c initialize eb and emb
      n = tmpnmb
      se_bond1=0.d0
      se_bond2=0.d0

c calculation of bond-energies in unligated case
      snb1(0) = 0
      initial = snb1(n-1) + 1
c initialize loop over bonds in ichunk chunks
        do 100 iv=initial,snb1(n),ichunk

c check for chunk size, since the last one is (usually) smaller than
c ichunk
	rv=min(ichunk,snb1(n)-iv+1)

	do 200 nn=1,rv

c maintain the pointer to the "real" vector
		mm=iv-1+nn

                ii = sib11(mm)
                jj = sib12(mm)
      		i=newpoit(ii)
      		j=newpoit(jj)
      		rx=coor(1,i)-coor(1,j)
      		ry=coor(2,i)-coor(2,j)
      		rz=coor(3,i)-coor(3,j)
      		r2=rx*rx + ry*ry + rz*rz
      		s=dsqrt(r2)
      		r=2.d0/s
      		db=s-sreq1(mm)
      		df=skbond1(mm)*db
		e=df*db
		se_bond1 = se_bond1 + e
		df=df*r
      		xtemp1(nn)=rx*df
      		ytemp1(nn)=ry*df
      		ztemp1(nn)=rz*df
200            continue

               do 300 nn=1,rv

                mm=nn-1+iv
                ii = sib11(mm)
                jj = sib12(mm)
      		i=newpoit(ii)
      		j=newpoit(jj)

		dx1(ii)=dx1(ii)+xtemp1(nn)
		dy1(ii)=dy1(ii)+ytemp1(nn)
		dz1(ii)=dz1(ii)+ztemp1(nn)
		dx1(jj)=dx1(jj)-xtemp1(nn)
		dy1(jj)=dy1(jj)-ytemp1(nn)
		dz1(jj)=dz1(jj)-ztemp1(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

		 xtemp1(nn) = f(n) *xtemp1(nn)
		 ytemp1(nn) = f(n) *ytemp1(nn)
		 ztemp1(nn) = f(n) *ztemp1(nn)
		 dfx(n)=tempdfx(n)*e
		 dfy(n)=tempdfy(n)*e
		 dfz(n)=tempdfz(n)*e
		 tempdx1(ii) = xtemp1(nn) + factr(ii) * dfx(n)
		 tempdy1(ii) = ytemp1(nn) + factr(ii) * dfy(n)
		 tempdz1(ii) = ztemp1(nn) + factr(ii) * dfz(n)
		 tempdx1(jj) = -xtemp1(nn) + factr(jj) * dfx(n)
		 tempdy1(jj) = -ytemp1(nn) + factr(jj) * dfy(n)
		 tempdz1(jj) = -ztemp1(nn) + factr(jj) * dfz(n)
7               continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     		sdx(ii)=sdx(ii)+tempdx1(ii)
      		sdy(ii)=sdy(ii)+tempdy1(ii)
      		sdz(ii)=sdz(ii)+tempdz1(ii)
      		sdx(jj)=sdx(jj)+tempdx1(jj)
      		sdy(jj)=sdy(jj)+tempdy1(jj)
      		sdz(jj)=sdz(jj)+tempdz1(jj)
c =======================================
300	continue

100   continue

c calculation of the bond-energy in ligated case
      snb2(0) = 0
      initial = snb2(n-1) + 1
      do 1000 iv=initial,snb2(n),ichunk

c check for chunk size, since the last one is (usually) smaller than
c ichunk
	rv=min(ichunk,snb2(n)-iv+1)

	do 2000 nn=1,rv

c maintain the pointer to the "real" vector
		mm=iv-1+nn

		ii = sib21(mm)
		jj = sib22(mm)
      		i=newpoit(ii)
      		j=newpoit(jj)
      		rx=coor(1,i)-coor(1,j)
      		ry=coor(2,i)-coor(2,j)
      		rz=coor(3,i)-coor(3,j)
      		r2=rx*rx + ry*ry + rz*rz
      		s=dsqrt(r2)
      		r=2.d0/s
      		db=s-sreq2(mm)
      		df=skbond2(mm)*db
		e=df*db
		se_bond2 = se_bond2 + e
		df=df*r
      		xtemp2(nn)=rx*df
      		ytemp2(nn)=ry*df
      		ztemp2(nn)=rz*df

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

		  dfx(n)=tempdfx(n)*e
		  dfy(n)=tempdfy(n)*e
		  dfz(n)=tempdfz(n)*e
		  xtemp2(nn)= (1-f(n)) * xtemp2(nn)
		  ytemp2(nn)= (1-f(n)) * ytemp2(nn)
		  ztemp2(nn)= (1-f(n)) * ztemp2(nn)
		  tempdx2(ii) = xtemp2(nn) - factr(ii) * dfx(n)
		  tempdy2(ii) = ytemp2(nn) - factr(ii) * dfy(n)
		  tempdz2(ii) = ztemp2(nn) - factr(ii) * dfz(n)
		  tempdx2(jj) = -xtemp2(nn) - factr(jj) * dfx(n)
		  tempdy2(jj) = -ytemp2(nn) - factr(jj) * dfy(n)
		  tempdz2(jj) = -ztemp2(nn) - factr(jj) * dfz(n)
     		sdx(ii)=sdx(ii)+tempdx2(ii)
      		sdy(ii)=sdy(ii)+tempdy2(ii)
      		sdz(ii)=sdz(ii)+tempdz2(ii)
      		sdx(jj)=sdx(jj)+tempdx2(jj)
      		sdy(jj)=sdy(jj)+tempdy2(jj)
      		sdz(jj)=sdz(jj)+tempdz2(jj)
c =======================================
2000	continue

1000   continue

      return
      end
