	subroutine prd2v(unit)
c
c	printing out the second derivative matrix
c	by constructing first a full 3N*3N matrix
c
	implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SCNDRV.BLOCK'
	include 'COMMON/SC2.BLOCK'

	call prd2v1(d2vec,unit)

	return
	end

	subroutine prd2v1(d2v,unit)
	implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SCNDRV.BLOCK'

	double precision d2v(3*npt,*)

	integer unit
	integer i,j,k,l,m,n,jstart1,jstart2,jstart3,jstart4
	integer ii,jj,kk,ll,mm,nn,kkk,jend1,jend2,jend3,jend4
	integer jjj,ith,iphi,iiphi,iith
	
	do 1 i=1,3*npt
	 do 1 j=1,3*npt
	  d2v(i,j) = 0.d0
1	continue

c
c copy from the special packing of the non bonded second derivative 
c matrix to the regular 3N*3N second derivative matrix
c
	jstart1 = 1
	jstart2 = 1
	jstart3 = 1
	jstart4 = 1
	do 41 i=1,npt
		j = 6*(i-1)+1
		k = 3*(i-1)+1
		d2v(k,k)     = diag(j)
		d2v(k+1,k)   = diag(j+1)
		d2v(k,k+1)   = diag(j+1)
		d2v(k+2,k)   = diag(j+2)
		d2v(k,k+2)   = diag(j+2)
		d2v(k+1,k+1) = diag(j+3)
		d2v(k+1,k+2) = diag(j+4)
		d2v(k+2,k+1) = diag(j+4)
		d2v(k+2,k+2) = diag(j+5)
		if (eelyes.or.evdyes) then
		jend1 = point1(i)
		do 23 jj=jstart1,jend1
		 kk = list1(jj)
		 kkk=  3*(kk-1) + 1
		 jjj=  6*(jj-1) + 1
c element xi,xj which is the same as xj,xi (i,j, particle indices)
c
		 d2v(k,kkk)     = d2nb_1(jjj)
		 d2v(kkk,k)     = d2nb_1(jjj)
c element xi,yj which is the same as yj,xi , xj,yi and yi,xj
c the last two are the same only for function of distances
c
		 d2v(k+1,kkk)   = d2nb_1(jjj+1)
		 d2v(kkk+1,k)   = d2nb_1(jjj+1)
		 d2v(kkk,k+1)   = d2nb_1(jjj+1)
		 d2v(k,kkk+1)   = d2nb_1(jjj+1)
c elements zi,xj xj,zi zj,xi xi,zj
		 d2v(k+2,kkk)   = d2nb_1(jjj+2)
		 d2v(kkk,k+2)   = d2nb_1(jjj+2)
		 d2v(kkk+2,k)   = d2nb_1(jjj+2)
		 d2v(k,kkk+2)   = d2nb_1(jjj+2)
c element yi,yj and yj,yi
		 d2v(k+1,kkk+1) = d2nb_1(jjj+3)
		 d2v(kkk+1,k+1) = d2nb_1(jjj+3)
c elements yi,zj zj,yi yj,zi zi,yj
		 d2v(k+1,kkk+2) = d2nb_1(jjj+4)
		 d2v(kkk+2,k+1) = d2nb_1(jjj+4)
		 d2v(kkk+1,k+2) = d2nb_1(jjj+4)
		 d2v(k+2,kkk+1) = d2nb_1(jjj+4)
c elements zi,zj zj,zi
		 d2v(k+2,kkk+2) = d2nb_1(jjj+5)
		 d2v(kkk+2,k+2) = d2nb_1(jjj+5)
23		continue
		jstart1 = jend1 + 1
	        jend2 = point2(i)
		do 24 jj=jstart2,jend2
		 kk = list2(jj)
		 kkk=  3*(kk-1) + 1
		 jjj=  6*(jj-1) + 1
c element xi,xj which is the same as xj,xi (i,j, particle indices)
c
		 d2v(k,kkk)     = d2nb_2(jjj)
		 d2v(kkk,k)     = d2nb_2(jjj)
c element xi,yj which is the same as yj,xi , xj,yi and yi,xj
c the last two are the same only for function of distances
c
		 d2v(k+1,kkk)   = d2nb_2(jjj+1)
		 d2v(kkk+1,k)   = d2nb_2(jjj+1)
		 d2v(kkk,k+1)   = d2nb_2(jjj+1)
		 d2v(k,kkk+1)   = d2nb_2(jjj+1)
c elements zi,xj xj,zi zj,xi xi,zj
		 d2v(k+2,kkk)   = d2nb_2(jjj+2)
		 d2v(kkk,k+2)   = d2nb_2(jjj+2)
		 d2v(kkk+2,k)   = d2nb_2(jjj+2)
		 d2v(k,kkk+2)   = d2nb_2(jjj+2)
c element yi,yj and yj,yi
		 d2v(k+1,kkk+1) = d2nb_2(jjj+3)
		 d2v(kkk+1,k+1) = d2nb_2(jjj+3)
c elements yi,zj zj,yi yj,zi zi,yj
		 d2v(k+1,kkk+2) = d2nb_2(jjj+4)
		 d2v(kkk+2,k+1) = d2nb_2(jjj+4)
		 d2v(kkk+1,k+2) = d2nb_2(jjj+4)
		 d2v(k+2,kkk+1) = d2nb_2(jjj+4)
c elements zi,zj zj,zi
		 d2v(k+2,kkk+2) = d2nb_2(jjj+5)
		 d2v(kkk+2,k+2) = d2nb_2(jjj+5)
24		continue
		jstart2 = jend2 + 1
	        jend3 = point3(i)
		do 25 jj=jstart3,jend3
		 kk = list3(jj)
		 kkk=  3*(kk-1) + 1
		 jjj=  6*(jj-1) + 1
c element xi,xj which is the same as xj,xi (i,j, particle indices)
c
		 d2v(k,kkk)     = d2nb_3(jjj)
		 d2v(kkk,k)     = d2nb_3(jjj)
c element xi,yj which is the same as yj,xi , xj,yi and yi,xj
c the last two are the same only for function of distances
c
		 d2v(k+1,kkk)   = d2nb_3(jjj+1)
		 d2v(kkk+1,k)   = d2nb_3(jjj+1)
		 d2v(kkk,k+1)   = d2nb_3(jjj+1)
		 d2v(k,kkk+1)   = d2nb_3(jjj+1)
c elements zi,xj xj,zi zj,xi xi,zj
		 d2v(k+2,kkk)   = d2nb_3(jjj+2)
		 d2v(kkk,k+2)   = d2nb_3(jjj+2)
		 d2v(kkk+2,k)   = d2nb_3(jjj+2)
		 d2v(k,kkk+2)   = d2nb_3(jjj+2)
c element yi,yj and yj,yi
		 d2v(k+1,kkk+1) = d2nb_3(jjj+3)
		 d2v(kkk+1,k+1) = d2nb_3(jjj+3)
c elements yi,zj zj,yi yj,zi zi,yj
		 d2v(k+1,kkk+2) = d2nb_3(jjj+4)
		 d2v(kkk+2,k+1) = d2nb_3(jjj+4)
		 d2v(kkk+1,k+2) = d2nb_3(jjj+4)
		 d2v(k+2,kkk+1) = d2nb_3(jjj+4)
c elements zi,zj zj,zi
		 d2v(k+2,kkk+2) = d2nb_3(jjj+5)
		 d2v(kkk+2,k+2) = d2nb_3(jjj+5)
25		continue
		jstart3 = jend3 + 1
	        jend4 = spec1(i)
		do 26 jj=jstart4,jend4
		 kk = spec2(jj)
		 kkk=  3*(kk-1) + 1
		 jjj=  6*(jj-1) + 1
c element xi,xj which is the same as xj,xi (i,j, particle indices)
c
		 d2v(k,kkk)     = d2v(k,kkk) + d2spec(jjj)
		 d2v(kkk,k)     = d2v(kkk,k) + d2spec(jjj)
c element xi,yj which is the same as yj,xi , xj,yi and yi,xj
c the last two are the same only for function of distances
c
		 d2v(k+1,kkk)   = d2v(k+1,kkk) + d2spec(jjj+1)
		 d2v(kkk+1,k)   = d2v(kkk+1,k) + d2spec(jjj+1)
		 d2v(kkk,k+1)   = d2v(kkk,k+1) + d2spec(jjj+1)
		 d2v(k,kkk+1)   = d2v(k,kkk+1) + d2spec(jjj+1)
c elements zi,xj xj,zi zj,xi xi,zj
		 d2v(k+2,kkk)   = d2v(k+2,kkk) + d2spec(jjj+2)
		 d2v(kkk,k+2)   = d2v(kkk,k+2) + d2spec(jjj+2)
		 d2v(kkk+2,k)   = d2v(kkk+2,k) + d2spec(jjj+2)
		 d2v(k,kkk+2)   = d2v(k,kkk+2) + d2spec(jjj+2)
c element yi,yj and yj,yi
		 d2v(k+1,kkk+1) = d2v(k+1,kkk+1) + d2spec(jjj+3)
		 d2v(kkk+1,k+1) = d2v(kkk+1,k+1) + d2spec(jjj+3)
c elements yi,zj zj,yi yj,zi zi,yj
		 d2v(k+1,kkk+2) = d2v(k+1,kkk+2) + d2spec(jjj+4)
		 d2v(kkk+2,k+1) = d2v(kkk+2,k+1) + d2spec(jjj+4)
		 d2v(kkk+1,k+2) = d2v(kkk+1,k+2) + d2spec(jjj+4)
		 d2v(k+2,kkk+1) = d2v(k+2,kkk+1) + d2spec(jjj+4)
c elements zi,zj zj,zi
		 d2v(k+2,kkk+2) = d2v(k+2,kkk+2) + d2spec(jjj+5)
		 d2v(kkk+2,k+2) = d2v(kkk+2,k+2) + d2spec(jjj+5)
26		continue
		jstart4 = jend4 + 1
		end if
41	continue
	if (ebyes) then
	do 42 i=1,nb
		j  = 3*(ib1(i)-1)+1
		k  = 3*(ib2(i)-1)+1
		kk = 6*(i-1) + 1
c xi xj ; xj xi
		d2v(j,k)     = d2v(j,k)     + d2bond(kk)
		d2v(k,j)     = d2v(k,j)     + d2bond(kk)
c yi xj ; xj yi ; yj xi ; xi yj
		d2v(j+1,k)   = d2v(j+1,k)   + d2bond(kk+1)
		d2v(k,j+1)   = d2v(k,j+1)   + d2bond(kk+1)
		d2v(j,k+1)   = d2v(j,k+1)   + d2bond(kk+1)
		d2v(k+1,j)   = d2v(k+1,j)   + d2bond(kk+1)
c zi xj ; xi zj ; xj zi ; zj xi
		d2v(j+2,k)   = d2v(j+2,k)   + d2bond(kk+2)
		d2v(j,k+2)   = d2v(j,k+2)   + d2bond(kk+2)
		d2v(k+2,j)   = d2v(k+2,j)   + d2bond(kk+2)
		d2v(k,j+2)   = d2v(k,j+2)   + d2bond(kk+2)
c yi yj ; yj yi
		d2v(j+1,k+1) = d2v(j+1,k+1) + d2bond(kk+3)
		d2v(k+1,j+1) = d2v(k+1,j+1) + d2bond(kk+3)
c zi yj ; yi zj ; zj yi ; yj zi
		d2v(j+2,k+1) = d2v(j+2,k+1) + d2bond(kk+4)
		d2v(j+1,k+2) = d2v(j+1,k+2) + d2bond(kk+4)
		d2v(k+1,j+2) = d2v(k+1,j+2) + d2bond(kk+4)
		d2v(k+2,j+1) = d2v(k+2,j+1) + d2bond(kk+4)
c zi zj ; zj zi
		d2v(j+2,k+2) = d2v(j+2,k+2) + d2bond(kk+5)
		d2v(k+2,j+2) = d2v(k+2,j+2) + d2bond(kk+5)
42	continue
	end if
	if (ethyes) then
	do 43 ith=1,nangl
		i=iangl1(ith)
		j=iangl2(ith)
		k=iangl3(ith)
		ii = 3*(i-1)+1
		jj = 3*(j-1)+1
		kk = 3*(k-1)+1
		iith = 27*(ith-1)+1
c pair i,j
		d2v(ii,jj)     =  d2v(ii,jj) + d2theta(iith)
		d2v(jj,ii)     =  d2v(jj,ii) + d2theta(iith)
		d2v(ii,jj+1)   =  d2v(ii,jj+1)   + d2theta(iith+1)
		d2v(jj+1,ii)   =  d2v(jj+1,ii)   + d2theta(iith+1)
		d2v(ii,jj+2)   =  d2v(ii,jj+2)   + d2theta(iith+2)
		d2v(jj+2,ii)   =  d2v(jj+2,ii)   + d2theta(iith+2)
		d2v(ii+1,jj)   =  d2v(ii+1,jj)   + d2theta(iith+3)
		d2v(jj,ii+1)   =  d2v(jj,ii+1)   + d2theta(iith+3)
		d2v(ii+1,jj+1) =  d2v(ii+1,jj+1) + d2theta(iith+4)
		d2v(jj+1,ii+1) =  d2v(jj+1,ii+1) + d2theta(iith+4)
		d2v(ii+1,jj+2) =  d2v(ii+1,jj+2) + d2theta(iith+5)
		d2v(jj+2,ii+1) =  d2v(jj+2,ii+1) + d2theta(iith+5)
		d2v(ii+2,jj)   =  d2v(ii+2,jj)   + d2theta(iith+6)
		d2v(jj,ii+2)   =  d2v(jj,ii+2)   + d2theta(iith+6)
		d2v(ii+2,jj+1) =  d2v(ii+2,jj+1) + d2theta(iith+7)
		d2v(jj+1,ii+2) =  d2v(jj+1,ii+2) + d2theta(iith+7)
		d2v(ii+2,jj+2) =  d2v(ii+2,jj+2) + d2theta(iith+8)
		d2v(jj+2,ii+2) =  d2v(jj+2,ii+2) + d2theta(iith+8)

c pair i,k
		d2v(ii,kk)     =  d2v(ii,kk)     + d2theta(iith+9)
		d2v(kk,ii)     =  d2v(kk,ii)     + d2theta(iith+9)
		d2v(ii,kk+1)   =  d2v(ii,kk+1)   + d2theta(iith+10)
		d2v(kk+1,ii)   =  d2v(kk+1,ii)   + d2theta(iith+10)
		d2v(ii,kk+2)   =  d2v(ii,kk+2)   + d2theta(iith+11)
		d2v(kk+2,ii)   =  d2v(kk+2,ii)   + d2theta(iith+11)
		d2v(ii+1,kk)   =  d2v(ii+1,kk)   + d2theta(iith+12)
		d2v(kk,ii+1)   =  d2v(kk,ii+1)   + d2theta(iith+12)
		d2v(ii+1,kk+1) =  d2v(ii+1,kk+1) + d2theta(iith+13)
		d2v(kk+1,ii+1) =  d2v(kk+1,ii+1) + d2theta(iith+13)
		d2v(ii+1,kk+2) =  d2v(ii+1,kk+2) + d2theta(iith+14)
		d2v(kk+2,ii+1) =  d2v(kk+2,ii+1) + d2theta(iith+14)
		d2v(ii+2,kk)   =  d2v(ii+2,kk)   + d2theta(iith+15)
		d2v(kk,ii+2)   =  d2v(kk,ii+2)   + d2theta(iith+15)
		d2v(ii+2,kk+1) =  d2v(ii+2,kk+1) + d2theta(iith+16)
		d2v(kk+1,ii+2) =  d2v(kk+1,ii+2) + d2theta(iith+16)
		d2v(ii+2,kk+2) =  d2v(ii+2,kk+2) + d2theta(iith+17)
		d2v(kk+2,ii+2) =  d2v(kk+2,ii+2) + d2theta(iith+17)

c pair j,k
		d2v(jj,kk)     =  d2v(jj,kk)     + d2theta(iith+18)
		d2v(kk,jj)     =  d2v(kk,jj)     + d2theta(iith+18)
		d2v(jj,kk+1)   =  d2v(jj,kk+1)   + d2theta(iith+19)
		d2v(kk+1,jj)   =  d2v(kk+1,jj)   + d2theta(iith+19)
		d2v(jj,kk+2)   =  d2v(jj,kk+2)   + d2theta(iith+20)
		d2v(kk+2,jj)   =  d2v(kk+2,jj)   + d2theta(iith+20)
		d2v(jj+1,kk)   =  d2v(jj+1,kk)   + d2theta(iith+21)
		d2v(kk,jj+1)   =  d2v(kk,jj+1)   + d2theta(iith+21)
		d2v(jj+1,kk+1) =  d2v(jj+1,kk+1) + d2theta(iith+22)
		d2v(kk+1,jj+1) =  d2v(kk+1,jj+1) + d2theta(iith+22)
		d2v(jj+1,kk+2) =  d2v(jj+1,kk+2) + d2theta(iith+23)
		d2v(kk+2,jj+1) =  d2v(kk+2,jj+1) + d2theta(iith+23)
		d2v(jj+2,kk)   =  d2v(jj+2,kk)   + d2theta(iith+24)
		d2v(kk,jj+2)   =  d2v(kk,jj+2)   + d2theta(iith+24)
		d2v(jj+2,kk+1) =  d2v(jj+2,kk+1) + d2theta(iith+25)
		d2v(kk+1,jj+2) =  d2v(kk+1,jj+2) + d2theta(iith+25)
		d2v(jj+2,kk+2) =  d2v(jj+2,kk+2) + d2theta(iith+26)
		d2v(kk+2,jj+2) =  d2v(kk+2,jj+2) + d2theta(iith+26)

43	continue
	end if
	if (etoyes) then
	do 44 iphi=1,ntors
		k=itor1(iphi)
		l=itor2(iphi)
		m=itor3(iphi)
		n=itor4(iphi)
		kk=3*(k-1)+1
		ll=3*(l-1)+1
		mm=3*(m-1)+1
		nn=3*(n-1)+1
		iiphi=54*(iphi-1)+1
c pair kl
		d2v(kk,ll) = d2v(kk,ll) + d2phi(iiphi)
		d2v(ll,kk) = d2v(ll,kk) + d2phi(iiphi)
		d2v(kk+1,ll) = d2v(kk+1,ll) + d2phi(iiphi+1)
		d2v(ll,kk+1) = d2v(ll,kk+1) + d2phi(iiphi+1)
		d2v(kk+2,ll) = d2v(kk+2,ll) + d2phi(iiphi+2)
		d2v(ll,kk+2) = d2v(ll,kk+2) + d2phi(iiphi+2)
		d2v(kk,ll+1) = d2v(kk,ll+1) + d2phi(iiphi+3)
		d2v(ll+1,kk) = d2v(ll+1,kk) + d2phi(iiphi+3)
		d2v(kk+1,ll+1) = d2v(kk+1,ll+1) + d2phi(iiphi+4)
		d2v(ll+1,kk+1) = d2v(ll+1,kk+1) + d2phi(iiphi+4)
		d2v(kk+2,ll+1) = d2v(kk+2,ll+1) + d2phi(iiphi+5)
		d2v(ll+1,kk+2) = d2v(ll+1,kk+2) + d2phi(iiphi+5)
		d2v(kk,ll+2) = d2v(kk,ll+2) + d2phi(iiphi+6)
		d2v(ll+2,kk) = d2v(ll+2,kk) + d2phi(iiphi+6)
		d2v(kk+1,ll+2) = d2v(kk+1,ll+2) + d2phi(iiphi+7)
		d2v(ll+2,kk+1) = d2v(ll+2,kk+1) + d2phi(iiphi+7)
		d2v(kk+2,ll+2) = d2v(kk+2,ll+2) + d2phi(iiphi+8)
		d2v(ll+2,kk+2) = d2v(ll+2,kk+2) + d2phi(iiphi+8)

c pair km
		d2v(kk,mm) = d2v(kk,mm) + d2phi(iiphi+9)
		d2v(mm,kk) = d2v(mm,kk) + d2phi(iiphi+9)
		d2v(kk+1,mm) = d2v(kk+1,mm) + d2phi(iiphi+10)
		d2v(mm,kk+1) = d2v(mm,kk+1) + d2phi(iiphi+10)
		d2v(kk+2,mm) = d2v(kk+2,mm) + d2phi(iiphi+11)
		d2v(mm,kk+2) = d2v(mm,kk+2) + d2phi(iiphi+11)
		d2v(kk,mm+1) = d2v(kk,mm+1) + d2phi(iiphi+12)
		d2v(mm+1,kk) = d2v(mm+1,kk) + d2phi(iiphi+12)
		d2v(kk+1,mm+1) = d2v(kk+1,mm+1) + d2phi(iiphi+13)
		d2v(mm+1,kk+1) = d2v(mm+1,kk+1) + d2phi(iiphi+13)
		d2v(kk+2,mm+1) = d2v(kk+2,mm+1) + d2phi(iiphi+14)
		d2v(mm+1,kk+2) = d2v(mm+1,kk+2) + d2phi(iiphi+14)
		d2v(kk,mm+2) = d2v(kk,mm+2) + d2phi(iiphi+15)
		d2v(mm+2,kk) = d2v(mm+2,kk) + d2phi(iiphi+15)
		d2v(kk+1,mm+2) = d2v(kk+1,mm+2) + d2phi(iiphi+16)
		d2v(mm+2,kk+1) = d2v(mm+2,kk+1) + d2phi(iiphi+16)
		d2v(kk+2,mm+2) = d2v(kk+2,mm+2) + d2phi(iiphi+17)
		d2v(mm+2,kk+2) = d2v(mm+2,kk+2) + d2phi(iiphi+17)

c pair kn
		d2v(kk,nn) = d2v(kk,nn) + d2phi(iiphi+18)
		d2v(nn,kk) = d2v(nn,kk) + d2phi(iiphi+18)
		d2v(kk+1,nn) = d2v(kk+1,nn) + d2phi(iiphi+19)
		d2v(nn,kk+1) = d2v(nn,kk+1) + d2phi(iiphi+19)
		d2v(kk+2,nn) = d2v(kk+2,nn) + d2phi(iiphi+20)
		d2v(nn,kk+2) = d2v(nn,kk+2) + d2phi(iiphi+20)
		d2v(kk,nn+1) = d2v(kk,nn+1) + d2phi(iiphi+21)
		d2v(nn+1,kk) = d2v(nn+1,kk) + d2phi(iiphi+21)
		d2v(kk+1,nn+1) = d2v(kk+1,nn+1) + d2phi(iiphi+22)
		d2v(nn+1,kk+1) = d2v(nn+1,kk+1) + d2phi(iiphi+22)
		d2v(kk+2,nn+1) = d2v(kk+2,nn+1) + d2phi(iiphi+23)
		d2v(nn+1,kk+2) = d2v(nn+1,kk+2) + d2phi(iiphi+23)
		d2v(kk,nn+2) = d2v(kk,nn+2) + d2phi(iiphi+24)
		d2v(nn+2,kk) = d2v(nn+2,kk) + d2phi(iiphi+24)
		d2v(kk+1,nn+2) = d2v(kk+1,nn+2) + d2phi(iiphi+25)
		d2v(nn+2,kk+1) = d2v(nn+2,kk+1) + d2phi(iiphi+25)
		d2v(kk+2,nn+2) = d2v(kk+2,nn+2) + d2phi(iiphi+26)
		d2v(nn+2,kk+2) = d2v(nn+2,kk+2) + d2phi(iiphi+26)

c pair lm
		d2v(ll,mm) = d2v(ll,mm) + d2phi(iiphi+27)
		d2v(mm,ll) = d2v(mm,ll) + d2phi(iiphi+27)
		d2v(ll+1,mm) = d2v(ll+1,mm) + d2phi(iiphi+28)
		d2v(mm,ll+1) = d2v(mm,ll+1) + d2phi(iiphi+28)
		d2v(ll+2,mm) = d2v(ll+2,mm) + d2phi(iiphi+29)
		d2v(mm,ll+2) = d2v(mm,ll+2) + d2phi(iiphi+29)
		d2v(ll,mm+1) = d2v(ll,mm+1) + d2phi(iiphi+30)
		d2v(mm+1,ll) = d2v(mm+1,ll) + d2phi(iiphi+30)
		d2v(ll+1,mm+1) = d2v(ll+1,mm+1) + d2phi(iiphi+31)
		d2v(mm+1,ll+1) = d2v(mm+1,ll+1) + d2phi(iiphi+31)
		d2v(ll+2,mm+1) = d2v(ll+2,mm+1) + d2phi(iiphi+32)
		d2v(mm+1,ll+2) = d2v(mm+1,ll+2) + d2phi(iiphi+32)
		d2v(ll,mm+2) = d2v(ll,mm+2) + d2phi(iiphi+33)
		d2v(mm+2,ll) = d2v(mm+2,ll) + d2phi(iiphi+33)
		d2v(ll+1,mm+2) = d2v(ll+1,mm+2) + d2phi(iiphi+34)
		d2v(mm+2,ll+1) = d2v(mm+2,ll+1) + d2phi(iiphi+34)
		d2v(ll+2,mm+2) = d2v(ll+2,mm+2) + d2phi(iiphi+35)
		d2v(mm+2,ll+2) = d2v(mm+2,ll+2) + d2phi(iiphi+35)

c pair ln
		d2v(ll,nn) = d2v(ll,nn) + d2phi(iiphi+36)
		d2v(nn,ll) = d2v(nn,ll) + d2phi(iiphi+36)
		d2v(ll+1,nn) = d2v(ll+1,nn) + d2phi(iiphi+37)
		d2v(nn,ll+1) = d2v(nn,ll+1) + d2phi(iiphi+37)
		d2v(ll+2,nn) = d2v(ll+2,nn) + d2phi(iiphi+38)
		d2v(nn,ll+2) = d2v(nn,ll+2) + d2phi(iiphi+38)
		d2v(ll,nn+1) = d2v(ll,nn+1) + d2phi(iiphi+39)
		d2v(nn+1,ll) = d2v(nn+1,ll) + d2phi(iiphi+39)
		d2v(ll+1,nn+1) = d2v(ll+1,nn+1) + d2phi(iiphi+40)
		d2v(nn+1,ll+1) = d2v(nn+1,ll+1) + d2phi(iiphi+40)
		d2v(ll+2,nn+1) = d2v(ll+2,nn+1) + d2phi(iiphi+41)
		d2v(nn+1,ll+2) = d2v(nn+1,ll+2) + d2phi(iiphi+41)
		d2v(ll,nn+2) = d2v(ll,nn+2) + d2phi(iiphi+42)
		d2v(nn+2,ll) = d2v(nn+2,ll) + d2phi(iiphi+42)
		d2v(ll+1,nn+2) = d2v(ll+1,nn+2) + d2phi(iiphi+43)
		d2v(nn+2,ll+1) = d2v(nn+2,ll+1) + d2phi(iiphi+43)
		d2v(ll+2,nn+2) = d2v(ll+2,nn+2) + d2phi(iiphi+44)
		d2v(nn+2,ll+2) = d2v(nn+2,ll+2) + d2phi(iiphi+44)

c pair mn
		d2v(mm,nn) = d2v(mm,nn) + d2phi(iiphi+45)
		d2v(nn,mm) = d2v(nn,mm) + d2phi(iiphi+45)
		d2v(mm+1,nn) = d2v(mm+1,nn) + d2phi(iiphi+46)
		d2v(nn,mm+1) = d2v(nn,mm+1) + d2phi(iiphi+46)
		d2v(mm+2,nn) = d2v(mm+2,nn) + d2phi(iiphi+47)
		d2v(nn,mm+2) = d2v(nn,mm+2) + d2phi(iiphi+47)
		d2v(mm,nn+1) = d2v(mm,nn+1) + d2phi(iiphi+48)
		d2v(nn+1,mm) = d2v(nn+1,mm) + d2phi(iiphi+48)
		d2v(mm+1,nn+1) = d2v(mm+1,nn+1) + d2phi(iiphi+49)
		d2v(nn+1,mm+1) = d2v(nn+1,mm+1) + d2phi(iiphi+49)
		d2v(mm+2,nn+1) = d2v(mm+2,nn+1) + d2phi(iiphi+50)
		d2v(nn+1,mm+2) = d2v(nn+1,mm+2) + d2phi(iiphi+50)
		d2v(mm,nn+2) = d2v(mm,nn+2) + d2phi(iiphi+51)
		d2v(nn+2,mm) = d2v(nn+2,mm) + d2phi(iiphi+51)
		d2v(mm+1,nn+2) = d2v(mm+1,nn+2) + d2phi(iiphi+52)
		d2v(nn+2,mm+1) = d2v(nn+2,mm+1) + d2phi(iiphi+52)
		d2v(mm+2,nn+2) = d2v(mm+2,nn+2) + d2phi(iiphi+53)
		d2v(nn+2,mm+2) = d2v(nn+2,mm+2) + d2phi(iiphi+53)
44	continue
	end if

	if (eimyes) then
	do 45 iphi=1,nimp
		k=iimp1(iphi)
		l=iimp2(iphi)
		m=iimp3(iphi)
		n=iimp4(iphi)
		kk=3*(k-1)+1
		ll=3*(l-1)+1
		mm=3*(m-1)+1
		nn=3*(n-1)+1
		iiphi=54*(iphi-1)+1
c pair kl
		d2v(kk,ll) = d2v(kk,ll) + d2imp(iiphi)
		d2v(ll,kk) = d2v(ll,kk) + d2imp(iiphi)
		d2v(kk+1,ll) = d2v(kk+1,ll) + d2imp(iiphi+1)
		d2v(ll,kk+1) = d2v(ll,kk+1) + d2imp(iiphi+1)
		d2v(kk+2,ll) = d2v(kk+2,ll) + d2imp(iiphi+2)
		d2v(ll,kk+2) = d2v(ll,kk+2) + d2imp(iiphi+2)
		d2v(kk,ll+1) = d2v(kk,ll+1) + d2imp(iiphi+3)
		d2v(ll+1,kk) = d2v(ll+1,kk) + d2imp(iiphi+3)
		d2v(kk+1,ll+1) = d2v(kk+1,ll+1) + d2imp(iiphi+4)
		d2v(ll+1,kk+1) = d2v(ll+1,kk+1) + d2imp(iiphi+4)
		d2v(kk+2,ll+1) = d2v(kk+2,ll+1) + d2imp(iiphi+5)
		d2v(ll+1,kk+2) = d2v(ll+1,kk+2) + d2imp(iiphi+5)
		d2v(kk,ll+2) = d2v(kk,ll+2) + d2imp(iiphi+6)
		d2v(ll+2,kk) = d2v(ll+2,kk) + d2imp(iiphi+6)
		d2v(kk+1,ll+2) = d2v(kk+1,ll+2) + d2imp(iiphi+7)
		d2v(ll+2,kk+1) = d2v(ll+2,kk+1) + d2imp(iiphi+7)
		d2v(kk+2,ll+2) = d2v(kk+2,ll+2) + d2imp(iiphi+8)
		d2v(ll+2,kk+2) = d2v(ll+2,kk+2) + d2imp(iiphi+8)

c pair km
		d2v(kk,mm) = d2v(kk,mm) + d2imp(iiphi+9)
		d2v(mm,kk) = d2v(mm,kk) + d2imp(iiphi+9)
		d2v(kk+1,mm) = d2v(kk+1,mm) + d2imp(iiphi+10)
		d2v(mm,kk+1) = d2v(mm,kk+1) + d2imp(iiphi+10)
		d2v(kk+2,mm) = d2v(kk+2,mm) + d2imp(iiphi+11)
		d2v(mm,kk+2) = d2v(mm,kk+2) + d2imp(iiphi+11)
		d2v(kk,mm+1) = d2v(kk,mm+1) + d2imp(iiphi+12)
		d2v(mm+1,kk) = d2v(mm+1,kk) + d2imp(iiphi+12)
		d2v(kk+1,mm+1) = d2v(kk+1,mm+1) + d2imp(iiphi+13)
		d2v(mm+1,kk+1) = d2v(mm+1,kk+1) + d2imp(iiphi+13)
		d2v(kk+2,mm+1) = d2v(kk+2,mm+1) + d2imp(iiphi+14)
		d2v(mm+1,kk+2) = d2v(mm+1,kk+2) + d2imp(iiphi+14)
		d2v(kk,mm+2) = d2v(kk,mm+2) + d2imp(iiphi+15)
		d2v(mm+2,kk) = d2v(mm+2,kk) + d2imp(iiphi+15)
		d2v(kk+1,mm+2) = d2v(kk+1,mm+2) + d2imp(iiphi+16)
		d2v(mm+2,kk+1) = d2v(mm+2,kk+1) + d2imp(iiphi+16)
		d2v(kk+2,mm+2) = d2v(kk+2,mm+2) + d2imp(iiphi+17)
		d2v(mm+2,kk+2) = d2v(mm+2,kk+2) + d2imp(iiphi+17)

c pair kn
		d2v(kk,nn) = d2v(kk,nn) + d2imp(iiphi+18)
		d2v(nn,kk) = d2v(nn,kk) + d2imp(iiphi+18)
		d2v(kk+1,nn) = d2v(kk+1,nn) + d2imp(iiphi+19)
		d2v(nn,kk+1) = d2v(nn,kk+1) + d2imp(iiphi+19)
		d2v(kk+2,nn) = d2v(kk+2,nn) + d2imp(iiphi+20)
		d2v(nn,kk+2) = d2v(nn,kk+2) + d2imp(iiphi+20)
		d2v(kk,nn+1) = d2v(kk,nn+1) + d2imp(iiphi+21)
		d2v(nn+1,kk) = d2v(nn+1,kk) + d2imp(iiphi+21)
		d2v(kk+1,nn+1) = d2v(kk+1,nn+1) + d2imp(iiphi+22)
		d2v(nn+1,kk+1) = d2v(nn+1,kk+1) + d2imp(iiphi+22)
		d2v(kk+2,nn+1) = d2v(kk+2,nn+1) + d2imp(iiphi+23)
		d2v(nn+1,kk+2) = d2v(nn+1,kk+2) + d2imp(iiphi+23)
		d2v(kk,nn+2) = d2v(kk,nn+2) + d2imp(iiphi+24)
		d2v(nn+2,kk) = d2v(nn+2,kk) + d2imp(iiphi+24)
		d2v(kk+1,nn+2) = d2v(kk+1,nn+2) + d2imp(iiphi+25)
		d2v(nn+2,kk+1) = d2v(nn+2,kk+1) + d2imp(iiphi+25)
		d2v(kk+2,nn+2) = d2v(kk+2,nn+2) + d2imp(iiphi+26)
		d2v(nn+2,kk+2) = d2v(nn+2,kk+2) + d2imp(iiphi+26)

c pair lm
		d2v(ll,mm) = d2v(ll,mm) + d2imp(iiphi+27)
		d2v(mm,ll) = d2v(mm,ll) + d2imp(iiphi+27)
		d2v(ll+1,mm) = d2v(ll+1,mm) + d2imp(iiphi+28)
		d2v(mm,ll+1) = d2v(mm,ll+1) + d2imp(iiphi+28)
		d2v(ll+2,mm) = d2v(ll+2,mm) + d2imp(iiphi+29)
		d2v(mm,ll+2) = d2v(mm,ll+2) + d2imp(iiphi+29)
		d2v(ll,mm+1) = d2v(ll,mm+1) + d2imp(iiphi+30)
		d2v(mm+1,ll) = d2v(mm+1,ll) + d2imp(iiphi+30)
		d2v(ll+1,mm+1) = d2v(ll+1,mm+1) + d2imp(iiphi+31)
		d2v(mm+1,ll+1) = d2v(mm+1,ll+1) + d2imp(iiphi+31)
		d2v(ll+2,mm+1) = d2v(ll+2,mm+1) + d2imp(iiphi+32)
		d2v(mm+1,ll+2) = d2v(mm+1,ll+2) + d2imp(iiphi+32)
		d2v(ll,mm+2) = d2v(ll,mm+2) + d2imp(iiphi+33)
		d2v(mm+2,ll) = d2v(mm+2,ll) + d2imp(iiphi+33)
		d2v(ll+1,mm+2) = d2v(ll+1,mm+2) + d2imp(iiphi+34)
		d2v(mm+2,ll+1) = d2v(mm+2,ll+1) + d2imp(iiphi+34)
		d2v(ll+2,mm+2) = d2v(ll+2,mm+2) + d2imp(iiphi+35)
		d2v(mm+2,ll+2) = d2v(mm+2,ll+2) + d2imp(iiphi+35)

c pair ln
		d2v(ll,nn) = d2v(ll,nn) + d2imp(iiphi+36)
		d2v(nn,ll) = d2v(nn,ll) + d2imp(iiphi+36)
		d2v(ll+1,nn) = d2v(ll+1,nn) + d2imp(iiphi+37)
		d2v(nn,ll+1) = d2v(nn,ll+1) + d2imp(iiphi+37)
		d2v(ll+2,nn) = d2v(ll+2,nn) + d2imp(iiphi+38)
		d2v(nn,ll+2) = d2v(nn,ll+2) + d2imp(iiphi+38)
		d2v(ll,nn+1) = d2v(ll,nn+1) + d2imp(iiphi+39)
		d2v(nn+1,ll) = d2v(nn+1,ll) + d2imp(iiphi+39)
		d2v(ll+1,nn+1) = d2v(ll+1,nn+1) + d2imp(iiphi+40)
		d2v(nn+1,ll+1) = d2v(nn+1,ll+1) + d2imp(iiphi+40)
		d2v(ll+2,nn+1) = d2v(ll+2,nn+1) + d2imp(iiphi+41)
		d2v(nn+1,ll+2) = d2v(nn+1,ll+2) + d2imp(iiphi+41)
		d2v(ll,nn+2) = d2v(ll,nn+2) + d2imp(iiphi+42)
		d2v(nn+2,ll) = d2v(nn+2,ll) + d2imp(iiphi+42)
		d2v(ll+1,nn+2) = d2v(ll+1,nn+2) + d2imp(iiphi+43)
		d2v(nn+2,ll+1) = d2v(nn+2,ll+1) + d2imp(iiphi+43)
		d2v(ll+2,nn+2) = d2v(ll+2,nn+2) + d2imp(iiphi+44)
		d2v(nn+2,ll+2) = d2v(nn+2,ll+2) + d2imp(iiphi+44)

c pair mn
		d2v(mm,nn) = d2v(mm,nn) + d2imp(iiphi+45)
		d2v(nn,mm) = d2v(nn,mm) + d2imp(iiphi+45)
		d2v(mm+1,nn) = d2v(mm+1,nn) + d2imp(iiphi+46)
		d2v(nn,mm+1) = d2v(nn,mm+1) + d2imp(iiphi+46)
		d2v(mm+2,nn) = d2v(mm+2,nn) + d2imp(iiphi+47)
		d2v(nn,mm+2) = d2v(nn,mm+2) + d2imp(iiphi+47)
		d2v(mm,nn+1) = d2v(mm,nn+1) + d2imp(iiphi+48)
		d2v(nn+1,mm) = d2v(nn+1,mm) + d2imp(iiphi+48)
		d2v(mm+1,nn+1) = d2v(mm+1,nn+1) + d2imp(iiphi+49)
		d2v(nn+1,mm+1) = d2v(nn+1,mm+1) + d2imp(iiphi+49)
		d2v(mm+2,nn+1) = d2v(mm+2,nn+1) + d2imp(iiphi+50)
		d2v(nn+1,mm+2) = d2v(nn+1,mm+2) + d2imp(iiphi+50)
		d2v(mm,nn+2) = d2v(mm,nn+2) + d2imp(iiphi+51)
		d2v(nn+2,mm) = d2v(nn+2,mm) + d2imp(iiphi+51)
		d2v(mm+1,nn+2) = d2v(mm+1,nn+2) + d2imp(iiphi+52)
		d2v(nn+2,mm+1) = d2v(nn+2,mm+1) + d2imp(iiphi+52)
		d2v(mm+2,nn+2) = d2v(mm+2,nn+2) + d2imp(iiphi+53)
		d2v(nn+2,mm+2) = d2v(nn+2,mm+2) + d2imp(iiphi+53)
45	continue
	end if

	write(unit,*)' Second derivative matrix '
	do 29 i=1,3*npt
	 write(unit,100)(d2v(i,j),j=1,3*maxpt2d)
100	 format(12(f5.1,1x))
29	continue

	return 
	end
