	subroutine matvec(vector)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c On input the vector "vector" is provided by the user and on
c output it is modified to the product of the matrix by the
c vector.
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
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

	double precision vector(*)
	double precision vectmp(3*maxpt)
	integer i,j,k,m,l,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
	integer jbeg4,jend4,iphi,iiphi,ith,iith
	integer ii,jj,kk,ll,mm,nn

c	<<*********************** added by jing ***********************
c       For using two data files act as two big array d2nb1()
c       and d2nb_2() respectively. All stuff added by jing are
c       included in two "****************" comment lines, and old
c       statement will add a "c*" in each of them.

        integer fid1, fid2

        parameter (fid1 = 32)
        parameter (fid2 = 33)

        double precision d1, d2, d3, d4, d5, d6

c       open two data files use fid1, fid2 as units respectively
        open (unit=fid1, file='d2nb_1.data', access='sequential',
     $        form='unformatted', status='old')
        open (unit=fid2, file='d2nb_2.data', access='sequential',
     $        form='unformatted', status='old')
	rewind fid1
	rewind fid2

c       *************************************************************>>

	do 2 i=1,npt
c Contribution of diagonal terms
c Note that the diag matrix is organized as:
c i=xx i+1=xy=yx i+2=xz=zx i+3=yy i+4=yz=zy i+5=zz
c
	 ii = 3*(i-1)+1
	 jj = 6*(i-1)+1
	 vectmp(ii) = diag(jj)*vector(ii) + diag(jj+1)*vector(ii+1)
     1		+ diag(jj+2)*vector(ii+2) 
	 vectmp(ii+1) = diag(jj+1)*vector(ii) + diag(jj+3)*vector(ii+1)
     1		+ diag(jj+4)*vector(ii+2)
	 vectmp(ii+2) = diag(jj+2)*vector(ii) + diag(jj+4)*vector(ii+1)
     1		+ diag(jj+5)*vector(ii+2)
2	continue
c Non-bonded contribution
c Start loop on short range van der Waals and elec. interactions.
c
	jbeg1 = 1
	jbeg2 = 1
	jbeg3 = 1
	jbeg4 = 1

	if (evdyes.or.eelyes) then
	do 7 i=1,npt-1
	 jend1 = point1(i)
	 do 3 k=jbeg1,jend1
c Off diagonal elements: for each pair i and j, the storage is
c kk=xx kk+1=xy=yx kk+2=xz=zx kk+3=yy kk+4=yz=zy kk+5=zz
c
	  j=list1(k)
	  kk = 6*(k-1)+1
	  ii = 3*(i-1)+1
	  jj = 3*(j-1)+1

c*	  vectmp(ii) = vectmp(ii) + d2nb_1(kk)*vector(jj) +
c*     1	  	d2nb_1(kk+1)*vector(jj+1) + d2nb_1(kk+2)*vector(jj+2)
c*	  vectmp(ii+1) = vectmp(ii+1) + d2nb_1(kk+1)*vector(jj) +
c*     1	  	d2nb_1(kk+3)*vector(jj+1) + d2nb_1(kk+4)*vector(jj+2)
c*	  vectmp(ii+2) = vectmp(ii+2) + d2nb_1(kk+2)*vector(jj) +
c*     1	  	d2nb_1(kk+4)*vector(jj+1) + d2nb_1(kk+5)*vector(jj+2)
c*	  vectmp(jj) = vectmp(jj) + d2nb_1(kk)*vector(ii) +
c*     1	  	d2nb_1(kk+1)*vector(ii+1) + d2nb_1(kk+2)*vector(ii+2)
c*	  vectmp(jj+1) = vectmp(jj+1) + d2nb_1(kk+1)*vector(ii) +
c*     1	  	d2nb_1(kk+3)*vector(ii+1) + d2nb_1(kk+4)*vector(ii+2)
c*	  vectmp(jj+2) = vectmp(jj+2) + d2nb_1(kk+2)*vector(ii) +
c*     1	  	d2nb_1(kk+4)*vector(ii+1) + d2nb_1(kk+5)*vector(ii+2)

c	<<*************************************************************

c	  read d1-d6 from data file d2nb_1.data
	  read (fid1, end=51) d1, d2, d3, d4, d5, d6

          vectmp(ii) = vectmp(ii) + d1*vector(jj) +
     1          d2*vector(jj+1) + d3*vector(jj+2)
          vectmp(ii+1) = vectmp(ii+1) + d2*vector(jj) +
     1          d4*vector(jj+1) + d5*vector(jj+2)
          vectmp(ii+2) = vectmp(ii+2) + d3*vector(jj) +
     1          d5*vector(jj+1) + d6*vector(jj+2)
          vectmp(jj) = vectmp(jj) + d1*vector(ii) +
     1          d2*vector(ii+1) + d3*vector(ii+2)
          vectmp(jj+1) = vectmp(jj+1) + d2*vector(ii) +
     1          d4*vector(ii+1) + d5*vector(ii+2)
          vectmp(jj+2) = vectmp(jj+2) + d3*vector(ii) +
     1          d5*vector(ii+1) + d6*vector(ii+2)

c	*************************************************************>>

3	 continue

	 jbeg1 = jend1 + 1

	 jend2 = point2(i)
	 do 4 k=jbeg2,jend2
	  j=list2(k)
	  kk = 6*(k-1)+1
	  ii = 3*(i-1)+1
	  jj = 3*(j-1)+1

c*	  vectmp(ii) = vectmp(ii) + d2nb_2(kk)*vector(jj) +
c*     1	  	d2nb_2(kk+1)*vector(jj+1) + d2nb_2(kk+2)*vector(jj+2)
c*	  vectmp(ii+1) = vectmp(ii+1) + d2nb_2(kk+1)*vector(jj) +
c*     1	  	d2nb_2(kk+3)*vector(jj+1) + d2nb_2(kk+4)*vector(jj+2)
c*	  vectmp(ii+2) = vectmp(ii+2) + d2nb_2(kk+2)*vector(jj) +
c*     1	  	d2nb_2(kk+4)*vector(jj+1) + d2nb_2(kk+5)*vector(jj+2)
c*	  vectmp(jj) = vectmp(jj) + d2nb_2(kk)*vector(ii) +
c*     1	  	d2nb_2(kk+1)*vector(ii+1) + d2nb_2(kk+2)*vector(ii+2)
c*	  vectmp(jj+1) = vectmp(jj+1) + d2nb_2(kk+1)*vector(ii) +
c*     1	  	d2nb_2(kk+3)*vector(ii+1) + d2nb_2(kk+4)*vector(ii+2)
c* 	  vectmp(jj+2) = vectmp(jj+2) + d2nb_2(kk+2)*vector(ii) +
c*     1	  	d2nb_2(kk+4)*vector(ii+1) + d2nb_2(kk+5)*vector(ii+2)

c	<<*******************************************************************

c	  read d1-d6 from data file d2_nb_2.data
	  read (fid2, end=52) d1, d2, d3, d4, d5, d6

	  vectmp(ii) = vectmp(ii) + d1*vector(jj) +
     1          d2*vector(jj+1) + d3*vector(jj+2)
          vectmp(ii+1) = vectmp(ii+1) + d2*vector(jj) +
     1          d4*vector(jj+1) + d5*vector(jj+2)
          vectmp(ii+2) = vectmp(ii+2) + d3*vector(jj) +
     1          d5*vector(jj+1) + d6*vector(jj+2)
          vectmp(jj) = vectmp(jj) + d1*vector(ii) +
     1          d2*vector(ii+1) + d3*vector(ii+2)
          vectmp(jj+1) = vectmp(jj+1) + d2*vector(ii) +
     1          d4*vector(ii+1) + d5*vector(ii+2)
          vectmp(jj+2) = vectmp(jj+2) + d3*vector(ii) +
     1          d5*vector(ii+1) + d6*vector(ii+2)

c	*******************************************************************>>

4	 continue

	 jbeg2 = jend2 + 1

	 jend3 = point3(i)
	 do 5 k=jbeg3,jend3
	  j=list3(k)
	  kk = 6*(k-1)+1
	  ii = 3*(i-1)+1
	  jj = 3*(j-1)+1
	  vectmp(ii) = vectmp(ii) + d2nb_3(kk)*vector(jj) +
     1	  	d2nb_3(kk+1)*vector(jj+1) + d2nb_3(kk+2)*vector(jj+2)
	  vectmp(ii+1) = vectmp(ii+1) + d2nb_3(kk+1)*vector(jj) +
     1	  	d2nb_3(kk+3)*vector(jj+1) + d2nb_3(kk+4)*vector(jj+2)
	  vectmp(ii+2) = vectmp(ii+2) + d2nb_3(kk+2)*vector(jj) +
     1	  	d2nb_3(kk+4)*vector(jj+1) + d2nb_3(kk+5)*vector(jj+2)
	  vectmp(jj) = vectmp(jj) + d2nb_3(kk)*vector(ii) +
     1	  	d2nb_3(kk+1)*vector(ii+1) + d2nb_3(kk+2)*vector(ii+2)
	  vectmp(jj+1) = vectmp(jj+1) + d2nb_3(kk+1)*vector(ii) +
     1	  	d2nb_3(kk+3)*vector(ii+1) + d2nb_3(kk+4)*vector(ii+2)
	  vectmp(jj+2) = vectmp(jj+2) + d2nb_3(kk+2)*vector(ii) +
     1	  	d2nb_3(kk+4)*vector(ii+1) + d2nb_3(kk+5)*vector(ii+2)
5	 continue
	 jbeg3 = jend3 + 1

c@ lable 111 added for debuging
c@111	continue
	 

	 jend4 = spec1(i)
	 do 6 k=jbeg4,jend4
	  j=spec2(k)
	  kk = 6*(k-1)+1
	  ii = 3*(i-1)+1
	  jj = 3*(j-1)+1
	  vectmp(ii) = vectmp(ii) + d2spec(kk)*vector(jj) +
     1	  	d2spec(kk+1)*vector(jj+1) + d2spec(kk+2)*vector(jj+2)
	  vectmp(ii+1) = vectmp(ii+1) + d2spec(kk+1)*vector(jj) +
     1	  	d2spec(kk+3)*vector(jj+1) + d2spec(kk+4)*vector(jj+2)
	  vectmp(ii+2) = vectmp(ii+2) + d2spec(kk+2)*vector(jj) +
     1	  	d2spec(kk+4)*vector(jj+1) + d2spec(kk+5)*vector(jj+2)
	  vectmp(jj) = vectmp(jj) + d2spec(kk)*vector(ii) +
     1	  	d2spec(kk+1)*vector(ii+1) + d2spec(kk+2)*vector(ii+2)
	  vectmp(jj+1) = vectmp(jj+1) + d2spec(kk+1)*vector(ii) +
     1	  	d2spec(kk+3)*vector(ii+1) + d2spec(kk+4)*vector(ii+2)
	  vectmp(jj+2) = vectmp(jj+2) + d2spec(kk+2)*vector(ii) +
     1	  	d2spec(kk+4)*vector(ii+1) + d2spec(kk+5)*vector(ii+2)
6	 continue
	 jbeg4 = jend4 + 1

7	continue
	end if

c Bond contribution

	if (ebyes) then
	do 8 k=1,nb
	  kk = 6*(k-1)+1
	  ii = 3*(ib1(k)-1)+1
	  jj = 3*(ib2(k)-1)+1
	  vectmp(ii) = vectmp(ii) + d2bond(kk)*vector(jj) +
     1	  	d2bond(kk+1)*vector(jj+1) + d2bond(kk+2)*vector(jj+2)
	  vectmp(ii+1) = vectmp(ii+1) + d2bond(kk+1)*vector(jj) +
     1	  	d2bond(kk+3)*vector(jj+1) + d2bond(kk+4)*vector(jj+2)
	  vectmp(ii+2) = vectmp(ii+2) + d2bond(kk+2)*vector(jj) +
     1	  	d2bond(kk+4)*vector(jj+1) + d2bond(kk+5)*vector(jj+2)
	  vectmp(jj) = vectmp(jj) + d2bond(kk)*vector(ii) +
     1	  	d2bond(kk+1)*vector(ii+1) + d2bond(kk+2)*vector(ii+2)
	  vectmp(jj+1) = vectmp(jj+1) + d2bond(kk+1)*vector(ii) +
     1	  	d2bond(kk+3)*vector(ii+1) + d2bond(kk+4)*vector(ii+2)
	  vectmp(jj+2) = vectmp(jj+2) + d2bond(kk+2)*vector(ii) +
     1	  	d2bond(kk+4)*vector(ii+1) + d2bond(kk+5)*vector(ii+2)
8	continue
	end if

c Angle contribution

	if (ethyes) then
	do 9 ith=1,nangl
	 iith= 27*(ith-1)+1
	 ii  = 3*(iangl1(ith)-1)+1
	 jj  = 3*(iangl2(ith)-1)+1
	 kk  = 3*(iangl3(ith)-1)+1
c pair i,j
	 vectmp(ii) = vectmp(ii) + d2theta(iith)*vector(jj) +
     1	  d2theta(iith+1)*vector(jj+1) + d2theta(iith+2)*vector(jj+2)
	 vectmp(ii+1) = vectmp(ii+1) + d2theta(iith+3)*vector(jj) +
     1	  d2theta(iith+4)*vector(jj+1) + d2theta(iith+5)*vector(jj+2)
	 vectmp(ii+2) = vectmp(ii+2) + d2theta(iith+6)*vector(jj) +
     1	  d2theta(iith+7)*vector(jj+1) + d2theta(iith+8)*vector(jj+2)
	 vectmp(jj) = vectmp(jj) + d2theta(iith)*vector(ii) +
     1	  d2theta(iith+3)*vector(ii+1) + d2theta(iith+6)*vector(ii+2)
	 vectmp(jj+1) = vectmp(jj+1) + d2theta(iith+1)*vector(ii) +
     1	  d2theta(iith+4)*vector(ii+1) + d2theta(iith+7)*vector(ii+2)
	 vectmp(jj+2) = vectmp(jj+2) + d2theta(iith+2)*vector(ii) +
     1	  d2theta(iith+5)*vector(ii+1) + d2theta(iith+8)*vector(ii+2)
c pair i,k
	 vectmp(ii) = vectmp(ii) + d2theta(iith+9)*vector(kk) +
     1	  d2theta(iith+10)*vector(kk+1) + d2theta(iith+11)*vector(kk+2)
	 vectmp(ii+1) = vectmp(ii+1) + d2theta(iith+12)*vector(kk) +
     1	  d2theta(iith+13)*vector(kk+1) + d2theta(iith+14)*vector(kk+2)
	 vectmp(ii+2) = vectmp(ii+2) + d2theta(iith+15)*vector(kk) +
     1	  d2theta(iith+16)*vector(kk+1) + d2theta(iith+17)*vector(kk+2)
	 vectmp(kk) = vectmp(kk) + d2theta(iith+9)*vector(ii) +
     1	  d2theta(iith+12)*vector(ii+1) + d2theta(iith+15)*vector(ii+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2theta(iith+10)*vector(ii) +
     1	  d2theta(iith+13)*vector(ii+1) + d2theta(iith+16)*vector(ii+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2theta(iith+11)*vector(ii) +
     1	  d2theta(iith+14)*vector(ii+1) + d2theta(iith+17)*vector(ii+2)
c pair j,k
	 vectmp(jj) = vectmp(jj) + d2theta(iith+18)*vector(kk) +
     1	  d2theta(iith+19)*vector(kk+1) + d2theta(iith+20)*vector(kk+2)
	 vectmp(jj+1) = vectmp(jj+1) + d2theta(iith+21)*vector(kk) +
     1	  d2theta(iith+22)*vector(kk+1) + d2theta(iith+23)*vector(kk+2)
	 vectmp(jj+2) = vectmp(jj+2) + d2theta(iith+24)*vector(kk) +
     1	  d2theta(iith+25)*vector(kk+1) + d2theta(iith+26)*vector(kk+2)
	 vectmp(kk) = vectmp(kk) + d2theta(iith+18)*vector(jj) +
     1	  d2theta(iith+21)*vector(jj+1) + d2theta(iith+24)*vector(jj+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2theta(iith+19)*vector(jj) +
     1	  d2theta(iith+22)*vector(jj+1) + d2theta(iith+25)*vector(jj+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2theta(iith+20)*vector(jj) +
     1	  d2theta(iith+23)*vector(jj+1) + d2theta(iith+26)*vector(jj+2)
9	continue
	end if

c contribution of torsions

	if (etoyes) then
	do 10 iphi=1,ntors
	 iiphi = 54*(iphi-1)+1
	 kk    = 3*(itor1(iphi)-1)+1
	 ll    = 3*(itor2(iphi)-1)+1
	 mm    = 3*(itor3(iphi)-1)+1
	 nn    = 3*(itor4(iphi)-1)+1
c k,l pair
	 vectmp(kk)   = vectmp(kk)   + d2phi(iiphi)*vector(ll) + 
     1    d2phi(iiphi+3)*vector(ll+1) + d2phi(iiphi+6)*vector(ll+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2phi(iiphi+1)*vector(ll) +
     1    d2phi(iiphi+4)*vector(ll+1) + d2phi(iiphi+7)*vector(ll+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2phi(iiphi+2)*vector(ll) +
     1	  d2phi(iiphi+5)*vector(ll+1) + d2phi(iiphi+8)*vector(ll+2)
	 vectmp(ll)   = vectmp(ll)   + d2phi(iiphi)*vector(kk) +
     1	  d2phi(iiphi+1)*vector(kk+1) + d2phi(iiphi+2)*vector(kk+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2phi(iiphi+3)*vector(kk) +
     1 	  d2phi(iiphi+4)*vector(kk+1) + d2phi(iiphi+5)*vector(kk+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2phi(iiphi+6)*vector(kk) +
     1	  d2phi(iiphi+7)*vector(kk+1) + d2phi(iiphi+8)*vector(kk+2)

c k,m
	 vectmp(kk)   = vectmp(kk)   + d2phi(iiphi+9)*vector(mm) +
     1	  d2phi(iiphi+12)*vector(mm+1) + d2phi(iiphi+15)*vector(mm+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2phi(iiphi+10)*vector(mm) +
     1	  d2phi(iiphi+13)*vector(mm+1) + d2phi(iiphi+16)*vector(mm+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2phi(iiphi+11)*vector(mm) +
     1	  d2phi(iiphi+14)*vector(mm+1) + d2phi(iiphi+17)*vector(mm+2)
	 vectmp(mm)   = vectmp(mm)   + d2phi(iiphi+9)*vector(kk) +
     1	  d2phi(iiphi+10)*vector(kk+1) + d2phi(iiphi+11)*vector(kk+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2phi(iiphi+12)*vector(kk) +
     1	  d2phi(iiphi+13)*vector(kk+1) + d2phi(iiphi+14)*vector(kk+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2phi(iiphi+15)*vector(kk) +
     1	  d2phi(iiphi+16)*vector(kk+1) + d2phi(iiphi+17)*vector(kk+2)

c k,n
	 vectmp(kk)   = vectmp(kk)   + d2phi(iiphi+18)*vector(nn) +
     1	  d2phi(iiphi+21)*vector(nn+1) + d2phi(iiphi+24)*vector(nn+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2phi(iiphi+19)*vector(nn) +
     1	  d2phi(iiphi+22)*vector(nn+1) + d2phi(iiphi+25)*vector(nn+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2phi(iiphi+20)*vector(nn) +
     1	  d2phi(iiphi+23)*vector(nn+1) + d2phi(iiphi+26)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2phi(iiphi+18)*vector(kk) +
     1	  d2phi(iiphi+19)*vector(kk+1) + d2phi(iiphi+20)*vector(kk+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2phi(iiphi+21)*vector(kk) +
     1	  d2phi(iiphi+22)*vector(kk+1) + d2phi(iiphi+23)*vector(kk+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2phi(iiphi+24)*vector(kk) +
     1	  d2phi(iiphi+25)*vector(kk+1) + d2phi(iiphi+26)*vector(kk+2)

c l,m
	 vectmp(ll)   = vectmp(ll)   + d2phi(iiphi+27)*vector(mm) +
     1	  d2phi(iiphi+30)*vector(mm+1) + d2phi(iiphi+33)*vector(mm+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2phi(iiphi+28)*vector(mm) +
     1	  d2phi(iiphi+31)*vector(mm+1) + d2phi(iiphi+34)*vector(mm+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2phi(iiphi+29)*vector(mm) +
     1	  d2phi(iiphi+32)*vector(mm+1) + d2phi(iiphi+35)*vector(mm+2)
	 vectmp(mm)   = vectmp(mm)   + d2phi(iiphi+27)*vector(ll) +
     1	  d2phi(iiphi+28)*vector(ll+1) + d2phi(iiphi+29)*vector(ll+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2phi(iiphi+30)*vector(ll) +
     1	  d2phi(iiphi+31)*vector(ll+1) + d2phi(iiphi+32)*vector(ll+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2phi(iiphi+33)*vector(ll) +
     1	  d2phi(iiphi+34)*vector(ll+1) + d2phi(iiphi+35)*vector(ll+2)

c l,n
	 vectmp(ll)   = vectmp(ll)   + d2phi(iiphi+36)*vector(nn) +
     1	  d2phi(iiphi+39)*vector(nn+1) + d2phi(iiphi+42)*vector(nn+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2phi(iiphi+37)*vector(nn) +
     1	  d2phi(iiphi+40)*vector(nn+1) + d2phi(iiphi+43)*vector(nn+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2phi(iiphi+38)*vector(nn) +
     1	  d2phi(iiphi+41)*vector(nn+1) + d2phi(iiphi+44)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2phi(iiphi+36)*vector(ll) +
     1	  d2phi(iiphi+37)*vector(ll+1) + d2phi(iiphi+38)*vector(ll+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2phi(iiphi+39)*vector(ll) +
     1	  d2phi(iiphi+40)*vector(ll+1) + d2phi(iiphi+41)*vector(ll+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2phi(iiphi+42)*vector(ll) +
     1	  d2phi(iiphi+43)*vector(ll+1) + d2phi(iiphi+44)*vector(ll+2)

c m,n
	 vectmp(mm)   = vectmp(mm)   + d2phi(iiphi+45)*vector(nn) +
     1	  d2phi(iiphi+48)*vector(nn+1) + d2phi(iiphi+51)*vector(nn+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2phi(iiphi+46)*vector(nn) +
     1	  d2phi(iiphi+49)*vector(nn+1) + d2phi(iiphi+52)*vector(nn+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2phi(iiphi+47)*vector(nn) +
     1	  d2phi(iiphi+50)*vector(nn+1) + d2phi(iiphi+53)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2phi(iiphi+45)*vector(mm) +
     1	  d2phi(iiphi+46)*vector(mm+1) + d2phi(iiphi+47)*vector(mm+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2phi(iiphi+48)*vector(mm) +
     1	  d2phi(iiphi+49)*vector(mm+1) + d2phi(iiphi+50)*vector(mm+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2phi(iiphi+51)*vector(mm) +
     1	  d2phi(iiphi+52)*vector(mm+1) + d2phi(iiphi+53)*vector(mm+2)

10	continue
	end if

c contribution of improper torsions

	if (eimyes) then
	do 11 iphi=1,nimp

	 kk = 3*(iimp1(iphi)-1)+1
	 ll = 3*(iimp2(iphi)-1)+1
	 mm = 3*(iimp3(iphi)-1)+1
	 nn = 3*(iimp4(iphi)-1)+1
	 iiphi = 54*(iphi-1)+1

c k,l pair
	 vectmp(kk)   = vectmp(kk)   + d2imp(iiphi)*vector(ll) + 
     1    d2imp(iiphi+3)*vector(ll+1) + d2imp(iiphi+6)*vector(ll+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2imp(iiphi+1)*vector(ll) +
     1    d2imp(iiphi+4)*vector(ll+1) + d2imp(iiphi+7)*vector(ll+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2imp(iiphi+2)*vector(ll) +
     1	  d2imp(iiphi+5)*vector(ll+1) + d2imp(iiphi+8)*vector(ll+2)
	 vectmp(ll)   = vectmp(ll)   + d2imp(iiphi)*vector(kk) +
     1	  d2imp(iiphi+1)*vector(kk+1) + d2imp(iiphi+2)*vector(kk+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2imp(iiphi+3)*vector(kk) +
     1 	  d2imp(iiphi+4)*vector(kk+1) + d2imp(iiphi+5)*vector(kk+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2imp(iiphi+6)*vector(kk) +
     1	  d2imp(iiphi+7)*vector(kk+1) + d2imp(iiphi+8)*vector(kk+2)

c k,m
	 vectmp(kk)   = vectmp(kk)   + d2imp(iiphi+9)*vector(mm) +
     1	  d2imp(iiphi+12)*vector(mm+1) + d2imp(iiphi+15)*vector(mm+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2imp(iiphi+10)*vector(mm) +
     1	  d2imp(iiphi+13)*vector(mm+1) + d2imp(iiphi+16)*vector(mm+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2imp(iiphi+11)*vector(mm) +
     1	  d2imp(iiphi+14)*vector(mm+1) + d2imp(iiphi+17)*vector(mm+2)
	 vectmp(mm)   = vectmp(mm)   + d2imp(iiphi+9)*vector(kk) +
     1	  d2imp(iiphi+10)*vector(kk+1) + d2imp(iiphi+11)*vector(kk+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2imp(iiphi+12)*vector(kk) +
     1	  d2imp(iiphi+13)*vector(kk+1) + d2imp(iiphi+14)*vector(kk+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2imp(iiphi+15)*vector(kk) +
     1	  d2imp(iiphi+16)*vector(kk+1) + d2imp(iiphi+17)*vector(kk+2)

c k,n
	 vectmp(kk)   = vectmp(kk)   + d2imp(iiphi+18)*vector(nn) +
     1	  d2imp(iiphi+21)*vector(nn+1) + d2imp(iiphi+24)*vector(nn+2)
	 vectmp(kk+1) = vectmp(kk+1) + d2imp(iiphi+19)*vector(nn) +
     1	  d2imp(iiphi+22)*vector(nn+1) + d2imp(iiphi+25)*vector(nn+2)
	 vectmp(kk+2) = vectmp(kk+2) + d2imp(iiphi+20)*vector(nn) +
     1	  d2imp(iiphi+23)*vector(nn+1) + d2imp(iiphi+26)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2imp(iiphi+18)*vector(kk) +
     1	  d2imp(iiphi+19)*vector(kk+1) + d2imp(iiphi+20)*vector(kk+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2imp(iiphi+21)*vector(kk) +
     1	  d2imp(iiphi+22)*vector(kk+1) + d2imp(iiphi+23)*vector(kk+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2imp(iiphi+24)*vector(kk) +
     1	  d2imp(iiphi+25)*vector(kk+1) + d2imp(iiphi+26)*vector(kk+2)

c l,m
	 vectmp(ll)   = vectmp(ll)   + d2imp(iiphi+27)*vector(mm) +
     1	  d2imp(iiphi+30)*vector(mm+1) + d2imp(iiphi+33)*vector(mm+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2imp(iiphi+28)*vector(mm) +
     1	  d2imp(iiphi+31)*vector(mm+1) + d2imp(iiphi+34)*vector(mm+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2imp(iiphi+29)*vector(mm) +
     1	  d2imp(iiphi+32)*vector(mm+1) + d2imp(iiphi+35)*vector(mm+2)
	 vectmp(mm)   = vectmp(mm)   + d2imp(iiphi+27)*vector(ll) +
     1	  d2imp(iiphi+28)*vector(ll+1) + d2imp(iiphi+29)*vector(ll+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2imp(iiphi+30)*vector(ll) +
     1	  d2imp(iiphi+31)*vector(ll+1) + d2imp(iiphi+32)*vector(ll+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2imp(iiphi+33)*vector(ll) +
     1	  d2imp(iiphi+34)*vector(ll+1) + d2imp(iiphi+35)*vector(ll+2)

c l,n
	 vectmp(ll)   = vectmp(ll)   + d2imp(iiphi+36)*vector(nn) +
     1	  d2imp(iiphi+39)*vector(nn+1) + d2imp(iiphi+42)*vector(nn+2)
	 vectmp(ll+1) = vectmp(ll+1) + d2imp(iiphi+37)*vector(nn) +
     1	  d2imp(iiphi+40)*vector(nn+1) + d2imp(iiphi+43)*vector(nn+2)
	 vectmp(ll+2) = vectmp(ll+2) + d2imp(iiphi+38)*vector(nn) +
     1	  d2imp(iiphi+41)*vector(nn+1) + d2imp(iiphi+44)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2imp(iiphi+36)*vector(ll) +
     1	  d2imp(iiphi+37)*vector(ll+1) + d2imp(iiphi+38)*vector(ll+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2imp(iiphi+39)*vector(ll) +
     1	  d2imp(iiphi+40)*vector(ll+1) + d2imp(iiphi+41)*vector(ll+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2imp(iiphi+42)*vector(ll) +
     1	  d2imp(iiphi+43)*vector(ll+1) + d2imp(iiphi+44)*vector(ll+2)

c m,n
	 vectmp(mm)   = vectmp(mm)   + d2imp(iiphi+45)*vector(nn) +
     1	  d2imp(iiphi+48)*vector(nn+1) + d2imp(iiphi+51)*vector(nn+2)
	 vectmp(mm+1) = vectmp(mm+1) + d2imp(iiphi+46)*vector(nn) +
     1	  d2imp(iiphi+49)*vector(nn+1) + d2imp(iiphi+52)*vector(nn+2)
	 vectmp(mm+2) = vectmp(mm+2) + d2imp(iiphi+47)*vector(nn) +
     1	  d2imp(iiphi+50)*vector(nn+1) + d2imp(iiphi+53)*vector(nn+2)
	 vectmp(nn)   = vectmp(nn)   + d2imp(iiphi+45)*vector(mm) +
     1	  d2imp(iiphi+46)*vector(mm+1) + d2imp(iiphi+47)*vector(mm+2)
	 vectmp(nn+1) = vectmp(nn+1) + d2imp(iiphi+48)*vector(mm) +
     1	  d2imp(iiphi+49)*vector(mm+1) + d2imp(iiphi+50)*vector(mm+2)
	 vectmp(nn+2) = vectmp(nn+2) + d2imp(iiphi+51)*vector(mm) +
     1	  d2imp(iiphi+52)*vector(mm+1) + d2imp(iiphi+53)*vector(mm+2)

11	continue
	end if


c---This is added by M. Kara-Ivanov for hydrophobic Interactions
c---i and j are numbers of beta atoms in general list
c---k -number of interaction
c---kk,jj number in the matrix



	k=0
	  do 20 l=1,nbeta-1
	     i=betap(l)
	     ii = 3*(i-1)+1
	     do 30 m=l+1,nbeta
		j=betap(m)
		jj = 3*(j-1)+1
		k=k+1
		kk = 6*(k-1)+1


	  vectmp(ii) = vectmp(ii) + d2hyd(kk)*vector(jj) +
     1	  	d2hyd(kk+1)*vector(jj+1) + d2hyd(kk+2)*vector(jj+2)
	  vectmp(ii+1) = vectmp(ii+1) + d2hyd(kk+1)*vector(jj) +
     1	  	d2hyd(kk+3)*vector(jj+1) + d2hyd(kk+4)*vector(jj+2)
	  vectmp(ii+2) = vectmp(ii+2) + d2hyd(kk+2)*vector(jj) +
     1	  	d2hyd(kk+4)*vector(jj+1) + d2hyd(kk+5)*vector(jj+2)
	  vectmp(jj) = vectmp(jj) + d2hyd(kk)*vector(ii) +
     1	  	d2hyd(kk+1)*vector(ii+1) + d2hyd(kk+2)*vector(ii+2)
	  vectmp(jj+1) = vectmp(jj+1) + d2hyd(kk+1)*vector(ii) +
     1	  	d2hyd(kk+3)*vector(ii+1) + d2hyd(kk+4)*vector(ii+2)
	  vectmp(jj+2) = vectmp(jj+2) + d2hyd(kk+2)*vector(ii) +
     1	  	d2hyd(kk+4)*vector(ii+1) + d2hyd(kk+5)*vector(ii+2)
 30	continue
 20	continue

c@ lable 111 added for debuging

c@ lable 111 added for debuging
c@111	continue
	do 12 i=1,3*npt
	 vector(i) = vectmp(i)
12	continue

c       <<**********************************************************

c       close temp data file d2nb_1.data and d2nb_2.data
        close (unit=fid1, status='delete')
        close (unit=fid2, status='delete')

c       **********************************************************>>

	return
c*	end

c	<<**********************************************************

51	print *, 'end of datafile d2nb_1.data'
	stop
52	print *, 'end of data file d2nb_2.data'
	end

c	**********************************************************>>
