	subroutine massweig()

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

	double precision sqtmp
	integer i,j,k,l,m,n
	integer idiag,id,ibond,ib,ith,itheta
	integer itors,iphi,jbeg,jend,kndex

c	this routine will take the second derivative matrix in condensed form
c	as included in the vector diag, d2bond, d2nb_1, d2nb_2, d2nb_3, 
c	d2spec, d2theta, d2phi and d2imp, mass weight it and send it back
	write(*,*)'we are in massweig'



c	do the diag vector
	do 100 idiag=1,npt
		id = 6*(idiag-1) + 1
		sqtmp = 1.d0/ptms(idiag)
		diag(id)   = diag(id)*sqtmp
		diag(id+1) = diag(id+1)*sqtmp
		diag(id+2) = diag(id+2)*sqtmp
		diag(id+3) = diag(id+3)*sqtmp
		diag(id+4) = diag(id+4)*sqtmp
		diag(id+5) = diag(id+5)*sqtmp
100	continue
	write(*,*)'out of diag'

c	now the bonds
	do 200 ibond=1,nb
		ib = 6*(ibond-1) + 1
		i=ib1(ibond)
		j=ib2(ibond)
		sqtmp = sqrt(ptms(i)*ptms(j))
		sqtmp = 1.d0/sqtmp
		d2bond(ib)   = d2bond(ib)*sqtmp
		d2bond(ib+1) = d2bond(ib+1)*sqtmp
		d2bond(ib+2) = d2bond(ib+2)*sqtmp
		d2bond(ib+3) = d2bond(ib+3)*sqtmp
		d2bond(ib+4) = d2bond(ib+4)*sqtmp
		d2bond(ib+5) = d2bond(ib+5)*sqtmp
200	continue
	write(*,*)'out of bonds'


c	now the angles
	do 300 itheta=1,nangl
		ith = 27*(itheta-1) + 1
		i = iangl1(itheta)
		j = iangl2(itheta)
		k = iangl3(itheta)

c	from element ith to ith+8 pair i-j
		sqtmp = sqrt(ptms(i)*ptms(j))
		sqtmp = 1.d0/sqtmp
		d2theta(ith)    = d2theta(ith)*sqtmp
		d2theta(ith+1)  = d2theta(ith+1)*sqtmp
		d2theta(ith+2)  = d2theta(ith+2)*sqtmp
		d2theta(ith+3)  = d2theta(ith+3)*sqtmp
		d2theta(ith+4)  = d2theta(ith+4)*sqtmp
		d2theta(ith+5)  = d2theta(ith+5)*sqtmp
		d2theta(ith+6)  = d2theta(ith+6)*sqtmp
		d2theta(ith+7)  = d2theta(ith+7)*sqtmp
		d2theta(ith+8)  = d2theta(ith+8)*sqtmp

		

c	from element ith+9 to ith+17 pair i-k
		sqtmp = sqrt(ptms(i)*ptms(k))
		sqtmp = 1.d0/sqtmp
		d2theta(ith+9)   = d2theta(ith+9)*sqtmp
		d2theta(ith+10)  = d2theta(ith+10)*sqtmp
		d2theta(ith+11)  = d2theta(ith+11)*sqtmp
		d2theta(ith+12)  = d2theta(ith+12)*sqtmp
		d2theta(ith+13)  = d2theta(ith+13)*sqtmp
		d2theta(ith+14)  = d2theta(ith+14)*sqtmp
		d2theta(ith+15)  = d2theta(ith+15)*sqtmp
		d2theta(ith+16)  = d2theta(ith+16)*sqtmp
		d2theta(ith+17)  = d2theta(ith+17)*sqtmp



c	from element ith+18 to ith+26 pair j-k
		sqtmp = sqrt(ptms(j)*ptms(k))
		sqtmp = 1.d0/sqtmp
		d2theta(ith+18)  = d2theta(ith+18)*sqtmp
		d2theta(ith+19)  = d2theta(ith+19)*sqtmp
		d2theta(ith+20)  = d2theta(ith+20)*sqtmp
		d2theta(ith+21)  = d2theta(ith+21)*sqtmp
		d2theta(ith+22)  = d2theta(ith+22)*sqtmp
		d2theta(ith+23)  = d2theta(ith+23)*sqtmp
		d2theta(ith+24)  = d2theta(ith+24)*sqtmp
		d2theta(ith+25)  = d2theta(ith+25)*sqtmp
		d2theta(ith+26)  = d2theta(ith+26)*sqtmp


300	continue
	write(*,*)'out of angles'



c	now the proper torsions
	do 400 itors=1,ntors
		iphi = 54*(itors-1) + 1
		k = itor1(itors)
		l = itor2(itors)
		m = itor3(itors)
		n = itor4(itors)

c	from element iphi to iphi+8 pair k-l
		sqtmp = sqrt(ptms(k)*ptms(l))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi)    = d2phi(iphi)*sqtmp
		d2phi(iphi+1)  = d2phi(iphi+1)*sqtmp
		d2phi(iphi+2)  = d2phi(iphi+2)*sqtmp
		d2phi(iphi+3)  = d2phi(iphi+3)*sqtmp
		d2phi(iphi+4)  = d2phi(iphi+4)*sqtmp
		d2phi(iphi+5)  = d2phi(iphi+5)*sqtmp
		d2phi(iphi+6)  = d2phi(iphi+6)*sqtmp
		d2phi(iphi+7)  = d2phi(iphi+7)*sqtmp
		d2phi(iphi+8)  = d2phi(iphi+8)*sqtmp

c	from element iphi+9 to iphi+17 pair k-m
		sqtmp = sqrt(ptms(k)*ptms(m))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi+9)   = d2phi(iphi+9)*sqtmp
		d2phi(iphi+10)  = d2phi(iphi+10)*sqtmp
		d2phi(iphi+11)  = d2phi(iphi+11)*sqtmp
		d2phi(iphi+12)  = d2phi(iphi+12)*sqtmp
		d2phi(iphi+13)  = d2phi(iphi+13)*sqtmp
		d2phi(iphi+14)  = d2phi(iphi+14)*sqtmp
		d2phi(iphi+15)  = d2phi(iphi+15)*sqtmp
		d2phi(iphi+16)  = d2phi(iphi+16)*sqtmp
		d2phi(iphi+17)  = d2phi(iphi+17)*sqtmp



c	from element iphi+18 to iphi+26 pair k-n 
		sqtmp = sqrt(ptms(k)*ptms(n))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi+18)  = d2phi(iphi+18)*sqtmp
		d2phi(iphi+19)  = d2phi(iphi+19)*sqtmp
		d2phi(iphi+20)  = d2phi(iphi+20)*sqtmp
		d2phi(iphi+21)  = d2phi(iphi+21)*sqtmp
		d2phi(iphi+22)  = d2phi(iphi+22)*sqtmp
		d2phi(iphi+23)  = d2phi(iphi+23)*sqtmp
		d2phi(iphi+24)  = d2phi(iphi+24)*sqtmp
		d2phi(iphi+25)  = d2phi(iphi+25)*sqtmp
		d2phi(iphi+26)  = d2phi(iphi+26)*sqtmp


c	from element iphi+27 to iphi+35 pair l-m 
		sqtmp = sqrt(ptms(l)*ptms(m))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi+27)  = d2phi(iphi+27)*sqtmp
		d2phi(iphi+28)  = d2phi(iphi+28)*sqtmp
		d2phi(iphi+29)  = d2phi(iphi+29)*sqtmp
		d2phi(iphi+30)  = d2phi(iphi+30)*sqtmp
		d2phi(iphi+31)  = d2phi(iphi+31)*sqtmp
		d2phi(iphi+32)  = d2phi(iphi+32)*sqtmp
		d2phi(iphi+33)  = d2phi(iphi+33)*sqtmp
		d2phi(iphi+34)  = d2phi(iphi+34)*sqtmp
		d2phi(iphi+35)  = d2phi(iphi+35)*sqtmp

c	from element iphi+36 to iphi+44 pair l-n 
		sqtmp = sqrt(ptms(l)*ptms(n))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi+36)  = d2phi(iphi+36)*sqtmp
		d2phi(iphi+37)  = d2phi(iphi+37)*sqtmp
		d2phi(iphi+38)  = d2phi(iphi+38)*sqtmp
		d2phi(iphi+39)  = d2phi(iphi+39)*sqtmp
		d2phi(iphi+40)  = d2phi(iphi+40)*sqtmp
		d2phi(iphi+41)  = d2phi(iphi+41)*sqtmp
		d2phi(iphi+42)  = d2phi(iphi+42)*sqtmp
		d2phi(iphi+43)  = d2phi(iphi+43)*sqtmp
		d2phi(iphi+44)  = d2phi(iphi+44)*sqtmp

c	from element iphi+45 to iphi+53 pair m-n 
		sqtmp = sqrt(ptms(m)*ptms(n))
		sqtmp = 1.d0/sqtmp
		d2phi(iphi+45)  = d2phi(iphi+45)*sqtmp
		d2phi(iphi+46)  = d2phi(iphi+46)*sqtmp
		d2phi(iphi+47)  = d2phi(iphi+47)*sqtmp
		d2phi(iphi+48)  = d2phi(iphi+48)*sqtmp
		d2phi(iphi+49)  = d2phi(iphi+49)*sqtmp
		d2phi(iphi+50)  = d2phi(iphi+50)*sqtmp
		d2phi(iphi+51)  = d2phi(iphi+51)*sqtmp
		d2phi(iphi+52)  = d2phi(iphi+52)*sqtmp
		d2phi(iphi+53)  = d2phi(iphi+53)*sqtmp

400	continue
	write(*,*)'out of proper torsions'


c	now the improper torsions
	do 500 itors=1,nimp
		iphi = 54*(itors-1) + 1
		k = iimp1(itors)
		l = iimp2(itors)
		m = iimp3(itors)
		n = iimp4(itors)

c	from element iphi to iphi+8 pair k-l
		sqtmp = sqrt(ptms(k)*ptms(l))
		sqtmp = 1.d0/sqtmp
		d2imp(iphi)    = d2imp(iphi)*sqtmp
		d2imp(iphi+1)  = d2imp(iphi+1)*sqtmp
		d2imp(iphi+2)  = d2imp(iphi+2)*sqtmp
		d2imp(iphi+3)  = d2imp(iphi+3)*sqtmp
		d2imp(iphi+4)  = d2imp(iphi+4)*sqtmp
		d2imp(iphi+5)  = d2imp(iphi+5)*sqtmp
		d2imp(iphi+6)  = d2imp(iphi+6)*sqtmp
		d2imp(iphi+7)  = d2imp(iphi+7)*sqtmp
		d2imp(iphi+8)  = d2imp(iphi+8)*sqtmp

		

c	from element iphi+9 to iphi+17 pair k-m
		sqtmp = sqrt(ptms(k)*ptms(m))
		sqtmp = 1.d0/sqtmp
		d2imp(iphi+9)   = d2imp(iphi+9)*sqtmp
		d2imp(iphi+10)  = d2imp(iphi+10)*sqtmp
		d2imp(iphi+11)  = d2imp(iphi+11)*sqtmp
		d2imp(iphi+12)  = d2imp(iphi+12)*sqtmp
		d2imp(iphi+13)  = d2imp(iphi+13)*sqtmp
		d2imp(iphi+14)  = d2imp(iphi+14)*sqtmp
		d2imp(iphi+15)  = d2imp(iphi+15)*sqtmp
		d2imp(iphi+16)  = d2imp(iphi+16)*sqtmp
		d2imp(iphi+17)  = d2imp(iphi+17)*sqtmp



c	from element iphi+18 to iphi+26 pair k-n 
		sqtmp = sqrt(ptms(k)*ptms(n))
		sqtmp = 1.d0/sqtmp
		d2imp(iphi+18)  = d2imp(iphi+18)*sqtmp
		d2imp(iphi+19)  = d2imp(iphi+19)*sqtmp
		d2imp(iphi+20)  = d2imp(iphi+20)*sqtmp
		d2imp(iphi+21)  = d2imp(iphi+21)*sqtmp
		d2imp(iphi+22)  = d2imp(iphi+22)*sqtmp
		d2imp(iphi+23)  = d2imp(iphi+23)*sqtmp
		d2imp(iphi+24)  = d2imp(iphi+24)*sqtmp
		d2imp(iphi+25)  = d2imp(iphi+25)*sqtmp
		d2imp(iphi+26)  = d2imp(iphi+26)*sqtmp


c	from element iphi+27 to iphi+35 pair l-m 
		sqtmp = sqrt(ptms(l)*ptms(m))
		sqtmp = 1.d0/ sqtmp
		d2imp(iphi+27)  = d2imp(iphi+27)*sqtmp
		d2imp(iphi+28)  = d2imp(iphi+28)*sqtmp
		d2imp(iphi+29)  = d2imp(iphi+29)*sqtmp
		d2imp(iphi+30)  = d2imp(iphi+30)*sqtmp
		d2imp(iphi+31)  = d2imp(iphi+31)*sqtmp
		d2imp(iphi+32)  = d2imp(iphi+32)*sqtmp
		d2imp(iphi+33)  = d2imp(iphi+33)*sqtmp
		d2imp(iphi+34)  = d2imp(iphi+34)*sqtmp
		d2imp(iphi+35)  = d2imp(iphi+35)*sqtmp

c	from element iphi+36 to iphi+44 pair l-n 
		sqtmp= sqrt(ptms(l)*ptms(n))
		sqtmp= 1.d0/sqtmp
		d2imp(iphi+36)  = d2imp(iphi+36)*sqtmp
		d2imp(iphi+37)  = d2imp(iphi+37)*sqtmp
		d2imp(iphi+38)  = d2imp(iphi+38)*sqtmp
		d2imp(iphi+39)  = d2imp(iphi+39)*sqtmp
		d2imp(iphi+40)  = d2imp(iphi+40)*sqtmp
		d2imp(iphi+41)  = d2imp(iphi+41)*sqtmp
		d2imp(iphi+42)  = d2imp(iphi+42)*sqtmp
		d2imp(iphi+43)  = d2imp(iphi+43)*sqtmp
		d2imp(iphi+44)  = d2imp(iphi+44)*sqtmp

c	from element iphi+45 to iphi+53 pair m-n 
		sqtmp = sqrt(ptms(m)*ptms(n))
		sqtmp = 1.d0/sqtmp
		d2imp(iphi+45)  = d2imp(iphi+45)*sqtmp
		d2imp(iphi+46)  = d2imp(iphi+46)*sqtmp
		d2imp(iphi+47)  = d2imp(iphi+47)*sqtmp
		d2imp(iphi+48)  = d2imp(iphi+48)*sqtmp
		d2imp(iphi+49)  = d2imp(iphi+49)*sqtmp
		d2imp(iphi+50)  = d2imp(iphi+50)*sqtmp
		d2imp(iphi+51)  = d2imp(iphi+51)*sqtmp
		d2imp(iphi+52)  = d2imp(iphi+52)*sqtmp
		d2imp(iphi+53)  = d2imp(iphi+53)*sqtmp


500	continue
	write(*,*)'out of improper torsions'


c	now do non-bonded lists.
c	there are three of them, as you may already now.
c	they have basically the same structure but different sets of
c	pointers.
c	also done here are the special 1-4 interactions.


	jbeg = 1
	do 610 i=1,npt-1
		jend = point1(i)
		if(jbeg.le.jend) then
		do 710 k=jbeg,jend
			j = list1(k)
			kndex = 6*(k-1) + 1
			sqtmp = dsqrt(ptms(i)*ptms(j))
			sqtmp = 1.d0/sqtmp
			d2nb_1(kndex)   = d2nb_1(kndex)*sqtmp
			d2nb_1(kndex+1) = d2nb_1(kndex+1)*sqtmp
			d2nb_1(kndex+2) = d2nb_1(kndex+2)*sqtmp
			d2nb_1(kndex+3) = d2nb_1(kndex+3)*sqtmp
			d2nb_1(kndex+4) = d2nb_1(kndex+4)*sqtmp
			d2nb_1(kndex+5) = d2nb_1(kndex+5)*sqtmp
710		continue
		endif
		jbeg = jend + 1
610	continue


	jbeg = 1
	do 620 i=1,npt-1
		jend = point2(i)
		if(jbeg.le.jend) then
		do 720 k=jbeg,jend
			j = list2(k)
			kndex = 6*(k-1) + 1
			sqtmp = dsqrt(ptms(i)*ptms(j))
			sqtmp = 1.d0/sqtmp
			d2nb_2(kndex)   = d2nb_2(kndex)*sqtmp
			d2nb_2(kndex+1) = d2nb_2(kndex+1)*sqtmp
			d2nb_2(kndex+2) = d2nb_2(kndex+2)*sqtmp
			d2nb_2(kndex+3) = d2nb_2(kndex+3)*sqtmp
			d2nb_2(kndex+4) = d2nb_2(kndex+4)*sqtmp
			d2nb_2(kndex+5) = d2nb_2(kndex+5)*sqtmp
720		continue
		endif
		jbeg = jend + 1
620	continue


	jbeg = 1
	do 630 i=1,npt-1
		jend = point3(i)
		if(jbeg.le.jend) then
		do 730 k=jbeg,jend
			j = list3(k)
			kndex = 6*(k-1) + 1
			sqtmp = dsqrt(ptms(i)*ptms(j))
			sqtmp = 1.d0/sqtmp
			d2nb_3(kndex)   = d2nb_3(kndex)*sqtmp
			d2nb_3(kndex+1) = d2nb_3(kndex+1)*sqtmp
			d2nb_3(kndex+2) = d2nb_3(kndex+2)*sqtmp
			d2nb_3(kndex+3) = d2nb_3(kndex+3)*sqtmp
			d2nb_3(kndex+4) = d2nb_3(kndex+4)*sqtmp
			d2nb_3(kndex+5) = d2nb_3(kndex+5)*sqtmp
730		continue
		endif
		jbeg = jend + 1
630	continue


	jbeg = 1
	do 640 i=1,npt-1
		jend = spec1(i)
		if(jbeg.le.jend) then
		do 740 k=jbeg,jend
			j = spec2(k)
			kndex = 6*(k-1) + 1
			sqtmp = dsqrt(ptms(i)*ptms(j))
			sqtmp = 1.d0/sqtmp
			d2spec(kndex)   = d2spec(kndex)*sqtmp
			d2spec(kndex+1) = d2spec(kndex+1)*sqtmp
			d2spec(kndex+2) = d2spec(kndex+2)*sqtmp
			d2spec(kndex+3) = d2spec(kndex+3)*sqtmp
			d2spec(kndex+4) = d2spec(kndex+4)*sqtmp
			d2spec(kndex+5) = d2spec(kndex+5)*sqtmp
740		continue
		endif
		jbeg = jend + 1
640	continue




c	hopefully the end and we did not forgot anything	

	return
	end
