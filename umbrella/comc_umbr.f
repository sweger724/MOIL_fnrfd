	subroutine comc_umbr(coor,grdp,divms,natom,pointr,
     1  nselec,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,dmass1,dmass2,
     2  coor1,sigma,debug)
c
c calculate gradient of rigid body constraints
c and then orthnormalize the constraints vectors.
c
c *** input
c coor     - the coordinates of the reference system for which
c            cenetr of mass translation and rotation are compared
c grdp     - the path slope at the position x0,y0,z0
c divms    - double precision vector of 1.d0 over masses
c debug    - logical. true=print a lot of debugging information
c
c *** output
c constraints value at initial configuration, should remain
c      		constants.
c grdcmx - gradient of the translation of center of mass (X)
c grdcmy                   and Y
c grdcmz                   and Z
c
c grdlx  - gradient of rotation (X)
c grdly                         (Y)
c grdlz                         (Z)
c
c NOTE THAT THE VECTORS ABOVE UNDERWENT ORTHONORMALIZTION WITH
c RESPCT TO EACH OTHER.
c
	include 'COMMON/UNITS.BLOCK'
	double precision coor(3,*)
	double precision grdp(3,*),grdcmx(3,*),grdcmy(3,*),grdcmz(3,*)
	double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
	double precision sigma(*)
	double precision divms(*),dmass1(*),dmass2(*),coor1(3,*)
	integer natom,nselec
	integer pointr(*)
	logical debug
c
c	local
c
	double precision norm,test(7)
	integer i,j,nsel3

	nsel3=3*nselec
c
c The zero order values of the constraints (which are not zero)
c are stored in sigma(1-7)
c sigma(1) - scalar product of path slope and initial coordinates
c sigma(2-4) - center of mass positions for x y z
c sigma(5-7) - Orientation value for x y & z ( zero before
c           orthonormalization)
c
	do 11 i=1,7
		sigma(i)=0.d0
11	continue
c
c calculate double precision mass and picking the
c selected atoms to a separate vector
c
	do 1 j=1,nselec
		i=pointr(j)
		dmass1(j) =1.d0/divms(i)
		dmass2(j) =divms(i)
		coor1(1,j)=coor(1,i)
		coor1(2,j)=coor(2,i)
		coor1(3,j)=coor(3,i)
1	continue
c
c calculate the gradient of the translation vector.
c Note the strange normalization: If a1(i) and a2(i) are
c vectors of constraint gradients and m(i) is the mass,
c the scalar product is defined by
c SUM a1(i)*a2(i)/m(i) 
c
	norm=0.d0
	do 2 j=1,nselec

		grdcmx(1,j)=dmass1(j)
		grdcmx(2,j)=0.d0
		grdcmx(3,j)=0.d0

		grdcmy(1,j)=0.d0
		grdcmy(2,j)=dmass1(j)
		grdcmy(3,j)=0.d0

		grdcmz(1,j)=0.d0
		grdcmz(2,j)=0.d0
		grdcmz(3,j)=dmass1(j)

		grdlx(1,j)=0.d0
		grdlx(2,j)=dmass1(j)*coor1(3,j)
		grdlx(3,j)=-dmass1(j)*coor1(2,j)

		grdly(1,j)=-grdlx(2,j)
		grdly(2,j)=0.d0
		grdly(3,j)=dmass1(j)*coor1(1,j)

		grdlz(1,j)=-grdlx(3,j)
		grdlz(2,j)=-grdly(3,j)
		grdlz(3,j)=0.d0
c
c actually it is: norm=norm+grdcm[xyz](i)*grdcm[xyz](i)/dmass(i)
c
		norm=norm+grdcmx(1,j)
c
c opportunity for constraints' values calculation
c
		sigma(1)=sigma(1)+coor1(1,j)*grdp(1,j)
     1			+coor1(2,j)*grdp(2,j)+coor1(3,j)*grdp(3,j)
		sigma(2)=sigma(2)+dmass1(j)*coor1(1,j)
		sigma(3)=sigma(3)+dmass1(j)*coor1(2,j)
		sigma(4)=sigma(4)+dmass1(j)*coor1(3,j)
2	continue
	norm=1.d0/dsqrt(norm)
c
c normalize grdcmx grdcmy grdcmz and the constraint functions
c
	sigma(2)=sigma(2)*norm
	sigma(3)=sigma(3)*norm
	sigma(4)=sigma(4)*norm
	do 3 i=1,nselec
	 do 3 j=1,3
		grdcmx(j,i)=grdcmx(j,i)*norm
		grdcmy(j,i)=grdcmy(j,i)*norm
		grdcmz(j,i)=grdcmz(j,i)*norm
3	continue
	if (debug) then
		do 31 i=1,7
			test(i)=0.d0
31		continue
		do 32 i=1,nselec
		 do 32 j=1,3
		  test(1)=test(1)+grdcmx(j,i)*grdcmx(j,i)*dmass2(i)
		  test(2)=test(2)+grdcmy(j,i)*grdcmy(j,i)*dmass2(i)
		  test(3)=test(3)+grdcmz(j,i)*grdcmz(j,i)*dmass2(i)
		  test(4)=test(4)+grdcmx(j,i)*grdlx(j,i)*dmass2(i)
		  test(5)=test(5)+grdcmx(j,i)*grdly(j,i)*dmass2(i)
		  test(6)=test(6)+grdcmx(j,i)*grdlz(j,i)*dmass2(i)
		  test(7)=test(7)+grdcmx(j,i)*grdp(j,i)*dmass2(i)
32		continue
		write(stdo,*)' before othog. norms of grdcm[xyz]'
		write(stdo,*)' scalar prod. of grdcmx & grdl[x-z] grdp'
		write(stdo,*)(test(i),i=1,7)
	end if
c
c orthogonalize the vectors with mass weighting. Obviously after
c orthonormaliztion, the meaning of rotation in a given direction
c may be lost.
c
c ** grdp
	call orthg(grdcmx,grdp,nselec,dmass2,sigma,2,1)
	call orthg(grdcmy,grdp,nselec,dmass2,sigma,3,1)
	call orthg(grdcmz,grdp,nselec,dmass2,sigma,4,1)
c ** grdlx
	call orthg(grdp  ,grdlx,nselec,dmass2,sigma,1,5)
	call orthg(grdcmx,grdlx,nselec,dmass2,sigma,2,5)
	call orthg(grdcmy,grdlx,nselec,dmass2,sigma,3,5)
	call orthg(grdcmz,grdlx,nselec,dmass2,sigma,4,5)
c ** grdly
	call orthg(grdp  ,grdly,nselec,dmass2,sigma,1,6)
	call orthg(grdcmx,grdly,nselec,dmass2,sigma,2,6)
	call orthg(grdcmy,grdly,nselec,dmass2,sigma,3,6)
	call orthg(grdcmz,grdly,nselec,dmass2,sigma,4,6)
	call orthg(grdlx ,grdly,nselec,dmass2,sigma,5,6)
c ** grdlz
	call orthg(grdp  ,grdlz,nselec,dmass2,sigma,1,7)
	call orthg(grdcmx,grdlz,nselec,dmass2,sigma,2,7)
	call orthg(grdcmy,grdlz,nselec,dmass2,sigma,3,7)
	call orthg(grdcmz,grdlz,nselec,dmass2,sigma,4,7)
	call orthg(grdlx ,grdlz,nselec,dmass2,sigma,5,7)
	call orthg(grdly ,grdlz,nselec,dmass2,sigma,6,7)

c
c check that everything is orthonormal
c
	if (debug) then
	do 61 i=1,7
		test(i)=0.d0
61	continue
	do 7 i=1,nselec
	 do 7 j=1,3
		test(1)=test(1)+grdp(j,i)*grdp(j,i)*dmass2(i)
		test(2)=test(2)+grdcmx(j,i)*grdcmx(j,i)*dmass2(i)
		test(3)=test(3)+grdcmy(j,i)*grdcmy(j,i)*dmass2(i)
		test(4)=test(4)+grdcmz(j,i)*grdcmz(j,i)*dmass2(i)
		test(5)=test(5)+grdlx(j,i)*grdlx(j,i)*dmass2(i)
		test(6)=test(6)+grdly(j,i)*grdly(j,i)*dmass2(i)
		test(7)=test(7)+grdlz(j,i)*grdlz(j,i)*dmass2(i)
7	continue
	write(*,*)' normalization of pgrd grdcm[x,y,z] grdl[x,y,z] '
	write(*,*)(test(i),i=1,7)
	do 8 i=1,7
		test(i)=0.d0
8	continue
	do 9 i=1,nselec
	 do 9 j=1,3
		test(1)=test(1)+grdp(j,i)*grdcmx(j,i)*dmass2(i)
		test(2)=test(2)+grdp(j,i)*grdcmy(j,i)*dmass2(i)
		test(3)=test(3)+grdp(j,i)*grdcmz(j,i)*dmass2(i)
		test(4)=test(4)+grdlx(j,i)*grdcmx(j,i)*dmass2(i)
		test(5)=test(5)+grdlx(j,i)*grdcmy(j,i)*dmass2(i)
		test(6)=test(6)+grdlx(j,i)*grdcmz(j,i)*dmass2(i)
		test(7)=test(7)+grdlx(j,i)*grdp(j,i)*dmass2(i)
9	continue
	write(*,*)'scalar products (should be zero) of:'
	write(*,*)'grdp-grdcmx grdp-grdcmy grdp-grdcmz grdlx-grdcmx'
	write(*,*)'grdlx-grdcmy grdlx-grdcmz grdlx-grdp'
	write(*,*)(test(i),i=1,7)
	end if
c --- end of debug statement
	return
	end
