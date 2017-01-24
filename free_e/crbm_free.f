		subroutine crbm_free(coor,scalar,sigma,
     1		grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,divms,grdp,
     2          npt,nselec,pointr,istp,ntest,debug,udata)
c
c Factor out rigid body motion from the current coordinate system
c
	double precision coor(3,*)
	double precision scalar(7),sigma(7),grdcmx(3,*),grdcmy(3,*)
	double precision grdp(3,*)
	double precision grdcmz(3,*),grdlx(3,*),grdly(3,*),grdlz(3,*)
	double precision divms(*)
	integer pointr(*)
	integer npt,nselec,istp,ntest,udata
	logical debug
c
c r          - coordinates
c scalar     - work vector used to store scalar products
c sigma      - fixed value of constraints (i.e. ideally scalar should be
c		equal to sigma)
c grdcm[x-z] - gradient of center of mass constraints
c grdl[x-z]  - gradient of rigid infitesimal rotation
c divms      - 1/m vector
c grdp       -
c pointr     - a pointr to atoms which required the constraints
c		for a given structure
c npt      - number of atoms
c nselec     - number of selected atoms
c debug      - if .true. print debugging info.
c
c
c Note the coordinate constraints are of the form
c
c constraint = sum a(i)*r(i) + constant
c
c the lagrange multiplier (which we call here scalar) is
c scalar = sum a(i)*r(i) + constant
c
c 
c calculate the scalar product of the constraints gradient and
c the " free " positions
c
	integer i,j

		do 11 i=1,7
		 scalar(i)=0.d0
11		continue

		do 12 j=1,nselec
		 i=pointr(j)
		 scalar(1)=scalar(1)+grdp(1,j)*coor(1,i)+
     1			grdp(2,j)*coor(2,i)+grdp(3,j)*coor(3,i)
		 scalar(2)=scalar(2)+grdcmx(1,j)*coor(1,i)
		 scalar(3)=scalar(3)+grdcmy(2,j)*coor(2,i)
		 scalar(4)=scalar(4)+grdcmz(3,j)*coor(3,i)
		 scalar(5)=scalar(5)+grdlx(1,j)*coor(1,i)+
     1			grdlx(2,j)*coor(2,i)+grdlx(3,j)*coor(3,i)
		 scalar(6)=scalar(6)+grdly(1,j)*coor(1,i)+
     1			grdly(2,j)*coor(2,i)+grdly(3,j)*coor(3,i)
		 scalar(7)=scalar(7)+grdlz(1,j)*coor(1,i)+
     1			grdlz(2,j)*coor(2,i)+grdlz(3,j)*coor(3,i)
12		continue
c
c add the constants to the scalar products
c
c NOTE: the constraints over the path direction and the rigid body
c motions underwet orthonormalization with respect to mass
c weighting. 
c
		do 13 i=1,7
			scalar(i)=scalar(i)-sigma(i)
13		continue
c
c calculate the "corrected" coordinates (finally)
		do 14 j=1,nselec
		 i=pointr(j)

		 coor(1,i)=coor(1,i)-divms(i)*(scalar(1)*grdp(1,j)+
     1			scalar(2)*grdcmx(1,j)+scalar(5)*grdlx(1,j)+
     2			scalar(6)*grdly(1,j)+scalar(7)*grdlz(1,j))

		 coor(2,i)=coor(2,i)-divms(i)*(scalar(1)*grdp(2,j)+
     1			scalar(3)*grdcmy(2,j)+scalar(5)*grdlx(2,j)+
     2			scalar(6)*grdly(2,j)+scalar(7)*grdlz(2,j))

		 coor(3,i)=coor(3,i)-divms(i)*(scalar(1)*grdp(3,j)+
     1			scalar(4)*grdcmz(3,j)+scalar(5)*grdlx(3,j)+
     2			scalar(6)*grdly(3,j)+scalar(7)*grdlz(3,j))
14		continue
c
c test that the new coordinates satisfies the constraints
c
c	if (debug .or. istp/ntest*ntest.eq.istp) then
	if (debug) then
	  do 15 i=1,7
		scalar(i)=0.d0
15	  continue
	  do 16 j=1,nselec
		i=pointr(j)
		scalar(1)=scalar(1)+coor(1,i)*grdp(1,j)+coor(2,i)
     1			 *grdp(2,j)+coor(3,i)*grdp(3,j)
		scalar(2)=scalar(2)+coor(1,i)*grdcmx(1,j)
		scalar(3)=scalar(3)+coor(2,i)*grdcmy(2,j)
		scalar(4)=scalar(4)+coor(3,i)*grdcmz(3,j)
		scalar(5)=scalar(5)+grdlx(1,j)*coor(1,i)+
     1			grdlx(2,j)*coor(2,i)
     1			 +grdlx(3,j)*coor(3,i)
		scalar(6)=scalar(6)+grdly(1,j)*coor(1,i)+
     1			grdly(2,j)*coor(2,i)
     1			 +grdly(3,j)*coor(3,i)
		scalar(7)=scalar(7)+grdlz(1,j)*coor(1,i)+
     1			grdlz(2,j)*coor(2,i)
     1			 +grdlz(3,j)*coor(3,i)

16		continue
	  write(udata,165)(scalar(i)-sigma(i),i=1,7)
165	  format(//,1x,' errors in constraints ',/,1x,
     1     (e14.5),/,3(e14.5),//)
	  write(udata,*)' Calculated constraints',(scalar(i),i=1,7)
	  write(udata,*)' Exact values',(sigma(i),i=1,7)
	end if
100     continue
	return 
	end

