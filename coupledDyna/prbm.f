	subroutine prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)
c
c calculate gradient of rigid body rotation
c and then project out the velocity component that is parallel 
c to the constraints
c

c
c *** input
c vx,vy,vz,x,y,z - velocities and coordinates
c ptms - mass vector
c npt  - mass vector
c
c *** output
c
c grdlx  - gradient of rotation (X)
c grdly                         (Y)
c grdlz                         (Z)
c vx, vy, vz - velocities without rotation (on output)
c
c NOTE THAT THE VECTORS ABOVE UNDERWENT ORTHONORMALIZTION WITH
c RESPCT TO EACH OTHER.
c
	double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
	double precision coor(3,*)
	double precision velo(3,*),ptms(*)
	integer npt
c
c	local
c
	double precision vscalar
	double precision norm,scalar(3)
	integer i,j,npt3

	npt3 = 3*npt
c
c correct center of mass velocity (set to zero)
c
	norm = 0.d0
	do 10 i=1,3
	 scalar(i) = 0.d0
10	continue
	do 11 i=1,npt
	 norm = norm + ptms(i)
	 scalar(1) = scalar(1) + ptms(i)*velo(1,i)
	 scalar(2) = scalar(2) + ptms(i)*velo(2,i)
	 scalar(3) = scalar(3) + ptms(i)*velo(3,i)
11	continue

	norm = 1.d0/norm

	scalar(1) = scalar(1)*norm
	scalar(2) = scalar(2)*norm
	scalar(3) = scalar(3)*norm

	do 12 i=1,npt
	 velo(1,i) = velo(1,i) - scalar(1)
	 velo(2,i) = velo(2,i) - scalar(2)
	 velo(3,i) = velo(3,i) - scalar(3)
12	continue


	do 1 j=1,npt

		grdlx(1,j) = 0.d0
		grdlx(2,j) = ptms(j)*coor(3,j)
		grdlx(3,j) = -ptms(j)*coor(2,j)

		grdly(1,j) = -grdlx(2,j)
		grdly(2,j) = 0.d0
		grdly(3,j) = ptms(j)*coor(1,j)

		grdlz(1,j) = -grdlx(3,j)
		grdlz(2,j) = -grdly(3,j)
		grdlz(3,j) = 0.d0
1	continue

c
c Use grdlx as a base vector for Gram-Schmidth, Normalize it first

	norm = 1.d0/dsqrt(vscalar(grdlx,grdlx,npt3))

	do 3 i=1,npt
	  grdlx(1,i) = grdlx(1,i)*norm
	  grdlx(2,i) = grdlx(2,i)*norm
	  grdlx(3,i) = grdlx(3,i)*norm
3	continue

	call ortho(grdlx,grdly,npt3)
	call ortho(grdlx,grdlz,npt3)
	call ortho(grdly,grdlz,npt3)

	scalar(1) = vscalar(velo,grdlx,npt3)
	scalar(2) = vscalar(velo,grdly,npt3)
	scalar(3) = vscalar(velo,grdlz,npt3)


	do 5 j=1,npt
	 velo(1,j) = velo(1,j) - scalar(1)*grdlx(1,j) - 
     1		scalar(2)*grdly(1,j) - scalar(3)*grdlz(1,j)
	 velo(2,j) = velo(2,j) - scalar(1)*grdlx(2,j) - 
     1		scalar(2)*grdly(2,j) - scalar(3)*grdlz(2,j)
	 velo(3,j) = velo(3,j) - scalar(1)*grdlx(3,j) - 
     1		scalar(2)*grdly(3,j) - scalar(3)*grdlz(3,j)
5	continue

	return
	end

	subroutine prbm_prll(velo,coor,ptms,grdlx,grdly,grdlz,
     1                         npt,pt_start,pt_end)
c
c calculate gradient of rigid body rotation
c and then project out the velocity component that is parallel 
c to the constraints
c

c
c *** input
c vx,vy,vz,x,y,z - velocities and coordinates
c ptms - mass vector
c npt  - mass vector
c
c *** output
c
c grdlx  - gradient of rotation (X)
c grdly                         (Y)
c grdlz                         (Z)
c vx, vy, vz - velocities without rotation (on output)
c
c NOTE THAT THE VECTORS ABOVE UNDERWENT ORTHONORMALIZTION WITH
c RESPCT TO EACH OTHER.
c
	double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
	double precision coor(3,*)
	double precision velo(3,*),ptms(*)
	integer npt,pt_start,pt_end
c
c	local
c
	double precision vscalar_prll
	double precision norm,scalar(3),temp
	integer i,j,npt3,pt_start3,pt_end3

	npt3 = 3*npt
	pt_start3 = (pt_start - 1)*3 + 1
	pt_end3 = pt_end * 3		
c
c correct center of mass velocity (set to zero)
c
	norm = 0.d0
	do 10 i=1,3
	 scalar(i) = 0.d0
10	continue
	do 11 i=pt_start,pt_end
	 norm = norm + ptms(i)
	 scalar(1) = scalar(1) + ptms(i)*velo(1,i)
	 scalar(2) = scalar(2) + ptms(i)*velo(2,i)
	 scalar(3) = scalar(3) + ptms(i)*velo(3,i)
11	continue

	call reduce_1(norm)
	call reduce_1(scalar(1))
	call reduce_1(scalar(2))   
	call reduce_1(scalar(3)) 

	norm = 1.d0/norm

	scalar(1) = scalar(1)*norm
	scalar(2) = scalar(2)*norm
	scalar(3) = scalar(3)*norm

	do 12 i=pt_start,pt_end
	 velo(1,i) = velo(1,i) - scalar(1)
	 velo(2,i) = velo(2,i) - scalar(2)
	 velo(3,i) = velo(3,i) - scalar(3)
12	continue


	do 1 j=1,npt

		grdlx(1,j) = 0.d0
		grdlx(2,j) = ptms(j)*coor(3,j)
		grdlx(3,j) = -ptms(j)*coor(2,j)

		grdly(1,j) = -grdlx(2,j)
		grdly(2,j) = 0.d0
		grdly(3,j) = ptms(j)*coor(1,j)

		grdlz(1,j) = -grdlx(3,j)
		grdlz(2,j) = -grdly(3,j)
		grdlz(3,j) = 0.d0
1	continue

c
c Use grdlx as a base vector for Gram-Schmidth, Normalize it first

	temp = vscalar_prll(grdlx,grdlx,pt_start3,pt_end3)
	call reduce_1(temp)

	norm = 1.d0/dsqrt(temp)

	do 3 i=1,npt
	  grdlx(1,i) = grdlx(1,i)*norm
	  grdlx(2,i) = grdlx(2,i)*norm
	  grdlx(3,i) = grdlx(3,i)*norm
3	continue

	call ortho(grdlx,grdly,npt3)
	call ortho(grdlx,grdlz,npt3)
	call ortho(grdly,grdlz,npt3)

	scalar(1) = vscalar_prll(velo,grdlx,pt_start3,pt_end3)
	scalar(2) = vscalar_prll(velo,grdly,pt_start3,pt_end3)
	scalar(3) = vscalar_prll(velo,grdlz,pt_start3,pt_end3)

	call reduce_1(scalar(1))
	call reduce_1(scalar(2))
	call reduce_1(scalar(3))

	do 5 j=pt_start,pt_end
	 velo(1,j) = velo(1,j) - scalar(1)*grdlx(1,j) - 
     1		scalar(2)*grdly(1,j) - scalar(3)*grdlz(1,j)
	 velo(2,j) = velo(2,j) - scalar(1)*grdlx(2,j) - 
     1		scalar(2)*grdly(2,j) - scalar(3)*grdlz(2,j)
	 velo(3,j) = velo(3,j) - scalar(1)*grdlx(3,j) - 
     1		scalar(2)*grdly(3,j) - scalar(3)*grdlz(3,j)
5	continue

	return
	end
