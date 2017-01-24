	subroutine sds(sener,dv,r,d0,e0,e1,nselec,pointr
     1		,ipick,npt,ndegf,gamma,rho,lambda,fixend
     2		,debug)
	implicit none
c
c  a subroutine to calculate the line integral using trapezoidal rule with a constraint
c of points equally distributed along the path
c

c  The polymer energy (sener) is giving by
c
c
c the second and the point before the last are somewhat special since they need 
c delta_first and delta_last
c

c first active structure (second point)
c
c Contribution from three parts:

c   S = 

c  FROM STRUCTURES THAT ARE NOT THE FIRST OR LAST
c
c    0.5 *   SUM { |grad(V(i))| + |grad(V(i+1))| }  d(i,i+1) +    gamma SUM (d(i,i+1) - <d>)^2 + 
c
c From the first point
c          + 0.5 * |grad(V(1))| d(1,2) + gamma (d(1,2)-<d>)^2 +

c From the last point
c	   + 0.5 * |grad(V(n)| d(n-1,n) + gamma (d(n-1,n)-<d>)^2
c
c Note that the norm2(gradV) is normalized with respect to the number
c of degrees of freedom
c

	double precision sener,gamma,rho,lambda,gam2
	integer nselec,npt,npt3,ndegf

	double precision dv(3,*),r(3,*),d0(*),e0(*),e1(*)
	integer pointr(*),ipick(*)
	logical fixend,debug

c
c common block for COORDinates and potential ENERGY derivatives
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SDP.BLOCK'

C      SGB MODIFICATIONS
	include 'COMMON/SGB.BLOCK'
C      END SGB MODIFICATIONS
c
c Local
c
	integer igrid,npt2,i,j,k,l
	integer k1,k2
	integer namel
	character*3 name
	double precision tmp,dave
	double precision gradf,gradf_prev
	double precision delta(3,maxpt)
	double precision delta_prev(3,maxpt)
	double precision dpot_next(3,maxpt)
	double precision dpot_save(3,maxpt)
	double precision epsilon
	double precision gradf_first,gradf_last
	double precision d_dist
	double precision vnorm
	logical debug1
	logical first
	data first/.true./
	save gradf_first,gradf_last,first
c
c initialiaztion
c
	namel=3
	name='sds'
	epsilon = 1.d-6
	debug1 = .false.
	npt2  = 2*npt
	npt3  = 3*npt
	igrid = ndegf/(3*npt)
	gam2  = 2.d0*gamma
	sener = 0.d0
c
c initialize vectors
c
	call vinit(dv,0.d0,(igrid)*npt3)
	call vinit(d0,0.d0,igrid-1)

	if (first) then
c first call to sds compute forces at the first and last (non moving) coordinate set
c note that the nonbonded list is computed inefficiently here (for every step)
c
		call vdcopy(r(1,1),coor(1,1),npt3)
		if (esymyes) call squeeze()
		call nbondm()
		if (esymyes) call syminit
		call eforce()
		e0(1) = e_total
		call wener(stdo)
		call force_norm(gradf_first)
		e1(1) = gradf_first
c		write(*,*)' First call: gradf_first = ',gradf_first
		call vdcopy(r(1,(igrid-1)*npt+1),coor(1,1),npt3)
		if (esymyes) call squeeze()
		call nbondm()
		if (esymyes) call syminit
		call eforce()
		e0(igrid) = e_total
		call wener(stdo)
		call force_norm(gradf_last)
		e1(igrid) = gradf_last
		first = .false.
	end if
c
c computing the distances and the derivative from the first point to the
c second 
c
	call vecmin(r(1,npt+1),r(1,1),npt3,delta_prev)
c@
c	write(*,*)' gradf_first = ',gradf_first
c	write(*,*)' delta_prev '
c	do i = 1,npt
c		write(*,*)delta_prev(1,i),delta_prev(2,i),delta_prev(3,i)
c	end do
	d0(1) = vnorm(delta_prev,npt3)
c
c@
c	write(*,*)' d0(1) = ',d0(1)
	if (d0(1).lt.1.d-12) then
		write(*,*)' First distance: '
		call alert(name,namel,'Distance is 0!',15,0)
	end if
	dave = d0(1)
	sener = sener + gradf_first*d0(1)
c	write(*,*)' on first step sener = ',sener
	gradf_prev = gradf_first

c This loop is for configurations along the chain that are not first or last
c
	do 1 j=2,igrid-1
		k = (j-1)*npt+1
		call vec_mul_add(dv(1,k),delta_prev,
     1			gradf_prev/(d0(j-1)*npt3),npt3)

		call vdcopy(r(1,k),coor(1,1),npt3)
		if (esymyes) call squeeze()
		call nbondm()
		if (esymyes) call syminit()
     		call eforce()
		e0(j) = e_total
		call force_norm(gradf)
		e1(j) = gradf
c
c delta = r(k+npt) - r(k)
c
		call vecmin(r(1,npt+k),r(1,k),npt3,delta)
c
		d0(j) = vnorm(delta,npt3)
c@
cc	write(*,*) ' r(3) '
c	do i = 0,npt-1
c		write(*,*)r(1,i+npt+k),delta(2,i+npt+k),delta(3,i+npt+k)
c	end do
c	write(*,*) ' r(2) '
c	do i = 1,npt
c		write(*,*)r(1,i+npt),delta(2,i+npt),delta(3,i+npt)
c	end do
c	write(*,*)' delta ',j
c	do i = 1,npt
c		write(*,*)delta(1,i),delta(2,i),delta(3,i)
c	end do
c
c@
c	write(*,*)' d0(j) = ',d0(j)
c@
c	write(*,*) ' j delta = ',j,((delta(k1,k2),k1=1,3),k2=1,npt)
c		write(*,*)' d0(',j,' ) = ',d0(j)
		if (d0(j).lt.1.d-12) then
			write(*,*) ' distance ',j,' is zero!'
			call alert(name,namel,'Distance is 0!',15,0)
		end if
		sener = sener + gradf*(d0(j)+d0(j-1))
c		write(*,*)' addition to sener is now ',j,gradf*(d0(j)+d0(j-1))
c		write(*,*)' gradf d0(j) d0(j-1) ',gradf,d0(j),d0(j-1)
c
c distance derivatives
c
c
c delta_prev = r(k) - r(k-npt)
c
c the term computed here is |gradU(r(k))|*(delta_prev/|r(k)-r(k-npt)|)/npt3
c
		call vec_mul_add(dv(1,k),delta_prev,gradf/(d0(j-1)*npt3),npt3)
c
c the term computed here is |gradU(r(k))|*(delta/|r(k)-r(k+npt)|)/npt3
c
		call vec_mul_add(dv(1,k),delta,-gradf/(d0(j)*npt3),npt3)
c
c add the length derivative times the force to former derivative (j>2)
c
		if (j.gt.2) then
		 call vec_mul_add(dv(1,k-npt),delta_prev,
     1-gradf/(d0(j-1)*npt3),npt3)
		end if
c
c save delta in delta_prev for next round

		call vdcopy(delta,delta_prev,npt3)
c
c save gradf in gradf_first for next round
c
		gradf_prev = gradf
c
c force derivatives
c
		call vdcopy(dpot,dpot_save,npt3)
c
c compute gradU[ coor + epsilon* gradU/|gradU| ]
c
		tmp = epsilon / gradf
		call vec_mul_add(coor,dpot_save,tmp,npt3)
		call eforce()
		call vdcopy(dpot,dpot_next,npt3)
c
c compute gradU[ coor - epsilon*gradU/|gradU| ]
c
		call vec_mul_add(coor,dpot_save,-2.d0*tmp,npt3)
		call eforce()

c		write(*,*)' dpot_next dpot_prev'
c		do k1=1,npt
c			write(*,*)dpot_next(1,k1),dpot(1,k1)
c			write(*,*)dpot_next(2,k1),dpot(2,k1)
c			write(*,*)dpot_next(3,k1),dpot(3,k1)
c		end do
c		write(*,*)
		call vsub_mul(dpot_next,dpot,1.d0/(2.d0*tmp),npt3)
c		write(*,*)' d2pot*dpot '
c		do k1=1,npt
c			write(*,*)dpot_next(1,k1),dpot_next(2,k1),dpot_next(3,k1) 
c		end do
c END @
		call vec_mul_add(coor,dpot_save,tmp,npt3)
		tmp = (d0(j) + d0(j-1))/(npt3*gradf)
		call vec_mul_add(dv(1,k),dpot_next,tmp,npt3)
c
c average distance
c
		dave = dave + d0(j)
c		write(*,*) ' dave d0(j) j ',dave,d0(j),j
c		write(*,*)' end of loop j = ',j
1	continue
c	write(*,*) ' out of loopp one '
c
c computing the ocntribution of the last point 
c
	sener = sener + gradf_last*d0(igrid-1)
c	write(*,*)' last addition to sener ',gradf_last*d0(igrid-1)
c	write(*,*) ' sener = ',sener
	k = (igrid-2)*npt+1
c	write(*,*)' d0 = ',(d0(j),j=1,igrid-1)
	call vec_mul_add(dv(1,k),delta,-gradf_last/
     1(npt3*d0(igrid-1)),npt3)

	dave = dave/(igrid-1)
c	write(*,*) ' dave = ',dave

c
c compute the equi-distance constraint. There is some waste here since 
c delta-s are re-computed for the different igrid. However, I am not sure how
c to avoid it (without using a lot of memory) since dave is require. RE.
c
	tmp = 0.d0
	do 2 j=1,igrid-1
		tmp = tmp + (d0(j) - dave)**2
2	continue
	sener = sener + gamma*tmp

	call vecmin(r(1,1+npt),r(1,1),npt3,delta_prev)
c	write(*,*) ' dave = ',dave
c	write(*,*) ' before loopp two '
	do 3 j=2,igrid-1
	 k = npt*(j-1) + 1
	 d_dist = d0(j-1) - dave
	 tmp = 2.d0*gamma*d_dist/(d0(j-1)*npt3)
	 call vec_mul_add(dv(1,k),delta_prev,tmp,npt3)
	 d_dist = d0(j) - dave
	 call vecmin(r(1,k+npt),r(1,k),npt3,delta)
	 tmp = -2.d0*gamma*d_dist/(d0(j)*npt3)
	 call vec_mul_add(dv(1,k),delta,tmp,npt3)
	 call vdcopy(delta,delta_prev,npt3)
3	continue
	return
	end
