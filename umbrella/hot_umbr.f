	function hot_umbr(velo,dmass,npt,ndegf,massw)
c
c calculate how "hot" the situation is (the temperature)
c vx vy vz - velocities
c npt    - number of atoms (dimension of the vectors)
c ndegf    - number of degrees of freedom (taking of constraints)
c
	double precision hot_umbr
	double precision velo(3,*)
	double precision dmass(*)
	integer npt,ndegf
	logical massw
c
c local
c
	integer i
	double precision factor
	data factor/1.98768d-3/

	hot_umbr=0.d0
	if (massw) then
	do 1 i=1,npt
		hot_umbr=hot_umbr+dmass(i)*(velo(1,i)*velo(1,i)
     1			+velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
1	continue
	else
	do 2 i=1,npt
		hot_umbr=hot_umbr+velo(1,i)*velo(1,i)+
     9 	   	velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i)
2	continue
	end if
	hot_umbr=hot_umbr/(factor*ndegf)
	return
	end
