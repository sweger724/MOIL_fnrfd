	function hotchaf(v,ptms,ndegf)
c
c calculate how "hot" the situation is (the temperature)
c vx vy vz - velocities
c natom    - number of atoms (dimension of the vectors)
c ndegf    - number of degrees of freedom (taking of constraints)
c
	double precision hotchaf
	double precision v(*)
	double precision ptms(*)
	integer ndegf
c
c local
c
	integer i
	double precision factor
	data factor/1.98768d-3/

	hotchaf=0.d0
	do 1 i=1,ndegf
	 hotchaf=hotchaf+v(i)*v(i)*ptms(i)
1	continue
	hotchaf=hotchaf/(factor*ndegf)
	return
	end
