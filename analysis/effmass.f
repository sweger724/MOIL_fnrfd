	subroutine effmass(rc_mass,inv2,invmass,n)
c
c compute the inverse of the effective mass along the reaction
c coordinbate - rc_mass
c
	integer n
	double precision rc_mass,inv2(n-1,*)
	double precision invmass(n,*)

c local
	integer i,j

	rc_mass = invmass(1,1)
	write(*,*)' invmass(1,1) = ',invmass(1,1)
	do 1 i=1,n-1
	 do 1 j=1,n-1
	  rc_mass = rc_mass - invmass(1,i+1)*inv2(i,j)*invmass(1,j+1)
1	continue
	rc_mass = 1.d0/rc_mass
	return
	end
