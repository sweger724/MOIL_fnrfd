	subroutine scalar_mu_vec2(scalar, vec1, vec2, length)
c
c multiplying vector by a scalar. On output the vector is changed
c
	double precision vec1(*),vec2(*)
	double precision scalar
	integer length

c local
	integer i
	
	do 1 i = 1,length
		vec2(i) = vec1(i) * scalar
1	continue

	return
	end
