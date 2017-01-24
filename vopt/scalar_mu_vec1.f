	subroutine scalar_mu_vec1(scalar, vec, length)
c
c multiplying vector by a scalar. On output the vector is changed
c
	double precision vec(*)
	double precision scalar
	integer length

c local
	integer i
	
	do 1 i = 1,length
		vec(i) = vec(i) * scalar
1	continue

	return
	end
