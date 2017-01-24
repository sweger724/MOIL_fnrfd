	subroutine vec_mul_add(v1,v2,scalar,length)
c
c Multiplying the vector v2 by the scalar scalar and adding to v1 (output)
c
	double precision v1(*),v2(*),scalar
	integer length

c local
	integer i 

	do 1 i=1,length
		v1(i) = v1(i) + scalar * v2(i)
1	continue

	return
	end
