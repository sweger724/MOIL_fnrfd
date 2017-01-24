	subroutine v_sub_mul_add(v1,v2,scalar,v3,length)
c
c Multiplying the vector v2 by the scalar scalar and adding to v1 (output)
c
	double precision v1(*),v2(*),v3(*),scalar

	integer length

c local
	integer i 

	do 1 i=1,length
		v3(i) = v3(i) + scalar*(v1(i)-v2(i))
1	continue

	return
	end
