	subroutine vinit(vector,scalar,length)
c
c initiatle all the length elements of a vector to a scalar (typically zero)
c
	double precision vector(*)
	double precision scalar
	integer length

c	local
	integer i
	
	do 1 i=1,length
		vector(i) = scalar
1	continue

	return	
	end
