	subroutine vadd(a,b,length)
c
c add b to a
c
	double precision a(*),b(*)
	integer length

	integer i

	do 1 i=1,length
	 a(i) = a(i) + b(i)
1	continue

	return
	end
