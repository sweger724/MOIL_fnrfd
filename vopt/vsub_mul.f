	subroutine vsub_mul(x,y,f,length)
	double precision x(*),y(*)
	double precision f
	integer length

c local
	integer i

	do 1 i=1,length
		x(i) = f*(x(i) - y(i))
1	continue

	return
	end
