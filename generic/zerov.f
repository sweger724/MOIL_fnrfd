	subroutine zerov(vector,length)
	double precision vector(*)
	integer length
	integer i

	do 1 i=1,length
		vector(i) =0.d0
1	continue
	return
	end
