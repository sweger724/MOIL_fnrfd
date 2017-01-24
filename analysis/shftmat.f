	subroutine shftmat(a,b,n,shft)

	integer n,shft
	double precision a(n,*),b(n+shft,*)

	integer i,j

	do 1 i=1,n
	 do 1 j=1,n
	  a(i,j) = b(i+shft,j+shft)
1	continue
	return
	end
