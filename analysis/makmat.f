	subroutine makmat(a,n)
	integer n
	double precision a(n,*)
	
	integer i,j

	do 2 i=2,n
		a(i,i) = 1.d0
2	continue
	
	return
	end
