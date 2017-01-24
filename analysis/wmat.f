	subroutine wmat(a,n)
	integer n
	double precision a(n,*)

	integer i,j

	do 1 i=1,n
	 write(*,2)(a(j,i),j=1,n)
1	continue
2	format(1x,10f7.3)

	return
	end
