	subroutine orthmat(a,n)
	integer n
	double precision a(n,*)

	integer i,j,k
	double precision tmp

	do 1 i=1,n-1
	 do 1 j=i+1,n
	  call ortho(a(1,i),a(1,j),n)
1	continue
	
c test that the vectors are orthonormal
c@
c	do 3 i=1,n
c	 do 3 j=1,n
c	  tmp = 0.d0
c	  do 2 k=1,n
c		tmp = tmp + a(k,j)*a(k,i)
c2	continue
c	write(*,*)' i j <i|j> ',i,j,tmp
c3	continue
	return
	end
