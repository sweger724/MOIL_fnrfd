	subroutine offdm(a,b,tmp,ptms,npt)
	integer npt
	double precision a(3*npt,*)
	double precision b(3*npt,*)
	double precision tmp(3*npt,*)
	double precision ptms(*)

	integer i,j,k,l

c compute off-diagonal mass matrix
c Let a(i,j) be the rotation matrix and ptms(i)
c is the vecotr with the diagonal elements of the mass matrix
c compute b(i,j) = a()ptm()a+()
c tmp is a buffer matrix
c
	do 1 i=1,3*npt
	 do 1 j=1,3*npt
	  b(i,j) = 0.d0
	  do 1 k=1,3*npt
	   l = (k-1)/3+1
	   b(i,j) = b(i,j) + a(k,i)*a(k,j)*ptms(l)
1	continue
	 
	return
	end
