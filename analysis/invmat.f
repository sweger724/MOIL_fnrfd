	subroutine invmat(a,b,c,vec,n)
	integer n
	double precision a(n,*),b(n,*),c(n,*),vec(*)

	integer i,j,k
	double precision tmp

	do 1 i=1,n
	 do 1 j=1,n
	  b(i,j) = 0.d0
1	continue
	do 2 k=1,n
	 tmp = 1.d0/vec(k)
	 do 2 i=1,n
	  do 2 j=1,n
		b(i,j) = b(i,j) + a(i,k)*a(j,k)*tmp
2	continue

c@
c	write(*,*) 'the matrix c'
c	do 25 i=1,n
c	 write(*,100)(c(i,j),j=1,n)
c25	continue
c	write(*,*) ' the inverse '
c	do 26 i=1,n
c	 write(*,100)(b(i,j),j=1,n)
c26	continue
c100	format(1x,6(f9.3))
c check that b is an inverse to c
c	do 3 i=1,n
c	 do 3 j=1,n
c	  a(i,j) = 0.d0
c	  do 3 k=1,n
c	  a(i,j) = a(i,j) + b(i,k)*c(j,k)
c3	continue

c	write(*,*)' --- check inverse '
c	do 4 i=1,n
c		write(*,*)(a(i,j),j=1,n)
c4	continue

	return
	end
