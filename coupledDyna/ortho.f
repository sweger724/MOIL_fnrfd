	subroutine ortho(X,Y,n)
c 
c Orthonormalization of two vectors. X is assumed to be normalized
c Y on output is normalized.
c
	double precision X(*),Y(*)
	integer n
c
c	local
c
	double precision scal,norm
	integer i

	scal=0.d0
	do 1 i=1,n
		scal=scal+X(i)*Y(i)
1	continue
	norm=0.d0
	do 2 i=1,n
		Y(i)=Y(i)-scal*X(i)
		norm=norm+Y(i)*Y(i)
2	continue
	norm=1.d0/dsqrt(norm)
	do 3 i=1,n
		Y(i)=Y(i)*norm
3	continue
	return
	end
