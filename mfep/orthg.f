	subroutine orthg(X,Y,n,amovr,sigma,n1,n2)
c 
c X is assumed to be normalized, Y on output is normalized.
c the orthnormalization includes also the scalar functions sigma
c
	double precision X(3,*),Y(3,*),amovr(*),sigma(*)
	integer n,n1,n2
c
c	local
c
	double precision scal,norm
	integer i

	scal=0.d0
	do 1 i=1,n
		scal=scal+amovr(i)*(X(1,i)*Y(1,i)
     1			+X(2,i)*Y(2,i)+X(3,i)*Y(3,i))
1	continue
	sigma(n2)=sigma(n2)-scal*sigma(n1)
	norm=0.d0
	do 2 i=1,n
		Y(1,i)=Y(1,i)-scal*X(1,i)
		Y(2,i)=Y(2,i)-scal*X(2,i)
		Y(3,i)=Y(3,i)-scal*X(3,i)
		norm=norm+amovr(i)*(Y(1,i)*Y(1,i)+
     1			Y(2,i)*Y(2,i)+Y(3,i)*Y(3,i))
2	continue
	norm=1.d0/dsqrt(norm)
	sigma(n2)=sigma(n2)*norm
	do 3 i=1,n
		Y(1,i)=Y(1,i)*norm
		Y(2,i)=Y(2,i)*norm
		Y(3,i)=Y(3,i)*norm
3	continue
	return
	end
