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
	double precision scal,norm,tmp
	integer i,j

	scal=0.d0
	do 1 i=1,n
		tmp = X(1,i)*Y(1,i) + X(2,i)*Y(2,i) + X(3,i)*Y(3,i)
		scal = scal + tmp*amovr(i)
1	continue
	sigma(n2)=sigma(n2)-scal*sigma(n1)
	norm=0.d0
	do 2 i=1,n
		Y(1,i)=Y(1,i)-scal*X(1,i)
		Y(2,i)=Y(2,i)-scal*X(2,i)
		Y(3,i)=Y(3,i)-scal*X(3,i)
		tmp =  Y(1,i)*Y(1,i) + Y(2,i)*Y(2,i) + Y(3,i)*Y(3,i)
		norm = norm + tmp*amovr(i)
2	continue
	norm=1.d0/dsqrt(norm)
	sigma(n2)=sigma(n2)*norm
	do 3 i=1,n
	 do 3 j=1,3
		Y(j,i)=Y(j,i)*norm
3	continue
	return
	end
