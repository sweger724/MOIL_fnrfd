	subroutine orthg(X,Y,n,amovr,sigma,n1,n2)
c 
c X is assumed to be normalized, Y on output is normalized.
c the orthnormalization includes also the scalar functions sigma
c
	double precision X(*),Y(*),amovr(*),sigma(*)
	integer n,n1,n2
c
c	local
c
	double precision scal,norm
	integer i,j,nselec

	nselec=n/3
	scal=0.d0
	do 1 i=1,n
		j=i-((i-1)/nselec)*nselec
		scal=scal+X(i)*Y(i)*amovr(j)
1	continue
	sigma(n2)=sigma(n2)-scal*sigma(n1)
	norm=0.d0
	do 2 i=1,n
		j=i-((i-1)/nselec)*nselec
		Y(i)=Y(i)-scal*X(i)
		norm=norm+Y(i)*Y(i)*amovr(j)
2	continue
	norm=1.d0/dsqrt(norm)
	sigma(n2)=sigma(n2)*norm
	do 3 i=1,n
		Y(i)=Y(i)*norm
3	continue
	return
	end
