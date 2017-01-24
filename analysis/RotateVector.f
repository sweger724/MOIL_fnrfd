	subroutine RotateVector(vector,matrix)
        implicit none
c
C    Rotates vector by a rotation defined by 3x3 matrix
C
	double precision vector(3),matrix(3,3)

c local
	integer i,l
	double precision tmp(3)
	
c rotate vector 
c
	do l=1,3
	  tmp(l) =  matrix(l,1)*vector(1) + matrix(l,2)*vector(2)
     1          + matrix(l,3)*vector(3)
	enddo

	do l=1,3
	  vector(l)=tmp(l)
	end do

	return
	end

