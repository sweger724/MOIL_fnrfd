	subroutine Find_Center(npt,coor,center)
        implicit none
c
C    Calculates the center of mass of coordinates coor
C    assumes the same masses, all = 1 
C
C

	double precision coor(3,*),center(3)
	integer npt

c local
	integer i,l
	
c calculate center of mass 
c
	do l=1,3
	  center(l)=0.d0
	enddo

	do i = 1, npt
	  do l = 1,3
	    center(l) = center(l) + coor(l,i)
	  end do
	end do

	do l=1,3
          center(l)= center(l)/npt
        enddo

C	write(6,*)"Center:",center(1),center(2),center(3)
	return
	end

