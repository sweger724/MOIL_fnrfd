	Subroutine eball()

C	Constraining existing water molecules to a sphere radius rball
C	centered at rcball (rcball is a vector of 3)
C

      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/EBALL.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'

	double precision xo,yo,zo,r2,r,s,e,de,dr,des

	integer i

c compute the constant shift of the energy
	
	e_ball = 0.d0
	rball2 = rball*rball
	do 400 i=1,npt
		xo  = coor(1,i) - rcball(1)
		yo  = coor(2,i) - rcball(2)
		zo  = coor(3,i) - rcball(3)
		r2   = (xo*xo+yo*yo+zo*zo)
		r  = dsqrt(r2)
c		write(*,*)' i r xo yo zo'
c		write(*,*)i,r,xo,yo,zo
		if (r2.gt.rball2) then
			r  = dsqrt(r2)
			s  = 1.d0/r
			dr = r - rball
			de = fball*dr
			des = de*s
			dpot(1,i) = dpot(1,i) + des*xo
			dpot(2,i) = dpot(2,i) + des*yo
			dpot(3,i) = dpot(3,i) + des*zo
			e = 0.5d0*de*dr
c			write(*,*)' i e dr r des xo yo zo'
c			write(*,*)i,e,dr,r,des,xo,yo,zo
			e_ball = e_ball + e
		end if
400	continue
	return
	end
