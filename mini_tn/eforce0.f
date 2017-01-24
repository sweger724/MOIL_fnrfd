	subroutine eforce0(idim,xx,etotal,g,efcall)
c front end routine for energy call. Used by the minimizer

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'

	integer idim
	double precision etotal
	double precision xx(*),g(*)
	integer efcall
	integer i

	efcall = efcall + 1

	etotal  = 0.d0
	e_total = 0.d0
	do 1 i=1,npt
	 coor(1,i)  = xx(i)
	 coor(2,i)  = xx(i+npt)
	 coor(3,i)  = xx(i+2*npt)
1	continue
	call eforce()

	etotal = e_total

	do 9 i=1,npt
	 g(i)       = dpot(1,i)
	 g(i+npt)   = dpot(2,i)
	 g(i+2*npt) = dpot(3,i)
9	continue

	return
	end
