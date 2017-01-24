         subroutine num_drv()
c
c a branch routine to accumulate the different energies and forces
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/LINE.BLOCK'
c	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/TETHER.BLOCK'
	include 'COMMON/VELOC.BLOCK'

	double precision delta , e_old
	integer i,j

	delta = 1.d-7
	e_total = 0.d0
	do 3 i=1,npt
	 velo(1,i) = 0.d0
	 velo(2,i) = 0.d0
	 velo(3,i) = 0.d0
3	continue

	do 1 i=1,npt
	 do 1 j=1,3
	   coor(j,i) = coor(j,i) + delta
	   call eforce()
c@	   call wener(stdo)
           e_old = e_total
	   coor(j,i) = coor(j,i) - 2.d0*delta
	   call eforce()
c@	   call wener(stdo)
	   coor2(j,i) = (e_old-e_total)/(2.d0*delta)
	   coor(j,i) = coor(j,i) + delta
1	continue
	   
	return
	end
