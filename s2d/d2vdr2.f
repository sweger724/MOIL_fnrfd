	subroutine d2vdr2()
c
c calling specific second derivative routines.
c
	
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	include 'COMMON/SCNDRV.BLOCK'


	character*6 name
	integer namel,i
	data name/'d2vdr2'/
	data namel/6/

	do 22 i=1,6*npt
		diag(i) = 0.d0
22	continue

	if (ebyes) then
		if (debug) write(*,*)' Calling ebond2 '
		call bond2()
	end if
	if (ethyes) then
		if (debug) write(*,*)' calling theta2 '
		call theta2()
	end if
	if (etoyes) then
		if (debug) write(*,*)' calling etors2 '
		call etors2()
	end if
	if (eimyes) then 
		if (debug) write(*,*)' calling eimphi2 '
		call eimphi2()
	end if
	if (evdyes.or.eelyes) then
		if (debug) write(*,*)' calling cdie2 '
		call cdie2()
	end if
	if (nbeta.gt.0) then
		if (debug) write(*,*)' calling hyd2'
		call dhyd2()
	end if

c	if (debug) call prd2v(stdo)

	return

	end
