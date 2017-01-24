	subroutine copy_angles()
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer i

	nangl = tnangl

	do 1 i=1,nangl

		iangl1(i) = tiangl1(i)
		iangl2(i) = tiangl2(i)
		iangl3(i) = tiangl3(i)
		kangl(i)  = tkangl(i)
		angleq(i) = tangleq(i)

1	continue

	return
	end
