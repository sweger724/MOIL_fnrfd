	subroutine fill_angle(k1,k2,k3,i,weight)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer i,k1,k2,k3,k4
	double precision weight

	tnangl = tnangl + 1

	tiangl1(tnangl)  = k1
	tiangl2(tnangl)  = k2
	tiangl3(tnangl)  = k3
	tkangl(tnangl)   = weight*kangl(i)
	tangleq(tnangl) = angleq(i)

	return
	end
