	subroutine fill_improper(k1,k2,k3,k4,i,weight)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer k1,k2,k3,k4,i
	double precision weight

	tnimp = tnimp + 1

	tiimp1(tnimp) = k1
	tiimp2(tnimp) = k2
	tiimp3(tnimp) = k3
	tiimp4(tnimp) = k4
	tkimp(tnimp)  = kimp(i)*weight
	timpeq(tnimp) = impeq(i)

	return
	end
