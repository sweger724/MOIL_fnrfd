	subroutine copy_improper()
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer i

	nimp = tnimp

	do 1 i=1,nimp

		iimp1(i)=tiimp1(i)
		iimp2(i)=tiimp2(i)
		iimp3(i)=tiimp3(i)
		iimp4(i)=tiimp4(i)
		kimp(i) =tkimp(i)
		impeq(i)=timpeq(i)

1	continue

	return
	end
