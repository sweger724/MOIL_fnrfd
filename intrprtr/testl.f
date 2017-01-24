	logical function testl()
c
c check if line (line) exists.
c
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
c
c verify that the line exists ...
c
	if (point(100).le.0 .or. nexp.le.0 ) then

		write(stdo,*)' ********************************* '
		write(stdo,*)' Hello there, you forgot the input '
		write(stdo,*)' ********************************* '
		testl =.false.
		return
	else
		testl =.true.
		return
	end if
	end
