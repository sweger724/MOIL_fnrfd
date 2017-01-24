	subroutine ralign(ualign)
c
c Read alignment data, the one stored in
c dyna_mpi common block: DYNA_MPI.BLCOK
c
	include	'COMMON/LENGTH.BLOCK'
	include 'COMMON/COUPLED_DYNA.BLOCK'

	integer ualign,namel

	character*5 name
	character*1 TILDA
	character*10 fff

	data name/'ralign'/
	data namel/6/

	rewind ualign

	read(ualign,100)TILDA
100	format(a1)
c Read length of the allignment alignN
	read(ualign,101)alignN
101	format(i6)

c@@@@@@@@
	write(*,*) 'alignN: ',alignN
c@@@@@@@@@

	read(ualign,100)TILDA
	read(ualign,100)TILDA
	
	write(fff,'(i3,a2)'), alignN, 'i4'
	write(*,*)'fff: ', fff
	read(ualign,fff)(align(i),i=1,alignN)
c103	format(fff)

	do 110 j=1,alignN
	    write(*,*)'align(',j,')= ',align(j)
110	continue

	return
	end
