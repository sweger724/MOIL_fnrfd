	subroutine wsubset(uwcrd,nfrm,ncoor,inofrz,nofreez)
c writing  asubset of dynamics coordinate file according to the
c nofreez vector
c currently the CHARM format is employed for compatability
c with graphic program (QUANTA)
c	
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/COORD.BLOCK'
	integer uwcrd,nfrm,ncoor
	integer inofrz,nofreez(*)

	integer vari(20)
	character*7 name
	character*80 dummy
	integer namel,i

	dummy = 'A dummy title to dynamics coordinates'
	name = 'wsubset'
	namel = 7

	if (ncoor.eq.1) then
		do 1 i=1,20
			vari(i) = 0
1		continue
c
c first writing of coordinates. A title and other junk
c is required
c
c vari(1) is the number of coordinate files to be written
c vari(9) is the new number of atoms
c
		vari(1) = nfrm
		vari(9) = inofrz
		write(uwcrd)'CORD',vari
		write(uwcrd)1,dummy
		write(uwcrd)inofrz
c First time write all coordinates, even if some of them
c are fixed
c
		write(uwcrd)(sngl(coor(1,nofreez(i))),i=1,inofrz)
		write(uwcrd)(sngl(coor(2,nofreez(i))),i=1,inofrz)
		write(uwcrd)(sngl(coor(3,nofreez(i))),i=1,inofrz)
	else
		 write(uwcrd) (sngl(coor(1,nofreez(i))),i=1,inofrz)
		 write(uwcrd) (sngl(coor(2,nofreez(i))),i=1,inofrz)
		 write(uwcrd) (sngl(coor(3,nofreez(i))),i=1,inofrz)
	end if
	return
	end

