	subroutine wdynvel(uwvel,nfrm,nvelo,inofrz,nofreez)
c writing dwn dynamics coordinate file
c currently the CHARM format is employed for compatability
c with graphic program (QUANTA)
c	
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	integer uwvel,nfrm,nvelo
	integer inofrz,nofreez(*)

	integer vari(20)
	character*7 name
	character*80 dummy
	integer namel,i

	dummy = 'A dummy title to dynamics coordinates'
	name = 'wdyncrd'
	namel = 7

	if (nvelo.eq.1) then
		do 1 i=1,20
			vari(i) = 0
1		continue
c
c first writing of coordinates. A title and other junk
c is required
c
c vari(1) is the number of coordinate files to be read
c vari(9) is the number of FIXED atoms
c
		vari(1) = nfrm
		vari(9) = npt - inofrz
		write(uwvel)'VELO',vari
		write(uwvel)1,dummy
		write(uwvel)npt
		if (inofrz.ne.npt) then
		 write(uwvel)(nofreez(i),i=1,inofrz)
		end if
c First time write all coordinates, even if some of them
c are fixed
c
		write(uwvel)(sngl(velo(1,i)),i=1,npt)
		write(uwvel)(sngl(velo(2,i)),i=1,npt)
		write(uwvel)(sngl(velo(3,i)),i=1,npt)
	else
		if (inofrz.ne.npt) then
		 write(uwvel) (sngl(velo(1,nofreez(i))),i=1,inofrz)
		 write(uwvel) (sngl(velo(2,nofreez(i))),i=1,inofrz)
		 write(uwvel) (sngl(velo(3,nofreez(i))),i=1,inofrz)
		else
		 write(uwvel) (sngl(velo(1,i)),i=1,npt)
		 write(uwvel) (sngl(velo(2,i)),i=1,npt)
		 write(uwvel) (sngl(velo(3,i)),i=1,npt)
		end if
	end if
	return
	end
