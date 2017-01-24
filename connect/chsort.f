	subroutine chsort()
c
c a subroutine to find number of charges particles and then to
c set logical flags for chanrged/uncharged particles
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'

c local
	integer i
        double precision chargetot
c find number of charged particles and set logical flags for charges

        open(155,file='total-charge.dat',status='replace')
        chargetot = 0.0
	do 1 i=1,npt
		if (dabs(ptchg(i)).gt.1.d-6) then
			nchg        = nchg + 1
                        chargetot = chargetot + ptchg(i)
			flagchr(i)  = .true.
		else
			flagchr(i)  = .false.
		end if
1	continue
	
	write(stdo,100)nchg
100	format('chsort> Number of charged particled',1x,i10)
        write(155,*) 'total charge of system =',chargetot

	return
	end
