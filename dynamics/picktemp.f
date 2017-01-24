	subroutine picktemp(tpo,tempi,tempf,tgroup,ntemp,nofreez,inofrz)
c
c Pick temperatures used in different velocity scaling in the
c system. The defualt is that all particles belong to temperature 1.
c Maximum number of temperature is itempg, a variable that is set up
c in dynamics. Currently is 10. It is doubtfull that it will be modified
c tpo - pointer to current temperature tpo(i) is the pointer to tempi & tempf
c	to get the assigned temperature of the i-th atom
c tempi - initial assigned temperature
c tempf - final assigned temperature
c ntemp - actual number of temperatures found
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/MSHAKE.BLOCK'
	integer tpo(*),nofreez(*),tgroup(*),ntemp,inofrz
	double precision tempi(ntemp),tempf(ntemp)
c local
	integer i,j,namel,level
	logical next
	double precision number
	character*8 name


	name = 'picktemp'
	namel = 8
c
c pick should pick the atoms for the different groups
	call pick(tpo,i)
	if (i.eq.0) then
		level = 1
		call alert(name,namel,'No group data',13,level)
	end if

	if (next('tmpi',4)) then
	 do 6 i=1,ntemp
	  tempi(i) = number()
6	 continue
	end if
	if (next('tmpf',4)) then
	 do 7 i=1,ntemp
	  tempf(i) = number()
7	 continue
	 return
	end if
c if this level was reached then there was a problem in reading the
c temperatures
c
	level = 1
	call alert(name,namel,'Value of temperature not found',30,level)
	return
	end
