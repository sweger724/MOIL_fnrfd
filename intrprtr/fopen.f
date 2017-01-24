	logical function fopen(unit)
c check that the unit is indeed opened
	integer unit
	inquire(unit=unit,opened=fopen)
	return
	end
