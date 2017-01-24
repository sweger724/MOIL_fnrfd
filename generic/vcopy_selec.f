
	subroutine vdcopy(v1,v2,pointr,nselec)
c
c This routine copy v1 (a double precision vector length - length)
c to v2
c
	double precision v1(*),v2(*)
	integer pointr(*),length
c
c local
	integer i

	do 1 i=1,length
	 v2(i) = v1(i)
1	continue

	return
	end


	subroutine vicopy(v1,v2,length)
c
c This routine copy v1 (an integer vector length - length) to v2
c
	integer v1(*),v2(*)
	integer length
c
c local
	integer i

	do 1 i=1,length
	 v2(i) = v1(i)
1	continue

	return
	end

	subroutine vccopy(v1,v2,length,size)
c
c This routine copy v1 (a character*size vector of length - "length")
c to v2
c
	character*4 v1(*),v2(*)
	integer length,size
c
c local
	integer i

	do 1 i=1,length
	 v2(i) = v1(i)
1	continue

	return
	end
