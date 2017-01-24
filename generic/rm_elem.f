	subroutine rm_elemi(vec,i,length)
c remove the i-th element of the integer vector "vec" of length - "length"
	integer vec(*),i,length
c local
	integer j

	do 1 j=i,length-1
		vec(j) = vec(j+1)
1	continue

	return
	end

	subroutine rm_elemd(vec,i,length)
c remove the i-th element of the double vector "vec" of length - "length"
	integer i,length
	double precision vec(*)
c local
	integer j

	do 1 j=i,length-1
		vec(j) = vec(j+1)
1	continue

	return
	end
