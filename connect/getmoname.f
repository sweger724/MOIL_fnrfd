	subroutine getmoname(name,imult,less,empty,copies
     1	,namelst,namevec,index,itype,numon,lestype,weight)
c
c get monomer names, check if it is of LES type, feed the
c monomers name array and identify the monomer type
c name - the name to be picked from line. If LOOS on output
c	 no monomer name was found
c imult - if LES particle then name = mult and imult increased by
c	  one and corresponds to the number of LES types.
c less   - a logical flag if set to true getmono identifed LES monomer
c empty  - a logical flag. if true no monomer was found
c copies - number of LES copies that is used
c namelst - list of possible monomer names
c namevec - the vector of monomer names in the macromolecule
c index   - the number of monomer in the molecular monomer list
c		is updated in the present routine.
c itype   - index of monomer type
c numon   - number of different monomers typr
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	character*4 name,namelst(*),namevec(*)
	integer imult,copies,index,lestype
	integer itype,numon,level
	double precision weight(*)
	double precision number
	logical next
	logical less,empty
c local
	character*9 subname
	integer i,namel
	integer intgr

	subname = 'getmoname'
	less  = .false.
	empty = .false.
	lestype = -1
	weight(1) = -1.d0

	call get4c(name,empty)
	if (empty) then
	 name = 'LOOS'
	 return
	end if
	if (name.eq.'mult') then
	 imult = imult + 1
	 less = .true.
	 call get4c(name,empty)
	 if (empty) then
	  level = 1
	  call alert(subname,namel,'mono mult. but no names',23,level)
	 end if
	 copies = intgr()
	 do 1 i=index+1,index+copies
	  namevec(i)=name
1	 continue 
	 index = index + copies
	 if (next('ltyp',4)) lestype = intgr()
	 if (next('weig',4)) then
	  do 15 i=1,copies
	   weight(i) = number()
15	  continue
	 end if
	else
	 index = index + 1
	 namevec(index)=name
	end if
	do 2 i=1,numon
		if (debug) then
		 write(stdo,*)' name namelst(i) ',name,' ',namelst(i)
		end if
		if (name.eq.namelst(i)) then
			itype = i
			return
		end if
2	continue
	level = 1
        write(*,*)' Unfound monomer is ',name(1:4)
	call alert(subname,namel,'Mono type not found',20,level)
	return
	end
