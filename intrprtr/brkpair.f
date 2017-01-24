	subroutine brkpair(exp1,l1,exp2,l2,special)
c
c find if the current line include an expression separated by "-"
c this corresponds (for example to bond between two atoms)
c the two expression separated by "-" will be assigned to exp1
c and exp2 with corresponding length l1 and l2
c On retrun the "-" that was read will be modified to "+"
c
	character*1 dash
	character*4 exp1,exp2
	integer l1,l2,istart
	logical special
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
c
c local
c
	integer i,j,k1,k2,level
	data dash/'-'/

	exp1 = ' '
	exp2 = ' '
	special = .false.

	if (point(100).le.0 .or. nexp.le.0 ) then
		level = 0
		call alert('brkpair',7,'No input line',13,level)
		return
	end if
	k1=1
	k2=point(1)
	do 3 i=1,nexp
	 do 2 j=k1,k2
	  if (line(j:j).eq.dash) then
	   if (i.eq.1) then
	    istart = 1
	   else
	    istart = point(i-1) + 1
	   end if
	   l1 = j - istart
	   l2 =point(i) - j
	   if (l1.gt.4 .or. l2.gt.5) then
		write(stdo,100)line(istart:point(i))
100		format(1x,' problem in expression',a)
		call alert('brkpair',7,' Atom names too long ',21,1)
	   end if
	   if (line(point(i):point(i)).eq.'*') then
	    special = .true.
	    l2 = l2 - 1
	   end if
	   if (l2.eq.5 .and. (.not. special)) then
	     write(stdo,100)line(istart:point(i))
	     call alert('brkpair',7,' Atom names too long ',21,1)
	   end if
	   if (istart.eq.1 .and. j.eq.1) then
	    call alert('brkpair',7,' Unreasonable test ',19,1)
	   else
	   exp1=line(istart:j-1)
	   end if
	   if (.not. special) then
	    exp2=line(j+1:point(i))
	   else
	    exp2=line(j+1:point(i)-1)
	   end if
	   line(j:j) = '+'
	   return
          end if
2	 continue
	 k1=point(i)
	 k2=point(i+1)
3	continue
	l1 = -1
	l2 = -1
	exp1 = ' '
	exp2 = ' '
	return
	end
