	subroutine exclu_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'

	integer poicpy(*),n_cpy(*)
c
c	local
c
	character*13 name
	integer namel
	integer i,j,k,l,ix
	integer i1,i4,i_old,j_old

	name =  'exclu_add_les'
	namel = 13
	ix = 0

	do 10 i=1,npt
	  do 1 j=1,nb
	   if (ib1(j).eq.i) then
	    ix = ix + 1
	    exc2(ix) = ib2(j)
	   end if
1	  continue
	  do 2 j=1,nangl
            if (iangl1(j).eq.i) then
		ix = ix + 1
		exc2(ix) = iangl3(j)
	    end if
2	  continue
	
	 do 3 j=1,ntors
	  i1 = itor1(j)
	  i4 = itor4(j)
	  if (i1.eq.i) then
	    ix = ix + 1
	    exc2(ix) = i4
           else if (i4.eq.i) then
	    ix = ix + 1
	    exc2(ix) = i1
	   end if
3	 continue

	i_old = i
c@
	write(*,*)' exclusion without LES', ix

c
c exclusions due to LES
c
	 if (lesid(i).ne.0) then
	  do 4 j=i+1,npt
	   if (lesid(j).eq.lesid(i) 
     1		.and. cplbl(i).ne.cplbl(j)) then
			ix = ix + 1
			exc2(ix) = j
	   end if
4	  continue
	 
	 end if
	 

	exc1(i) = ix

10	continue
	write(*,*) ' exclsuion with LES ',ix

c
c compress duplicates in the exclusion list. These duplication occurs because of multiple
c appearance of 1-4 pairs in the torsions list. Either because of AMBER torsions or becuase
c of ring torsions/angles
c
	do 20 i=1,npt
	 if (exc1(i-1)+1.ne.exc1(i)) then
	  j = exc1(i-1)
11	  j = j + 1
	  k = j
12	  k = k + 1
          if (k.le.exc1(i)) then
		if (exc2(j).eq.exc2(k)) then
		 call rm_elemi(exc2,k,ix)
		 ix = ix -1
		 do 13 l=i,npt
			exc1(l) = exc1(l) -1
13		 continue
		 j = j -1
		 go to 11
		end if
		go to 12
	  end if
	  if (j.lt.exc1(i)-1) go to 11
	 end if
20	continue
		 
	totex = ix
c 
c at this point there should be no repeats in the exclusion list
c so we do not add a check to that effect that is present in the old code.
c
	if (totex.gt.maxex) then
		write(stdo,*)' intermediate totex maxex ',totex,maxex
		call alert(name,namel,'Max excl exceeded',17,1)
	end if

	return
	end
