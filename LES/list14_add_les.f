	subroutine list14_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'

	integer poicpy(0:maxpt),n_cpy(maxpt)
c
c local
	character*14 name
	integer i,j,k,l,m,n,old_totspe
	integer jbeg,kbeg
	integer namel
	double precision weight

	name = 'list14_add_les'
	namel = 14
	old_totspe = totspe

c
c create first an uncompressed list of 1-4 interactions
c
	  totspe = 1
	  do 1 i=1,old_totspe
	    j = spec1(i)
	    k = spec2(i)
	    jbeg = poicpy(j-1)+1
	    kbeg = poicpy(k-1)+1
	    if ( ((lesid(j).eq.0) .and. (lesid(k).eq.0)) .or.
     1		(lesid(j).ne.lesid(k)) ) then
		weight = 1.d0/(n_cpy(j)*n_cpy(k))
		do 2 l=jbeg,poicpy(j)
		 do 2 m=kbeg,poicpy(k)
			tspec1(totspe) = l
			tspec2(totspe) = m
			do 11 n=1,3
			tp14(n,totspe)=p14(n,totspe)*weight
11			continue
			totspe         = totspe + 1
2		continue
	     else
	    	l = jbeg
		weight = 1.d0/n_cpy(j)
		do 3 m=kbeg,poicpy(k)
			tspec1(totspe) = l
			tspec2(totspe) = m
			do 21 n=1,3 
			tp14(n,totspe)=p14(n,totspe)*weight
21			continue
			totspe = totspe + 1
			l = l + 1
3		continue
	      end if
1	   continue
	   totspe = totspe - 1

	   do 4 i=1,totspe
		spec1(i)=tspec1(i)
		spec2(i)=tspec2(i)
		do 4 j=1,3
		p14(j,i)=tp14(j,i)
4	   continue
	   return
	   end
