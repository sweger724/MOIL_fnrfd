	subroutine molc_add_les(poicpy)

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:*)
	integer i,j,ix

	write(*,*)' BEFORE poipt ',(poipt(i),i=1,totmon)
	do 1 i=1,npt
	 texc1(i) = exc1(i)
1	continue
	do 2 i=1,totex
	 texc2(i) = exc2(i)
2	continue

	ix = 0
	do 2 i=1,npt
	 do 2 j=texc1(i-1)+1,texcc1(i))
	  k = texc2(i)
	 do 1 j = poicpy(k-1)+1,poicpy(k)
	  ix = ix + 1
	  exc2(ix) = j 
1	 continue
	
	do 1 i=1,totmon
		poipt(i) = poicpy(poipt(i))
1	continue
	write(*,*)' poipt ',(poipt(i),i=1,totmon)

	npt = tnpt


	do 2 i=1,npt
		poimon(i)   = tpoimon(i)
		if (mdivyes) then
		 poidmon(i) = tpoidmon(i)
		end if
                ptid(i)     = tptid(i)
                ptnm(i)     = tptnm(i)
                ptms(i)     = tptms(i)
                ptchg(i)    = tptchg(i)
                epsgm6(i)   = tepsgm6(i)
                epsgm12(i)  = tepsgm12(i)
                ptwei(i)    = tptwei(i)
                lesid(i)    = tlesid(i)
                flagchr(i)  = tflagchr(i)
2       continue

	pbulk(1) = poipt(totmon)
	return
	end

