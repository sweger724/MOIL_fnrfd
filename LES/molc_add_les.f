	subroutine molc_add_les(poicpy)

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:*)
	integer i,j,k,l,m,n,ix

	do 1 i=1,npt
	 texc1(i) = exc1(i)
1	continue
	do 2 i=1,totex
	 texc2(i) = exc2(i)
2	continue

	ix = 0
	do 5 i=1,npt
c n is the number of LES copies for no-LES n is eq 1
c
	 n = poicpy(i)-poicpy(i-1)
c run a loop over the LES copies of a "real" particle (before multiplication) i
	 do 5 m=poicpy(i-1)+1,poicpy(i)

c if LES add self copies to the exclusion list
c
	 if (n.gt.1) then
	 do 4 l=m+1,poicpy(i)
	   ix = ix + 1
	   exc2(ix) = l
4	 continue
	 end if
c
c expand no-LES exclusion list to account for extra LES particles
c
	 do 3 j=texc1(i-1)+1,texc1(i)
	  k = texc2(j)
	 do 3 l = poicpy(k-1)+1,poicpy(k)
	  ix = ix + 1
	  exc2(ix) = l 
3	 continue
		exc1(m) = ix
5	continue
	
	do 6 i=1,totmon
		poipt(i) = poicpy(poipt(i))
6	continue
	totex = ix
	npt = tnpt
	nchg = tnchg


	do 7 i=1,npt
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
7       continue

	pbulk(1) = poipt(totmon)
	return
	end

