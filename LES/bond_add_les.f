	subroutine bond_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:maxpt),n_cpy(maxpt)

c local
	integer i,j,k,kbeg,jbeg
	double precision weight
	integer m,n

c Creat tempotrary bond list
        tnb = 1
        do 400 i = 1,nb
                m = ib1(i)
                n = ib2(i)
                jbeg = poicpy(m - 1) + 1
                kbeg = poicpy(n - 1) + 1
                if (((lesid(m) .eq. 0) .and. (lesid(n) .eq. 0)) .or.
     *                (lesid(m) .ne. lesid(n))) then
                        weight = 1.d0/n_cpy(n) * 1.d0/n_cpy(m)
                        do 410 j = jbeg, poicpy(m)
                           do 420 k = kbeg, poicpy(n)
                                tib1(tnb)  = j
                                tib2(tnb)  = k
                                tkbond(tnb) = kbond(i) * weight
                                treq(tnb)   = req(i)
                                tnb = tnb + 1
420                        continue
410                     continue
                   else
                        k = kbeg
                        weight = 1.d0/n_cpy(n)
                        do 430 j = jbeg, poicpy(m)
                                tib1(tnb)  = j
                                tib2(tnb)  = k
                                tkbond(tnb) = kbond(i) * weight
                                treq(tnb)   = req(i)
                                tnb = tnb + 1
                                k = k + 1
430                     continue
                end if
400     continue
        nb = tnb - 1

	do 440 i=1,nb
		ib1(i)   = tib1(i)
		ib2(i)   = tib2(i)
		kbond(i) = tkbond(i)
		req(i)   = treq(i)
440	continue

	return
	end

