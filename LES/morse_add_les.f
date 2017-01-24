	subroutine morse_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:*),n_cpy(*)
c local
	integer i,j,k,n,m,jbeg,kbeg
	double precision weight

	tnmb = 1
        do 450 i = 1,nmb
                m = imb1(i)
                n = imb2(i)
                if (m .eq. 1) then
                        jbeg = 1
                   else
                        jbeg = poicpy(m - 1) + 1
                end if
                if (n .eq. 1) then
                        kbeg = 1
                   else
                        kbeg = poicpy(n - 1) + 1
                end if
                if (((lesid(m) .eq. 0) .and. (lesid(n) .eq. 0)) .or.
     *               (lesid(m) .ne. lesid(n))) then
                        weight = 1.d0/n_cpy(n) * 1.d0/n_cpy(m)
                        do 460 j = jbeg, poicpy(m)
                           do 470 k = kbeg, poicpy(n)
                                timb1(tnb)  = j
                                timb2(tnb)  = k
                                trmeq(tnb)  = rmeq(i)
                                tD(tnb)     = D(i) * weight
                                talpha(tnb) = alpha(i)
                                tbeta1(tnb) = beta1(i)
                                tnmb = tnmb + 1
470                        continue
460                     continue
                   else
                        k = kbeg
                        weight = 1.d0/n_cpy(n)
                        do 480 j = jbeg, poicpy(m)
                                timb1(tnb)  = j
                                timb2(tnb)  = k
                                trmeq(tnb)  = rmeq(i)
                                tD(tnb)     = D(i) * weight
                                talpha(tnb) = alpha(i)
                                tbeta1(tnb) = beta1(i)
                                tnmb = tnmb + 1
                                k = k + 1
480                     continue
                end if
450     continue
        nmb = tnmb - 1

	do 490 i=1,nmb
		imb1(i) = timb1(i)
		imb2(i) = timb2(i)
		rmeq(i) = trmeq(i)
		d(i)    = td(i)
		alpha(i)= talpha(i)
		beta1(i) = tbeta1(i)
490	continue

	return
	end

