	subroutine prtc_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:maxpt),n_cpy(maxpt)


c local
	integer i,j
	double precision weight

c Create temporary particle list.
        tnpt = 1
	tnchg = 0
        poicpy(0) = 0
        do 300 i = 1,npt
                weight = 1.d0/n_cpy(i)
                do 310 j = 1,n_cpy(i)
                        tpoimon(tnpt)   = poimon(i)
                        if (mdivyes) then
                                tpoidmon(tnpt) = poidmon(i)
                        end if
                        tptid(tnpt)     = ptid(i)
                        tptnm(tnpt)     = ptnm(i)
                        tptms(tnpt)     = ptms(i) * weight
                        tptchg(tnpt)    = ptchg(i) * weight
                        if (.not.arith) then
                         tepsgm6(tnpt)  = epsgm6(i) * weight
                        else
                         tepsgm6(tnpt)  = epsgm6(i)
                        end if
                        tepsgm12(tnpt)  = epsgm12(i) * weight
                        tptwei(tnpt)    = weight
                        tlesid(tnpt)    = lesid(i)
                        tflagchr(tnpt)  = flagchr(i)
			if (tflagchr(tnpt)) then
				tnchg = tnchg + 1
			end if
                        if (lesid(i) .eq. 0) then
                                cplbl(tnpt) = 0
                            else
                                cplbl(tnpt) = j
                        end if
                        tnpt = tnpt + 1
310             continue
                poicpy(i) = tnpt  - 1
300     continue
        tnpt = tnpt - 1
	return
	end
