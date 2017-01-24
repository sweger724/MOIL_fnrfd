	subroutine read_add_les(noles,
     1		cpyn,lid,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	integer cpyn,lid
	integer n_cpy(*)
	logical noles

c local
	character*12 name
	integer namel
	integer i,geti,igroup,ipick(maxpt)
	logical find
	integer level,l,l_temp

	name='read_add_les'
	namel = 12
1	continue
c@	debug = .true.
	call rline(name,namel,stdi)
		if (find('*EOD')) then
			if (noles) return
			go to 2
		end if
		noles = .false.
		if (find('MULT')) then
                   lestyp = lestyp + 1
                   lesflag = .true.
                   call pick(ipick,igroup)
                   cpyn = 0
                   cpyn = geti('#cpy',cpyn)
                   if (cpyn .eq. 0) then
                            level = 1
                            write(stdo,*) moname(poimon(i)),' ',
     *                          poimon(i),' ', ptnm(i)
                            call alert(name,namel,
     *                          'Number of copies in les',
     *                          23,level)
                    end if
                    if (cpyn .gt. maxcpy) then
                        level = 1
                        write(stdo,*) moname(poimon(i)),' ',
     *                               poimon(i),' ', ptnm(i)
                        call alert(name,namel,
     *                     'Particle with too many copies',29,level)
                     end if
                     l_temp= geti('lesi',-1)
                     if (l_temp .gt. 100) then
                                level = 1
                                call alert(name,namel,
     *                          'Explicit les id more then 100',
     *                                  28,level)
                        else
                                if (l_temp .eq. -1) then
                                        l = lid
                                        lid = lid + 1
                                    else
                                        l = l_temp
                                end if

                      end if

                      do 100 i = 1,npt
                        if (ipick(i) .eq. 1) then
                                lesid(i) = l
                                n_cpy(i) = cpyn
                        end if
100                    continue
		else
	write(stdo,*) 'read_add_les> illegal input to conn'
	stop
                end if
        go to 1

2       continue
	return
	end
