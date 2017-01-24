	subroutine test_many_ifs(case,j1,j2,j3,j4,lesid,weight,n_cpy)

	integer case,j1,j2,j3,j4,lesid(*),n_cpy(*)

c local
	integer l1,l2,l3,l4
	double precision weight
	logical correct1,correct2


	l1 = lesid(j1)
	l2 = lesid(j2)
	l3 = lesid(j3)
	l4 = lesid(j4)
	
	correct1 = (l1+l2+l3+l4 .eq. 0)
	correct2 = (l1.ne.l2) .and. (l1.ne.l3) .and. (l1.ne.l4) 
     *		.and. (l2.ne.l3) .and. (l2.ne.l4) .and. (l3.ne.l4)
        if (correct1 .or. correct2) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j3)*n_cpy(j4))
c
c case 0 - no les particles OR all particles of different LES types.
c
			case = 0
	else if ((l1.ne.l2) .and. (l2.eq.l3) .and. (l3.eq.l4) ) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j2))
c
c case 11 - two group sof size 3 and size 1 . The different LES particle is first
c
			case = 11
	else if ((l1.ne.l2) .and. (l1.eq.l3) .and. (l1.eq.l4) ) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j2))
c
c case 12 twp groups pf size 3 and size 1. The different LES particle is second
c
			case = 12
	else if ((l1.eq.l2) .and. (l1.ne.l3) .and. (l1.eq.l4) ) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j3))
c
c case 13 twp groups pf size 3 and size 1. The different LES particle is third
c
			case = 13
	else if ((l1.eq.l2) .and. (l1.eq.l3) .and. (l1.ne.l4) ) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j4))
c
c case 14 twp groups pf size 3 and size 1. The different LES particle is fourth
c
			case = 14
	else if ((l1.eq.l2) .and. (l1.ne.l3) .and. (l1.ne.l4) .and.
     1			 (l3.ne.l4)) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j3)*n_cpy(j4))
			case = 212
c
c case 212 two LES particles (1 & 2) with the same lesid, 3 & 4 are different
c
	else if ((l1.eq.l3) .and. (l1.ne.l2) .and. (l1.ne.l4) .and.
     1			(l2.ne.l4)) then
			weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j4))
			case = 213
c
c case 213 two LES particles (1 & 3) with the same lesid, 2 & 4 are different
c
	else if ((l1.eq.l4) .and. (l1.ne.l2) .and. (l1.ne.l3) .and.
     1                  (l2.ne.l3)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j3))
                        case = 214
c
c case 214 two LES particles (1 & 4) with the same lesid, 2 & 3 are different
c
        else if ((l2.eq.l3) .and. (l1.ne.l2) .and. (l2.ne.l4) .and.
     1                  (l1.ne.l4)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j4))
                        case = 223
c
c case 223 two LES particles (2 & 3) with the same lesid, 1 & 4 are different
c
         else if ((l2.eq.l4) .and. (l1.ne.l2) .and. (l1.ne.l3) .and.
     1                  (l3.ne.l4)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j3))
                        case = 224
c
c case 224 two LES particles (2 & 4) with the same lesid, 1 & 3 are different
c
	else if ((l3.eq.l4) .and. (l1.ne.l2) .and. (l2.ne.l4) .and.
     1                  (l1.ne.l4)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j4))
                        case = 234
c
c case 234 two LES particles (3 & 4) with the same lesid, 1 & 2 are different
c
        else if ((l1.eq.l2) .and. (l3.eq.l4) .and. (l1.ne.l3)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j3))
                        case = 2012
c
c case 2012 two groups of the same lesid but different from each other (1 2) & (3 4)
c
        else if ((l1.eq.l3) .and. (l2.eq.l4) .and. (l1.ne.l3)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2))
                        case = 2013
c
c case 2013 two groups of the same lesid but different from each other (1 3) & (2 4)
c
        else if ((l1.eq.l4) .and. (l2.eq.l3) .and. (l1.ne.l2)) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2))
                        case = 2014
c
c case 10: one group all identical
c
	else if ((l1.ne.0) .and. (l1.eq.l2) .and. (l2.eq.l3) .and.
     1		 (l3.eq.l4)) then
			weight=1.d0/n_cpy(j1)
			case = 10
c
c case 2014 two groups of the same lesid but different from each other (1 3) & (2 4)
c
		else
			write(*,*)' Problem cannot match lesid to known arrangement'
			write(*,*)' lesids = ',l1,l2,l3,l4
			write(*,*)' partcl ',j1,j2,j3,j4
			call alert('many_ifs',8,'no match for lesid',18,1)
		end if
		return
		end
