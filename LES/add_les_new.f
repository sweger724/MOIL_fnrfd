	subroutine add_les()
c Subroutine for multiplying les particles.

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/PROPERT.BLOCK'
	include 'COMMON/MONOMERS.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	character*7 name
	integer namel
	integer n_cpy(maxpt),poicpy(0:maxpt)
	integer lid,cpyn
	logical noles
	integer igroup

	name = 'add_les'
	namel=7

	call init_add_les(name,namel,igroup,noles,
     1		lesflag,n_cpy,lesid,npt,lestyp,lid)
	call read_add_les(noles,cpyn,lid,n_cpy)
	if (noles) return
	call mdiv_add_les(totmon,mdivlist,n_cpy)
	call prtc_add_les(poicpy,n_cpy)
	call bond_add_les(poicpy,n_cpy)
	call list14_add_les(poicpy,n_cpy)
	call angle_add_les(poicpy,n_cpy)
	call tors_add_les(poicpy,n_cpy)
	call impr_add_les(poicpy,n_cpy)
	call morse_add_les(poicpy,n_cpy)
	call molc_add_les(poicpy)

c
	return
	end

