	subroutine init_add_les(name,namel,igroup,noles,lesflag,
     1		n_cpy,lesid,npt,
     2		lestyp,lid)
	character*7 name
	integer lid,lestyp
	integer namel,igroup,n_cpy(*),lesid(*),npt
	logical noles,lesflag

c local
	integer i,level

	name = 'add_les'
	namel = 7
	igroup = 1
	noles = .true.
	if (lesflag) then
                level = 1
                call alert(name,namel,
     *                     'Old fashion les exists',20,level)
        end if
	lestyp = 0
        do 101 i = 1,npt
                n_cpy(i) = 1
                lesid(i) = 0
101     continue
	lestyp = 0

        lid = 100

	return 
	end
