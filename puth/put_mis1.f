       subroutine put_mis1()
c
c Identifying missing hydrogens and placing them in
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        character*8 name
        integer namel
        integer i,cycles

        integer level
        logical success,progress

        name   = 'put_mis1'
        namel  = 8
	cycles = 10
	do 10 j=1,cycles
	progress = .false.
c
c check for unidentifed coordinates
c
        do 1 i=1,npt
         if (coor(1,i).gt.9998.d0) then
	  write(*,*)' in put_mis1 '
          write(stdo,100)i,ptnm(i)
100       format(1x,' Coordinates of particle ',i6,' name ',a4,
     1' are not defined')
          if (ptnm(i)(1:1).eq.'H' .or. ptnm(i)(1:2).eq.'1H'
     1.or. ptnm(i)(1:2).eq.'2H') then
                call pos(success,i)
		if (success) progress = .true.
                if (.not.success) then
                 level = 1
                 call alert(name,namel,'Cannot position H',17,level)
                end if
          else
            call put_mis2(success,i)
	    if (success) progress = .true.
            if (.not.success) then
             call alert(name,namel,'Missing crd for none H',22,level)
            else
             write(*,*)' atom ',i,' ',ptnm(i),' was built '
            end if
          end if
         end if
1       continue
	if (.not.progress) return
10	continue
        end
