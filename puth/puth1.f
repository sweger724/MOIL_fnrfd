       subroutine puth1()
c
c Identifying missing hydrogens and placing them in
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        character*5 name
        integer namel
        integer i

        integer level
        logical success

        name   = 'puth1'
        namel  = 5
c
c check for unidentifed coordinates
c
        do 1 i=1,npt
         if (coor(1,i).gt.9998.d0) then
          write(stdo,100)i,ptnm(i)
100       format(1x,' Coordinates of particle ',i6,' name ',a4,
     1' are not defined')
          if (ptnm(i)(1:1).eq.'H' .or. ptnm(i)(1:2).eq.'1H'
     1.or. ptnm(i)(1:2).eq.'2H') then
                call pos(success,i)
                if (.not.success) then
                 level = 3
                 call alert(name,namel,'Cannot position H',17,level)
                end if
           else
	    call put_mis2(success,i)
	    if (.not.success) then
             level = 3
             call alert(name,namel,'Unable to build non H',21,level)
            end if
          end if
         end if
1       continue
        return
        end
