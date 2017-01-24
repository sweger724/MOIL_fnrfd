        subroutine wpath(upath)
c
c subroutine to write PATH format coordinate files
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer upath

c local
        character*5 name
        integer namel,i,j,level
        double precision e0

        name  = 'wpath'
        namel = 5
        e0    = e_total
        write(upath,err=999)e0,((coor(j,i),i=1,npt),j=1,3)
         write(stdo,100)e0
100      format(1x,' *** WRITING PATH FILE,  ENERGY = ',f15.5)
        return
999     continue
        level = 1
        call alert(name,namel,'Error while writing Path file',29,level)
        return
        end

c------------------------------------------------------------------
        subroutine wpathform(upath)
c
c subroutine to write PATH format coordinate files
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer upath

c local
        character*5 name
        integer namel,i,j,level
        double precision e0

        name  = 'wpath'
        namel = 5
        e0    = e_total
 
        write(stdo,*)' *** WRITING FORMATTED PATH FILE ***'

        write(upath,*,err=999)(coor(1,i),coor(2,i),coor(3,i),i=1,npt)

        return
999     continue
        level = 1
        call alert(name,namel,'Error while writing Path file',29,level)
        return
        end


c------------------------------------------------------------------
        subroutine wxyzpath(upath)
c
c subroutine to write PATH format coordinate files
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer upath

c local
        character*5 name
        integer namel,i,j,level
        double precision e0

100     format(1x,' *** WRITING XYZ PATH FILE,  ENERGY = ',f15.5)

        name  = 'wpath'
        namel = 5
        e0    = e_total
        do i=1,npt
          write(upath,*)coor(1,i),coor(2,i),coor(3,i)
        enddo
        write(stdo,100)e0
        return

999     continue
        level = 1
        call alert(name,namel,'Error while writing Path file',29,level)
        return
        end
