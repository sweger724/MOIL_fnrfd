        subroutine rpath(upath,nstru)
c
c subroutine to read PATH format coordinate files
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer upath,nstru

c local
        character*5 name
        integer namel,level,i,j,k
        double precision e0

        name  = 'rpath'
        namel = 5

                 do i=1,npt
                 chain(i) = '*'
                 enddo

        rewind upath

        do 1 k=1,nstru
         read(upath,err=999)e0,((coor(j,i),i=1,npt),j=1,3)
1       continue

        write(stdo,100)nstru,e0
100     format(1x,' *** READING PATH FILE ',i5,' ENERGY = ',e8.3)
        return
999     continue
        level = 1
        call alert(name,namel,'Error while reading Path file',29,level)
        return
        end
