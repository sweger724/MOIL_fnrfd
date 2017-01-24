        subroutine convertpath(urpath,uwpath,nstru,tobin)
c
c subroutine to read PATH format coordinate files
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer urpath,uwpath,nstru,tobin

c local
        character*11 name
        integer namel,level,i,j,k
        real e0

        name  = 'convertpath'
        namel = 12

        rewind urpath
        rewind uwpath


c        do 1 k=1,nstru
c         read(upath,err=999)e0,((coor(j,i),i=1,npt),j=1,3)
c1       continue

100     format(1x,' *** CONVERTING PATH FILE (binary to formatted)',
     &        i5,' ENERGY = ',e8.3)
101         format(1x,' *** CONVERTING PATH FILE (formatted to binary)',
     &        i5,' ENERGY = ',e8.3)

        if (tobin.eq.0) then ! convert to formatted
          do k=1,nstru
            read(urpath,err=998)   e0,((coor(j,i),i=1,npt),j=1,3)

            do i=1,npt
              write(uwpath,*,err=999)coor(1,i),coor(2,i),coor(3,i)
            enddo

            write(stdo,100)k,e0

          enddo

        else ! convert formatted to binary
          do k=1,nstru
            do i=1,npt
              write(urpath,*,err=999) coor(1,i),coor(2,i),coor(3,i)
            enddo

            write(uwpath,err=999) e0,((coor(j,i),i=1,npt),j=1,3)
            write(stdo,101)k,e0
          enddo
        endif
        return

998     continue
        level = 1
        call alert(name,namel,'Error while reading Path file',29,level)
999     continue
        level = 1
        call alert(name,namel,'Error while writing Path file',29,level)
        return
        end
