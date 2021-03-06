        subroutine putdata(ucrd,type,result)

        implicit none
c
c subroutine to read coordinates
c currently very simple (almost no checks) version
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer ucrd
        character*5 type
c local
        character*4 char1,char2,char3,char4
        character*7 name
        integer namel,i,j,level
        double precision xtmp,ytmp,ztmp,moretmp,result(3,*)
        data name/'putdata'/
        data namel/7/

        rewind ucrd
        if (type.eq.'CHARM') then
c
c               do not forget to backspace the file !!
c
        write(ucrd,100)
100     format('* title for CHARMM coordinates')
        write(ucrd,101)
101     format('*')
        write(ucrd,102) npt
102     format(i5)
        do 5 i=1,npt
c write one line of coordinate
          char1 = moname(poimon(i))
          char2 = ptnm(i)
          xtmp = result(1,i)
          ytmp = result(2,i)
          ztmp = result(3,i)
          moretmp = 0
          do 3 j=1,nbulk
           if (i.le.pbulk(j)) then
            char3 = BULK(j)
            go to 4
           end if
3         continue
4         continue
          char4 = 'FREE'
          write(ucrd,103,err=7)i,poimon(i),char1,char2,xtmp
     1          ,ytmp,ztmp,char3,char4,moretmp
103       format(i5,i5,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)   
5        continue
        return
7       continue
        level = 1
        call alert(name,namel,'Error during write',18,level)
        return
        end if
        end
