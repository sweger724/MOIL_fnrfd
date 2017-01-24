        subroutine getdata(uvel,type,result)

        implicit none
c
c subroutine to read coordinates
c currently very simple (almost no checks) version
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer uvel
        character*5 type
c local
        character*4 char1,char2,char3,char4
        character*7 name
        character*80 zevel
        integer leng,namel,i,j,readmo,readat
        integer istart,i1,j1,k,level
        double precision xtmp,ytmp,ztmp
        double precision moretmp, result(3,*)
        data name/'getdata'/
        data namel/7/

c initialize all data to 9999.0
        do i=1,npt
         result(1,i) = 9999.d0
         result(2,i) = 9999.d0
         result(3,i) = 9999.d0
        end do

        
1       continue
        rewind uvel
        if (type.eq.'CHARM') then
c
c               read the title and find out its length
c
        leng = 0
2       continue
        read(uvel,105,err=7,end=8)zevel
105     format(80a)
        if (zevel(1:1).eq.'*') then
         do 21 k=80,1,-1
           if (zevel(k:k).ne.' ') then
            call echo(name,namel,zevel,k)
            go to 22
           end if
21       continue
22       continue
         leng = leng + 1
         go to 2
        end if
c
c               check that title exist
c
        if (leng.eq.0) then
         level = 1
         call alert(name,namel,'vel file with no title',22,level)
         return
        end if
c
c               do not forget to backspace the file !!
c
        backspace uvel
        read(uvel,*,err=7,end=8) i
        if (i.ne.0 .and. i.ne.npt) then
         level = 1
         call alert(name,namel,'Number of particles do not
     1 match',32,level)
         return
        end if
        readat = 0
        readmo = 0
        istart = 1
        do 5 i=1,totmon
         do 4  j=istart,poipt(i)
c read one line of coordinate
          read(uvel,100,err=7,end=8)i1,j1,char1,char2,xtmp
     1          ,ytmp,ztmp,char3,char4,moretmp
100       format(i5,i5,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)   
c check that this is the right monomers name
          if (char1.ne.moname(i)) then
c may be we just miss some atoms (this is ok sometime..
c for example when hydrogens are not included
           if (i1-i.eq.1) then
            level = 0
            write(stdo,101)i,moname(i),j,poipt(i)-istart+1
101         format(1x,' Missing pts in monomer ',i5,' type =',
     1          a4,'  read ',i5,' need ',i5)
            call alert(name,namel,'Missing particles',17,level)
            backspace uvel
            go to 5
           else if (i1-i.eq.-1) then
            level = 1
            write(stdo,102)i1,moname(i1),i,moname(i)
102         format(1x,'Too many particles in monomer ',i5,
     1          ' name = ',a4,' expect ',i5,' name -',a4)
            call alert(name,namel,'Too many particles',18,level)
           else
            level = 1
            write(stdo,103)char1,moname(i)
103         format(1x,'Read = ',a4,' Expected = ',a4)
            call alert(name,namel,'Monomers do not match',21,level)
           end if
          else
            do 3 k=istart,poipt(i)
             if (char2.eq.ptnm(k)) then
              result(1,k) = xtmp
              result(2,k) = ytmp
              result(3,k) = ztmp
              more(k) = moretmp
             end if
3           continue
          end if
4         continue
          istart=poipt(i) + 1
5        continue
c check that there are not any unidentified atoms
         do 6 i=1,npt
          if (result(1,i).gt.9998.0) then
           write(stdo,104)i,ptnm(i)
104        format(1x,' Velocities of particle ',i5,2x,a4,
     1          ' not defined')
           level = 0
           call alert(name,namel,'Undefined coordinates',21,level)
          end if
6        continue
         return
7        continue
         level = 0
         call alert(name,namel,'Error during read',17,level)
         goto 1
         return
8        continue
         level = 0
         call alert(name,namel,'End of file encounter',21,level)
         goto 1
         return
        end if
        end
