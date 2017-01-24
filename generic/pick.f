        subroutine pick(ipick,ngrp)
c
c pick particles. On return ipick has integer values for all
c the atoms which were selected. If the group keyword
c was used ipick will have integer values from 1 up according
c to the group found. If only one group exists ipick is set to
c 1 for picked atoms. If no selection is made all particles are
c picked. It is assumed that rline was called BEFORE pick
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer ipick(maxpt),ipick1(maxpt)
        integer ngrp
        logical find,next

        integer i,j
        integer namel,level
        logical first
        character*4 name

        name='pick'
        namel=4

c
c initialize ipick to default value
        do 1 i=1,npt
         ipick(i)  = 1
         ipick1(i) = 0
1       continue
        ngrp = 0
        first = .true.
        if (.not.find('pick'))  then
                call alert(name,namel,'Did not find pick',17,0)
        end if

c start loop on pick-up line
c Go from left to right and pick up logical expressions
2       continue
        if (next('done',4)) go to 6

c .and.(&) or first
        if (next('&   ',1) .or.first) then
         first = .false.
         call pick2(ipick1,ngrp)
         do 3 i=1,npt
          if(ipick(i).gt.0 .and. ipick1(i).gt.0) then
           ipick(i)=ipick1(i)
          else
           ipick(i)=0
          end if
3        continue
         go to 2

c .or.(|)
        else if (next('|   ',1)) then
         call pick2(ipick1,ngrp)
         do 4 i=1,npt
          if (ipick1(i).gt.0)  ipick(i) = ipick1(i)
4        continue
         go to 2

c .not. (!=)
        else if (next('!=  ',2)) then
         call pick2(ipick1,ngrp)
         do 5 i=1,npt
          if (ipick1(i).gt.0) ipick(i) = 0
5        continue
         go to 2
        end if
        level = 1
        write(stdo,100)
100     format(1x,'Illegal logical pick syntax')
        call alert(name,namel,' Missing done',13,level)
6       continue
        j = 0
        do 7 i=1,npt
         if (ipick(i).ne.0) j = j + 1
7       continue
        write(stdo,101)j,npt
101     format(/,1x,'pick> ',i6,' particles picked out of ',i6,/)
        if (debug) then
         write(stdo,*)' Leaving... pick = ',(ipick(i),i=1,npt)
        end if
        return
        end
c----------------------------------------------------------


c-----------------------------------------------------------
        subroutine pick2(ipick,ngrp)
c
c pick atoms in a logical segment. I.e. expression enclosed
c between & (.and.) or | (.or.) signs.
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer ipick(*)
        integer ngrp
        integer intgr
        logical next

        character*5 name
        integer namel,level

        integer i,j,k,i1,i2,i3
        logical empty
        character*4 char1

        name='pick2'
        namel=5

        if (debug) write(stdo,*)' pick2 was called '
c initialize ipick
        do 1 i=1,npt
          ipick(i) = 0
1       continue

c check if groups exist
        if (next('grou',4)) then
          ngrp = ngrp + 1
          j = intgr()
        else
          ngrp = 1
          j    = 1
        end if

c start processing
        if (next('#prt',4)) then
          i1 = intgr()
          i2 = intgr()
          do 3 i=i1,i2
            ipick(i) = j
3         continue
        else if (next('#mon',4)) then
          i1 = intgr()
          i2 = intgr()
          if (i1.eq.1) then
           i3=1
          else
           i3 = poipt(i1-1)+1
          end if
          do 4 i=i3,poipt(i2)
            ipick(i) = j
4         continue
c chen
        else if (next('#les',4)) then
          i1 = intgr()
          do 77 i=1,npt
            if (i1 .eq. cplbl(i)) ipick(i) = j
77        continue
        else if (next('chem',4)) then
          call get4c(char1,empty)
          if (debug) write(stdo,*)' char1 = ',char1(1:4)
          if (empty) then
           level = 1
           call alert(name,namel,'Illegal chem pickup',19,level)
          end if
          if (char1.eq.'prtc') then
           call get4c(char1,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'missing particle name',21,level)
           end if
           do 5 i=1,npt
c
c check name match as well as possible wildcards
c
           if (ptnm(i).eq.char1 .or. char1(1:1).eq.'*') ipick(i) = j
           if (ptnm(i)(1:1).eq.char1(1:1).and.char1(2:2).eq.'*')
     &                ipick(i) = j
           if (ptnm(i)(1:2).eq.char1(1:2).and.char1(3:3).eq.'*')
     &                ipick(i) = j
           if (ptnm(i)(1:3).eq.char1(1:3).and.char1(4:4).eq.'*')
     &                ipick(i) = j
5          continue
          else if (char1.eq.'mono') then
           call get4c(char1,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'missing monomer name',21,level)
           end if
           if (char1(1:4).eq.'none') then
            level = 1
            call alert(name,namel,'Illegal chem pickup',19,level)
           end if
           do 7 i=1,totmon
            if ((moname(i).eq.char1 .or. char1(1:1).eq.'*').or.
     &       (moname(i)(1:1).eq.char1(1:1).and.char1(2:2).eq.'*').or.
     &       (moname(i)(1:2).eq.char1(1:2).and.char1(3:3).eq.'*').or.
     &       (moname(i)(1:3).eq.char1(1:3).and.char1(4:4).eq.'*')) then
             if (i.eq.1) then
              i1 = 1
             else
              i1 = poipt(i-1)+1
             end if
             do 6 k=i1,poipt(i)
              ipick(k) = j
6            continue
            end if
7          continue
          end if
         end if
         if (debug) then
          write(stdo,*)' ipick in pick2 =',(ipick(i),i=1,npt)
         end if
         return
         end
