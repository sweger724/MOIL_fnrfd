       subroutine edcon(uedit)
c
c Edit connectivity file (now used only to remove some angles)
c
       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/UNITS.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/MONOMERS.BLOCK'
       include 'COMMON/PROPERT.BLOCK'
       include 'COMMON/LINE.BLOCK'
       include 'COMMON/DEBUG.BLOCK'

       integer uedit
       logical find
       integer geti,intgr
c local
       character*5 name
       character*4 m1,m2,m3,pt1,pt2,pt3
       integer im1,im2,im3
       integer ist1,ist2,ist3,ien1,ien2,ien3
       integer namel,id1,id2,id3,level,i,j
       logical empty
       data name/'edcon'/
       data namel/5/

       if (uedit.ne.5) rewind uedit
1       continue
       call rline(name,namel,uedit)
       if (find('*EOD')) return
       if (find('remo')) then

c ================================================================
        if (find('bond')) then
         id1=geti('atm1',-1)
         id2=geti('atm2',-1)
         if (id1.gt.0 .and. id2. gt.0) go to 550
         if (id1.lt.0 .or. id2.lt.0) then
          if (find('chem')) then
           call get4c(m1,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'Missing mono name',17,level)
           end if
           im1 = intgr()
           call get4c(pt1,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'Missing atom name',17,level)
           end if
           call get4c(m2,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'Missing mono name',17,level)
           end if
           im2 = intgr()
           call get4c(pt2,empty)
           if (empty) then
            level = 1
            call alert(name,namel,'Missing atom name',17,level)
           end if 
c end if(CHEM)
        end if
           if (m1.ne.moname(im1) .or. m2.ne.moname(im2)) then
            level = 1
            write(stdo,100)m1,moname(im1),m2,moname(im2)
100         format(1x,'monomer name does not match, read ',a4,
     1          ' expected ',a4,' read ',a4,' expected ',a4)
            call alert(name,namel,'Mono do not match',17,level)
           end if
           ien1 = poipt(im1)
           if (im1.eq.1) then
            ist1 = 1
           else
            ist1 = poipt(im1-1)+1
           end if

           ien2 = poipt(im2)
           if (im2.eq.1) then
            ist2 = 1
           else
            ist2 = poipt(im2-1)+1
           end if

           do 200 i=ist1,ien1
            if(pt1.eq.ptnm(i)) then
             id1 = i
             go to 300
            end if
200        continue
           level = 1
           call alert(name,namel,'Particle not found',18,level)
300        continue

           do 400 i=ist2,ien2
            if(pt2.eq.ptnm(i)) then
             id2 =i
             go to 500
            end if
400        continue
           level = 1
           call alert(name,namel,'Particle not found',18,level)
500        continue
          else
          level=1
          call alert(name,namel,'Illegal remo bond input',18,level)
c end if(id)
          end if

550          continue
         do 50 i=1,nb
          if (ib1(i) .eq. id1 .and. ib2(i).eq.id2) then
           do 20 j=i,nb-1
            ib1(j)=ib1(j+1)
            ib2(j)=ib2(j+1)
            kbond(j)=kbond(j+1)
            req(j)=req(j+1)
20          continue
            nb=nb-1
            else if (ib1(i) .eq. id2 .and. ib2(i) .eq. id1) then
            do 30 j=i,nb-1
             ib1(j)=ib1(j+1)
             ib2(j)=ib2(j+1)
             kbond(j)=kbond(j+1)
             req(j)=req(j+1)
30           continue
            nb=nb-1 
            end if
50         continue
c end if (bond)
          end if
c =================================================================          

        if (find('angl')) then
        id1 = geti('atm1',-1)
        id2 = geti('atm2',-1)
        id3 = geti('atm3',-1)
        if (id1.lt.0 .or. id2.lt.0 .or. id3.lt.0) then
         if (find('chem')) then
          call get4c(m1,empty)
          if (empty) then
           level = 1
           call alert(name,namel,'Missing mono name',17,level)
          end if
          im1 = intgr()
          call get4c(pt1,empty)
          if (empty) then
           level = 1
           call alert(name,namel,'Missing atom name',17,level)
          end if
          call get4c(m2,empty)
          if (empty) then
           level = 1
           call alert(name,namel,'Missing mono name',17,level)
          end if
          im2 = intgr()
          call get4c(pt2,empty)
          if (empty) then
           level = 1
           call alert(name,namel,'Missing atom name',17,level)
          end if 
          call get4c(m3,empty)
          if (empty) then
           level =1
           call alert(name,namel,'Missing mono name',17,level)
          end if
          im3 = intgr()
          call get4c(pt3,empty)
          if (empty) then
           level = 1
           call alert(name,namel,'Missing atom name',17,level)
          end if
          if (m1.ne.moname(im1) .or. m2.ne.moname(im2) .or.
     1           m3.ne.moname(im3)) then
           level = 1
           write(stdo,101)m1,moname(im1),m2,moname(im2),m3,moname(3)
101           format(1x,'monomer name does not match, read ',a4,
     1             ' expected ',a4,' read ',a4,' expected ',a4,/,
     2           ' expected ',a4,' read ',a4)
           call alert(name,namel,'Mono do not match',17,level)
          end if

          ien1 = poipt(im1)
          if (im1.eq.1) then
           ist1 = 1
          else
           ist1 = poipt(im1-1)+1
          end if

          ien2 = poipt(im2)
          if (im2.eq.1) then
           ist2 = 1
          else
           ist2 = poipt(im2-1)+1
          end if

          ien3 = poipt(im3)
          if (im3.eq.1) then
           ist3 = 1
          else
           ist3 = poipt(im3-1)+1
          end if

          do 2 i=ist1,ien1
           if(pt1.eq.ptnm(i)) then
            id1 = i
            go to 3
           end if
2          continue
          level = 1
          call alert(name,namel,'Particle not found',18,level)
3          continue
          
          do 4 i=ist2,ien2
           if(pt2.eq.ptnm(i)) then
            id2 =i
            go to 5
           end if
4          continue
          level = 1
          call alert(name,namel,'Particle not found',18,level)
5          continue
          do 6 i=ist3,ien3
           if(pt3.eq.ptnm(i)) then
            id3 = i
            go to 7
           end if
6          continue
          level = 1
          call alert(name,namel,'Particle not found',18,level)
7          continue
        else
         level = 1
         call alert(name,namel,'Illegal remo input',18,level)
        end if
        end if
        do 10 i=1,nangl
         if (iangl2(i).eq.id2) then
          if (iangl1(i).eq.id1 .and. iangl3(i).eq.id3) then
           do 8 j=i,nangl-1
            iangl1(j) = iangl1(j+1)
            iangl2(j) = iangl2(j+1)
            iangl3(j) = iangl3(j+1)
            kangl(j)  = kangl(j+1)
            angleq(j) = angleq(j+1)
8           continue
           nangl = nangl - 1
           go to 11
          else if (iangl1(i).eq.id3 .and. iangl3(i).eq.id2) then
           do 9 j=i,nangl-1
            iangl1(j) = iangl1(j+1)
            iangl2(j) = iangl2(j+1)
            iangl3(j) = iangl3(j+1)
            kangl(j)  = kangl(j+1)
            angleq(j) = angleq(j+1)
9           continue
           nangl = nangl - 1
           go to 11
          end if
         end if
10        continue
        level = 1
        call alert(name,namel,'remo angle not found',20,level)
11        continue
        end if
       end if
       go to 1
       end
