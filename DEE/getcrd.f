      subroutine getcrd(ucrd,type)
c
c subroutine to read coordinates
c currently very simple (almost no checks) version, without LES
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      integer ucrd
      character*5 type
c local
      character*4 char1,char2,char3,char4
      character*6 name
      character*80 zevel
      integer leng,namel,i,j,imoncmp,imon
      integer istart,i1,j1,k,l,level
      double precision xtmp,ytmp,ztmp,moretmp
      data name/'getcrd'/
      data namel/6/

      rewind ucrd
c initialize all coordinates to 9999.0
      do 1 i=1,npt
         do 1 j=1,3
            coor(j,i) = 9999.d0
 1    continue
      if (type.eq.'CHARM') then
         write(stdo,*)' Reading CHARMm coordinate file '
c
c     read the title and find out its length
c
         leng = 0
 2       continue
            read(ucrd,105)zevel
 105        format(80a)
            if (zevel(1:1).eq.'*') then
               do 21 k=80,1,-1
                  if (zevel(k:k).ne.' ') then
                     call echo(name,namel,zevel,k)
                     go to 22
                  end if
 21            continue
 22            continue
               leng = leng + 1
               go to 2
            end if
c     
c     check that title exist
c
         if (leng.eq.0) then
            level = 1
            call alert(name,namel,'Crd file with no title',22,level)
            return
         end if
c
c     don't forget to backspace the file !!
c
         backspace ucrd
         read(ucrd,*) i
         if (i.ne.0 .and. i.ne.npt) then
            level = 1
            call alert(name,namel,
     1           'Number of particles do not match',32,level)
            return
         end if
         istart = 1
         if (debug) then
            write(stdo,*)' totmon npt poipt '
     1           ,totmon,npt,(poipt(i),i=1,totmon)
         end if
         do 5 i=1,totmon
            do 4  j=istart,poipt(i)
c
c     read one line of coordinate
c
               read(ucrd,100,err=7,end=55)i1,j1,char1,char2,xtmp
     1              ,ytmp,ztmp,char3,char4,moretmp
 100           format(i5,i5,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)
c
c     check that this is the right monomers name
c
               if (char1.ne.moname(i)) then
c     
c     may be we just miss some atoms (this is ok sometime..
c     for example when hydrogens are not included
c
                  if (j1.ne.i) then
                     level = 0
                     write(stdo,101)i,moname(i)
 101                 format(1x,' Missing pts for monomer ',i5,
     1                    ' type = ',a4)
                     call alert(name,namel,'Missing particles',17,level)
                     backspace ucrd
                     go to 4
                  else if (j1.eq.totmon) then
                     level = 0
                     write(stdo,102)i,moname(i)
 102                 format(1x,' Missing pts for last monomer ',i5,
     1                    ' type = ',a4)
                     call alert(name,namel,'Missing particles',17,level)
                     backspace ucrd
                     go to 5
                  else
                     level = 1
                     write(stdo,103)i,char1,moname(i)
 103                 format(1x,'mono ',i5,' Read = ',a4,' Expected = ',
     1                    a4)
                     call alert(name,namel,'Monomers do not match',21,
     1                    level)
                  end if
               else
                  do 3 k=istart,poipt(i)
                     if (char2 .eq.  ptnm(k))then 
                        coor(1,k) = xtmp
                        coor(2,k) = ytmp
                        coor(3,k) = ztmp
                        more(k) = moretmp
                        go to 35
                     end if
 3                continue
 35               continue
               end if
 4          continue
            istart=poipt(i) + 1
 5       continue
c check that there are not any unidentified atoms
 55      continue
         do 6 i=1,npt
            if (coor(1,i).gt.9998.0) then
               write(stdo,104)i,ptnm(i),moname(poimon(i))
 104           format(1x,' Coordinates of particle ',i5,2x,a4,
     1              ' not defined in monomer ',a4)
               level = 0
               call alert(name,namel,'Undefined coordinates',21,level)
            end if
 6       continue
         return
      else if (type(1:3).eq.'pdb') then
         write(stdo,*)' Reading pdb coordinate file '
c     
c check that the current line is not ATOM line. If not simply echo the
c line and keep reading. Stop this loop when ATOM is detected backspace
c the file and start read thee coordinates
c
 8       continue
            read(ucrd,105)zevel
            if (zevel(1:4).ne.'ATOM' .and. zevel(1:4).ne.'HETA') then
               do 9 k=80,1,-1
                  if (zevel(k:k).ne.' ') then
                     call echo(name,namel,zevel,k)
                     go to 8
                  end if
 9             continue
               go to 8
            end if
c
         backspace ucrd
         imoncmp = -1
         imon    = -1
         istart = 1
         do 11 i=1,totmon
            do 12 j=istart,poipt(i)
c read one line of coordinates
               read(ucrd,106,end=75,err=7)
     1              zevel(1:4),j1,char2,char1,i1,xtmp,ytmp,ztmp,moretmp
               if (imoncmp.ne.i1) then
                  imon    = imon + 1
                  imoncmp = i1
                  if (j.ne.istart) then
                     backspace ucrd
                     go to 125
                  end if
               end if
 106           format(a4,2x,i5,2x,a4,a4,1x,i4,4x,3(f8.4),6x,f6.2)
               if (char1.ne.moname(i)) then
                  if (imon.ne.i) then
                     level = 0
                     write(stdo,101)i,moname(i)
                     call alert(name,namel,'Missing particles',17,level)
                     backspace ucrd
                     go to 12
                  else if (imon.eq.totmon) then
                     level = 0
                     write(stdo,102)i,moname(i)
                     call alert(name,namel,'Missing particles',17,level)
                     backspace ucrd
                     go to 12
                  else
                     level = 1
                     write(stdo,103)i,char1,moname(i)
                     call alert(name,namel,'Monomers do not match',21,
     1                    level)
                  end if
               else
                  do 10 k=istart,poipt(i)
                     if (char2.eq.ptnm(k)) then
                        coor(1,k) = xtmp
                        coor(2,k) = ytmp
                        coor(3,k) = ztmp
                        more(k) = moretmp
                     end if
 10               continue
               end if
 12         continue
 125        continue
            istart = poipt(i) + 1
 11      continue
c check that ther are not any unidentified particles
 75      continue
         do 13 i=1,npt
            if (coor(1,i).gt.9998.0d0) then
               write(stdo,104)i,ptnm(i),moname(poimon(i))
               level = 0
               call alert(name,namel,'Undefined coordinates',21,level)
            end if
 13      continue
         return
      end if
 7    continue
      level = 1
      call alert(name,namel,'Error during read',17,level)
      return
      end
      
