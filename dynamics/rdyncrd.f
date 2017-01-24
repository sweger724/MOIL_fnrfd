        subroutine rdyncrd(urcrd,ncoor,inofrz,nofreez,rbin)
        implicit none
c reading dynamics coordinate file
c currently the CHARM format is employed for compatability
c with a graphic program (QUANTA). However since the options
c here are more limited than in Charmm, we shall not be able to
c read all possible dynamics files written on charmm
c       
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        integer urcrd,ncoor,nptx,rbin
        integer inofrz,nofreez(*)

        integer vari(20),i,j
        real r(maxpt)
        character*7 name
        character*4 head
        character*1 dummy
        integer namel,level


        name = 'rdyncrd'
        namel = 7
        write(stdo,*) ' READING DYNAMICS CRD, STRUCTURE # ',NCOOR

                 do i=1,npt
                 chain(i) = '*'
                 enddo

c
c first reading of coordinates. A title and other junk
c is required
c
c vari(1) is the number of coordinate files to be read
c vari(9) is the number of FIXED atoms
c
        if ((.not.norew).or.ncoor.eq.lpstr) then
          if (rbin.ne.0) then
                read(urcrd,err=900)head,vari
                read(urcrd,err=900)i,dummy
                read(urcrd,err=900)nptx
                if (npt.ne.nptx) then
                 level = 1
                 call alert(name,namel,'# of particles dont match'
     1             ,25,level)
                end if
                inofrz = npt - vari(9)
 
                if (inofrz.ne.npt .and. inofrz.ne.0) then
                 read(urcrd)(nofreez(i),i=1,inofrz)
                end if
c First time read all coordinates, even if some of them
c are fixed
c
                read(urcrd,err=901)(r(i),i=1,npt)
                do 1 i=1,npt
                 coor(1,i) = dble(r(i))
1               continue

                read(urcrd)(r(i),i=1,npt)
                do 2 i=1,npt
                 coor(2,i) = dble(r(i))
2               continue

                read(urcrd)(r(i),i=1,npt)
                do 3 i=1,npt
                 coor(3,i) = dble(r(i))
3               continue
         else
                read(urcrd,*,err=900)head,vari
                read(urcrd,*,err=900)i,dummy
                read(urcrd,*,err=900)nptx
                if (npt.ne.nptx) then
                 level = 1
                 write(*,*)' npt (wcon) nptx (dcd) ',npt,nptx
                 call alert(name,namel,'# of particles dont match'
     1             ,25,level)
                end if
c@                inofrz = npt - vari(9)
 
                if (inofrz.ne.npt .and. inofrz.ne.0) then
                 read(urcrd,*)(nofreez(i),i=1,inofrz)
                end if
c First time read all coordinates, even if some of them
c are fixed
c
                read(urcrd,*,err=901)(r(i),i=1,npt)
                do 11 i=1,npt
                 coor(1,i) = dble(r(i))
11               continue

                read(urcrd,*,err=901)(r(i),i=1,npt)
                do 12 i=1,npt
                 coor(2,i) = dble(r(i))
12               continue

                read(urcrd,*,err=901)(r(i),i=1,npt)
                do 13 i=1,npt
                 coor(3,i) = dble(r(i))
13               continue
         endif
        endif

                
        if (ncoor.eq.1) return

        if (.not.norew.or.lpstr.eq.ncoor) then
         if (rbin.ne.0) then
                do 35 i=1,ncoor-2
                 read(urcrd)r(1)
                 read(urcrd)r(1)
                 read(urcrd)r(1)
35              continue
         else
                do 135 i=1,ncoor-2
                 if (inofrz.ne.npt .and. inofrz.ne.0) then
                  read(urcrd,*)(r(j),j=1,inofrz)
                  read(urcrd,*)(r(j),j=1,inofrz)
                  read(urcrd,*)(r(j),j=1,inofrz)
                 else
                  read(urcrd,*)(r(j),j=1,npt)
                  read(urcrd,*)(r(j),j=1,npt)
                  read(urcrd,*)(r(j),j=1,npt) 
                 endif
135             continue
         endif 
        endif

       if (rbin.ne.0) then
         if (inofrz.ne.npt .and. inofrz.ne.0) then
                read(urcrd)(r(i),i=1,inofrz)
                do 4 i=1,inofrz
                 coor(1,nofreez(i)) = dble(r(i))
4               continue

                read(urcrd)(r(i),i=1,inofrz)
                do 5 i=1,inofrz
                 coor(2,nofreez(i)) = dble(r(i))
5               continue

                read(urcrd)(r(i),i=1,inofrz)
                do 6 i=1,inofrz
                 coor(3,nofreez(i)) = dble(r(i))
6               continue
          else
c
                read(urcrd)(r(i),i=1,npt)
                do 7 i=1,npt
                 coor(1,i) = dble(r(i))
7               continue

                read(urcrd)(r(i),i=1,npt)
                do 8 i=1,npt
                 coor(2,i) = dble(r(i))
8               continue

                read(urcrd)(r(i),i=1,npt)
                do 9 i=1,npt
                 coor(3,i) = dble(r(i))
9               continue
         end if
        else
         if (inofrz.ne.npt .and. inofrz.ne.0) then
                read(urcrd,*)(r(i),i=1,inofrz)
                do 14 i=1,inofrz
                 coor(1,nofreez(i)) = dble(r(i))
14               continue

                read(urcrd,*)(r(i),i=1,inofrz)
                do 15 i=1,inofrz
                 coor(2,nofreez(i)) = dble(r(i))
15               continue

                read(urcrd,*)(r(i),i=1,inofrz)
                do 16 i=1,inofrz
                 coor(3,nofreez(i)) = dble(r(i))
16               continue
          else
c
                read(urcrd,*)(r(i),i=1,npt)
                do 17 i=1,npt
                 coor(1,i) = dble(r(i))
17               continue

                read(urcrd,*)(r(i),i=1,npt)
                do 18 i=1,npt
                 coor(2,i) = dble(r(i))
18               continue

                read(urcrd,*)(r(i),i=1,npt)
                do 19 i=1,npt
                 coor(3,i) = dble(r(i))
19               continue
         end if
        endif

        return
900     level = 1
        call alert(name,namel,'Error in reading title',22,level)
        return
901     level = 1
        call alert(name,namel,'Error in reading crd#1',22,level)
        return
        end
