        subroutine wdyncrd(uwcrd,nfrm,ncoor,inofrz,nofreez,bin)
        implicit none
c writing dwn dynamics coordinate file
c currently the CHARM format is employed for compatability
c with graphic program (QUANTA)
c       
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        integer uwcrd,nfrm,ncoor
        integer inofrz,nofreez(*)
c        logical bin

        integer vari(20)
        character*7 name
        character*80 dummy
        integer namel,i,bin

        dummy = 'A dummy title to dynamics coordinates'
        name = 'wdyncrd'
        namel = 7


        if (ncoor.eq.1) then
                do 1 i=1,20
                        vari(i) = 0
1               continue
c
c first writing of coordinates. A title and other junk
c is required
c
c vari(1) is the number of coordinate files to be written
c vari(9) is the number of FIXED atoms
c
                vari(1) = nfrm
                vari(9) = npt - inofrz

		if (bin.ne.0) then
                  write(uwcrd)'CORD',vari
                  write(uwcrd)1,dummy
                  write(uwcrd)npt
                  if (inofrz.ne.npt) then
                   write(uwcrd)(nofreez(i),i=1,inofrz)
                  end if

c First time write all coordinates, even if some of them
c are fixed
c
                  write(uwcrd)(sngl(coor(1,i)),i=1,npt)
                  write(uwcrd)(sngl(coor(2,i)),i=1,npt)
                  write(uwcrd)(sngl(coor(3,i)),i=1,npt)

                else
                  write(uwcrd,*)'CORD',vari
                  write(uwcrd,*)1,dummy
                  write(uwcrd,*)npt
                  if (inofrz.ne.npt) then
                   write(uwcrd,*)(nofreez(i),i=1,inofrz)
                  end if
c First time write all coordinates, even if some of them
c are fixed
c
                  write(uwcrd,*)(sngl(coor(1,i)),i=1,npt)
                  write(uwcrd,*)(sngl(coor(2,i)),i=1,npt)
                  write(uwcrd,*)(sngl(coor(3,i)),i=1,npt)
                endif
        else
                if (inofrz.ne.npt) then
		 if (bin.ne.0) then
                  write(uwcrd) (sngl(coor(1,nofreez(i))),i=1,inofrz)
                  write(uwcrd) (sngl(coor(2,nofreez(i))),i=1,inofrz)
                  write(uwcrd) (sngl(coor(3,nofreez(i))),i=1,inofrz)
                 else
                  write(uwcrd,*) (sngl(coor(1,nofreez(i))),i=1,inofrz)
                  write(uwcrd,*) (sngl(coor(2,nofreez(i))),i=1,inofrz)
                  write(uwcrd,*) (sngl(coor(3,nofreez(i))),i=1,inofrz)
                 endif 

                else
		 if (bin.ne.0) then
                  write(uwcrd) (sngl(coor(1,i)),i=1,npt)
                  write(uwcrd) (sngl(coor(2,i)),i=1,npt)
                  write(uwcrd) (sngl(coor(3,i)),i=1,npt)
                 else
                  write(uwcrd,*) (sngl(coor(1,i)),i=1,npt)
                  write(uwcrd,*) (sngl(coor(2,i)),i=1,npt)
                  write(uwcrd,*) (sngl(coor(3,i)),i=1,npt)
                 endif
                end if
        end if
        return
        end
