        subroutine count_waters()
c
c counting the number of water molecules (TIP3) and storing the results
c in nwaters - a connect variable
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'

        integer i,j

        nwaters = 0

        if (mdivyes) then
        totdmon = 0

c first figure out the number of divisions

        do 2 i=1,totmon
         if (mdivlist(0,i).eq.-1) then
                totdmon = totdmon + 1
                realmono(totdmon) = i
         end if
         do 1 j=1,mdivlist(0,i)
                totdmon = totdmon + 1
                realmono(totdmon) = i   
1        continue
2       continue

        write(*,*)' Number of monomer division = ',totdmon
         

        do 3 i=1,totdmon
                if (moname(realmono(i)).eq.'TIP3' .or.
     1              moname(realmono(i))(1:3).eq.'SPC'.or.
     2              moname(realmono(i)).eq.'TIP4') then
                         nwaters = nwaters + 1
                         idxtip3(nwaters) = i
                end if
3       continue
        else

        do 4 i=1,totmon
                if (moname(i).eq.'TIP3' .or.
     1              moname(i)(1:3).eq.'SPC'.or.
     2              moname(i).eq.'TIP4') then
                        nwaters = nwaters + 1
                        idxtip3(nwaters) = i
                end if
4       continue
        end if

        return
        end
