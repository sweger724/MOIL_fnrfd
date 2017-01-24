      subroutine add_muta(ipick,igroup)
c Subroutine for give an index to mutated particles

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/PROPERT.BLOCK'
        include 'COMMON/MONOMERS.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/TMP_CONNECT.BLOCK'
        include 'COMMON/MUTA.BLOCK'

       integer ipick(maxpt),igroup,i

!       call pick(ipick,igroup)
       call rmute(ipick)
       do i=1,npt
        if (ipick(i).gt.0) then
         mutaid(i)=ipick(i)
        else
         mutaid(i)=0
        endif
       enddo
     
      return
      end
