      subroutine Communicate_Positions(npt,rr,positions)
     
        implicit none

        integer npt
        double precision rr(3,*)
        logical positions
    
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
    
        integer l,i
    
        if (first) then
           do i =1,npt
             do l=1,3
               if (positions) then
                  rr(l,i)=r_initial(l,i)
               else
                  rr(l,i)=0.d0
               endif
             end do
           end do
        endif

        if (last) then
         do i =1,npt
           do l=1,3
             if (positions) then
                rr(l,(pseg+1)*npt+i)=r_final(l,i)
             else
                rr(l,(pseg+1)*npt+i)=0.d0
             endif
           end do
         end do
      endif

      if (mod(procID,2) .eq. 0) then
        call Send_first(3*npt,rr)
        call Receive_first(3*npt,rr)
      else
        call Receive_first(3*npt,rr)
        call Send_first(3*npt,rr)
      endif
 
      return
      end

