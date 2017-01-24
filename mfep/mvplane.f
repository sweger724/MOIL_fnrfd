       subroutine mvplane(prev,next,pn,nselec)

c       Get structures i-1, i+1 and modify the plane for string i  
c       
c       INPUT:
c       prev(3,nselec) averge coordinates of the selected particles at the i-1th plane 
c       next(3,nselec) averge coordinates of the selected particles at i+1th plane 
c       
c       OUTPUT:
c       pp - milestone coordinate ("plane point")
c       pn - milestone normal ("plane normal")

 
       implicit none

       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/ENERGY.BLOCK'
       include 'COMMON/NBLIST.BLOCK'
       include 'COMMON/UNITS.BLOCK'
       include 'COMMON/COORD.BLOCK'
       include 'COMMON/VELOC.BLOCK'
       include 'COMMON/LINE.BLOCK'
       include 'COMMON/DEBUG.BLOCK'

       integer nselec,j
       double precision pn(3,*)
       double precision next(3,*),prev(3,*)


       
       call vecmin ( next, prev, nselec*3, pn )
       call normalize_vec( pn, 3*nselec )

c       do j=1,nselec      
c       write(stdo,14) pn(1,j),pn(2,j),pn(3,j)
c       enddo 
 
c14    format("mvplane:",3f14.8)
       return
       end
