       subroutine getmlst(prev,next,pp,pn,nselec)

c       Get a milestone point and its normal.
c
c       INPUT:
c       
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

       
       integer i,j,k,level,nselec
       double precision e, pn(3,*), pp(3,*)
       double precision next(3,maxpt),prev(3,maxpt)

c       do i = 1,nselec
c       write(6,*) i,prev(1,i),prev(2,i),prev(3,i)
c       write(6,*) i,next(1,i),next(2,i),next(3,i)
c       enddo


       call vecmin ( next, prev, nselec*3, pn )

       call normalize_vec( pn, 3*nselec )

       return
 999       continue
       write(*,*) '*** error while reading path coordinates (getmlst)'
       stop
       end
