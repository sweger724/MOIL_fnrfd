       subroutine getmlst(ucrd,pp,pn,mi,mt,nselec,pointr)

c       Get a milestone point and its normal. ucrd (path format) is
c       assumed to have the following structure: mlst_1_point,
c       mlst_1_normal, mlst_2_point, mlst_2_normal.
c
c       INPUT:
c       ucrd - a unit number of a file with path binary coordinates
c       mi - milestone number
c       mt - total number of milestones
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

       integer ucrd,mi,mt
       integer i,j,k,level,nselec,pointr(*)
       double precision e, pn(3,*), pp(3,*)
       double precision next(3,maxpt),prev(3,maxpt)


       write(stdo,*)' mi, mt: ',mi,mt

       if ( (mi .lt. 1) .or. (mi .gt. mt) ) then
          write(stdo,*)' *** error in getmlst ***'
          level = 1
          call alert('getmlst',10,'Invalid mlst num',20,level)
       end if
       rewind ucrd


       do 1 j = 1,mi-1
          read(ucrd,err=999,end=999) e,((prev(k,i),i=1,npt),k=1,3)
 1       continue
       read(ucrd,err=999,end=999) e,((pp(k,i),i=1,npt),k=1,3)
       if ( mi .eq. 1 ) call vdcopy(pp,prev,npt)
       if ( mi .eq. mt ) then
          call vdcopy(pp,next,npt)
       else
          read(ucrd,err=999,end=999) e,((next(k,i),i=1,npt),k=1,3)
       end if

       call vecmin ( next, prev, npt*3, pn )

       call deselect ( pp, pointr, nselec )
       call deselect ( pn, pointr, nselec )
       call normalize_vec( pn, 3*nselec )

       return
 999       continue
       write(*,*) '*** error while reading path coordinates (getmlst)'
       stop
       end
