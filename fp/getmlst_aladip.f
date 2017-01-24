	subroutine getmlst(ucrd,pp,pn,mi,mt,pointr,nselec)

c       !! WARNING !!: This function clobbers velo (in VELOC.BLOCK),
c       so don't call it when the contents of velo are meaningful
c       
C       Get a milestone point and its normal. ucrd (path format) is
C       assumed to have the following structure: mlst_1_point,
C       mlst_1_normal, mlst_2_point, mlst_2_normal.
c
c       INPUT:
c       ucrd - a unit number of a file with path binary coordinates
c       mi - milestone number
c       mt - total number of milestones
c       pointr - selection pointer
c       nselec - number of selected atoms (length of pointr)
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

	integer ucrd,mi,mt,nselec
	integer i,j,k,pointr(*),level
	double precision e, pn(3,*), pp(3,*)

c       ----------------------------------------------------------
c       preliminaries

	write(stdo,*)' mi, mt: ',mi,mt

	if ( (mi .lt. 1) .or. (mi .gt. mt) ) then
	   write(stdo,*)' *** error in getmlst ***'
	   level = 1
	   call alert('getmlst',10,'Invalid mlst num',20,level)
	end if

	rewind ucrd

c       ---------------------------------------------------------

	do 1 j = 1,mi-1
	   read(ucrd,err=999,end=999) e,((velo(k,i),i=1,npt),k=1,3)
 1	continue

	read(ucrd,err=999,end=999) e,((pp(k,i),i=1,npt),k=1,3)

	if ( mi .eq. 1 ) call vdcopy(pp,velo,npt)
	   
	if ( mi .eq. mt ) then
	   call vdcopy(pp,coor,npt)
	else
	   read(ucrd,err=999,end=999) e,((coor(k,i),i=1,npt),k=1,3)
	end if

	call deselect ( velo, pointr, nselec )
	call vecmin_fp ( coor, velo, pn, nselec, pointr )
	call deselect ( pp, pointr, nselec )
	
	return

 999	continue
	write(*,*) '*** error while reading path coordinates (getpd)'
	stop

	end
