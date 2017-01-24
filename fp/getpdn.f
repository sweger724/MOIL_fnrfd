	subroutine getpdn(ucrd,grdp,pstep,coor,coor1,coor2,mlsti,igrid,n)
       
C       get inital structure for path calculations, the derivative of
C       the path at current position and the step vector to next
C       position. NOTE: This version of getpd assumes that the
C       coordinate file has no water. The number of atoms per structure
C       is given by the argument n.
c
c       INPUT:
c       ucrd - a  unit number of a file with path binary coordinates
c       mlsti - milestone number
c       igrid  - the total number of path segements connecting igrid+1
c                structures *** NOTE THAT IN ucrd IGRID+1 POINTS ARE
c                REQUIRED (INCLUDING PRODUCTS) ***
c       n - number of atoms per structure
c       
c       OUTPUT:
c       coor, coor1, coor2 [3 n]: milestone structures. The nselec
c             adjustment is made at the end of the subroutine.  (see code)
c       grdp [3 n] - reaction coordinate slope at mlsti, estimated
c            here by finite difference, NOT NORMALIZED.
c       pstep [3 n] - a step vector translation from point mlsti to
c             mlsti+1


	implicit none

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	integer ucrd,mlsti,igrid,nselec,k,n
	double precision grdp(3,*),pstep(3,*)
	double precision coor(3,*)
	double precision coor1(3,*)
	double precision coor2(3,*)

c       local
	integer i,j,level
	double precision trash

	write(*,*) '(getpdn) mlsti, igrid, n:', mlsti, igrid, n

	if ( ( mlsti .lt. 1 ) .or. ( mlsti .gt. igrid ) ) then
	   write(stdo,*)' *** error in getpdn ***'
	   level = 1
	   call alert('getpdn',10,'Inconsistent data',16,level)
	end if

	rewind ucrd

c       first milestone
	if (mlsti.eq.1) then
	   read(ucrd,err=999,end=999)trash,(coor(1,i),i=1,n),
     1	(coor(2,i),i=1,n),(coor(3,i),i=1,n)
	   read(ucrd,err=999,end=999)trash,(coor1(1,i),i=1,n),
     1	(coor1(2,i),i=1,n),(coor1(3,i),i=1,n)

	   do 1 j=1,n
	      do 11 k = 1,3
		 pstep(k,j) = coor1(k,j)-coor(k,j)
		 grdp(k,j) = pstep(k,j)
 11	      continue
 1	   continue


c       all other milestones
	else
	   do 2 j=1,mlsti-1
	      read(ucrd,err=999)trash,(coor1(1,i),i=1,n),(coor1(2,i),i
     1	   =1,n),(coor1(3,i),i=1,n)
 2	   continue
	   read(ucrd,err=999)trash,(coor(1,i),i=1,n),
     1	(coor(2,i),i=1,n),(coor(3,i),i=1,n)
	   read(ucrd,err=999)trash,(coor2(1,i),i=1,n),
     1	(coor2(2,i),i=1,n),(coor2(3,i),i=1,n)

	   do 3 j = 1,nselec
	      do 14 k = 1,3
		 pstep(k,j) = coor2(k,j) - coor(k,j)
		 grdp(k,j)  = coor2(k,j) - coor1(k,j)
 14	      continue
 3	   continue
	end if


	return
	
 999	continue
	write(*,*) '*** error while reading path coordinates (getpdn)'
	stop
	end
