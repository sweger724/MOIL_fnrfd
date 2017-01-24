	subroutine getpd(ucrd,grdp,pstep,coor,coor1,coor2,pointr
     1    ,nselec,mlsti,igrid)
       
C       get inital structure for path calculations, the derivative of
C       the path at current position and the step vector to next
C       position. 
c
c       INPUT:
c       ucrd - a  unit number of a file with path binary coordinates
c       mlsti - milestone number
c       igrid  - the total number of path segements connecting igrid+1
c                structures *** NOTE THAT IN ucrd IGRID+1 POINTS ARE
c                REQUIRED (INCLUDING PRODUCTS) ***
c       nselec - number of selected atoms
c       pointr - is the pointer vector to the position of the selected
c                atoms
c       npt - number of atoms
c       
c       OUTPUT:
c       coor, coor1, coor2 [3 nselec]: milestone structures. The nselec
c             adjustment is made at the end of the subroutine.  (see code)
c       grdp [3 nselec] - reaction coordinate slope at mlsti, estimated
c            here by finite difference, NOT NORMALIZED.
c       pstep [3 nselec] - a step vector translation from point mlsti to
c             mlsti+1


	implicit none

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	integer ucrd,mlsti,igrid,nselec,k
	integer pointr(*)
	double precision grdp(3,*),pstep(3,*)
	double precision coor(3,*)
	double precision coor1(3,*)
	double precision coor2(3,*)
c	double precision pnorm, gnorm

c       local
	integer i,j,level
	double precision trash


	write(*,*) '(getpd) mlsti, igrid:', mlsti, igrid
c	write(*,*) '(getpd) npick:', nselec
c	write(*,*) '(getpd) pointr:', (pointr(i),i=1,npt)


	if ( ( mlsti .lt. 1 ) .or. ( mlsti .gt. igrid ) ) then
	   write(stdo,*)' *** error in getpd ***'
	   write(stdo,*)' mlsti, igrid: ',mlsti,igrid
	   level = 1
	   call alert('getpd',10,'Inconsistent data',16,level)
	end if


	rewind ucrd


c       first milestone
	if (mlsti.eq.1) then
	   read(ucrd,err=999,end=999)trash,(coor(1,i),i=1,npt),
     1	       (coor(2,i),i=1,npt),(coor(3,i),i=1,npt)
	   read(ucrd,err=999,end=999)trash,(coor1(1,i),i=1,npt),
     1	       (coor1(2,i),i=1,npt),(coor1(3,i),i=1,npt)

	   do 1 j=1,nselec
	      i = pointr(j)
	      do 11 k = 1,3
		 pstep(k,j) = coor1(k,i)-coor(k,i)
		 grdp(k,j) = pstep(k,j)
 11	      continue
 1	   continue


c       all other milestones
	else
	   do 2 j=1,mlsti-1
	      read(ucrd,err=999)trash,(coor1(1,i),i=1,npt),
     1            (coor1(2,i),i=1,npt),(coor1(3,i),i=1,npt)
 2	   continue
	   read(ucrd,err=999)trash,(coor(1,i),i=1,npt),
     1	       (coor(2,i),i=1,npt),(coor(3,i),i=1,npt)
	   read(ucrd,err=999)trash,(coor2(1,i),i=1,npt),
     1	       (coor2(2,i),i=1,npt),(coor2(3,i),i=1,npt)

	   do 3 j = 1,nselec
	      i = pointr(j)
	      do 14 k = 1,3
		 pstep(k,j) = coor2(k,i) - coor(k,i)
		 grdp(k,j)  = coor2(k,i) - coor1(k,i)
 14	      continue
 3	   continue
	end if



c       now shift the peptide coordinates to the beginning of
c       coor, coor1, coor2, as advertised
	do 171 j = 1,nselec
	   i = pointr(j)
	   do 172 k = 1,3
	      coor(k,j) = coor(k,i)
	      coor1(k,j) = coor1(k,i)
	      coor2(k,j) = coor2(k,i)
 172	   continue
 171	continue



c	normalize pstep and grdp
c	pnorm = 0.d0
c	gnorm = 0.d0
c	do 194 j = 1,nselec
c	   do 195 k = 1,3
c	      pnorm = pnorm + pstep(k,j) * pstep(k,j)
c	      gnorm = gnorm + grdp(k,j) * grdp(k,j)
c 195	   continue
c 194	continue
c	pnorm = 1.d0 / sqrt(pnorm)
c	gnorm = 1.d0 / sqrt(gnorm)
c	do 294 j = 1,nselec
c	   do 295 k = 1,3
c	      pstep(k,j) = pstep(k,j) * pnorm
c	      grdp(k,j) = grdp(k,j) * gnorm
c 295	   continue
c 294	continue

	
	
	return
	
 999	continue
	write(*,*) '*** error while reading path coordinates (getpd)'
	stop
	end
