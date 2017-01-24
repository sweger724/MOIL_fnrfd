	subroutine getpd_umbr(ucrd,grdp,pstep,coor,
     1		coor1,coor2,pointr,nselec,
     2		numpth,igrid,rvrs)
c
c get inital structure for path calculations, the derivative of the
c path at current position and the step vector to next position.
c variables:
c input- 
c ucrd - a unit number of a file with path binary coordinates
c numpth - the number of path structure of interest
c igrid  - the total number of path segements connecting igrid+1
c	   structures
c	*** NOTE THAT IN ucrd IGRID+1 POINTS ARE REQUIRED (INCLUDING
c		PRODUCTS) ***
c nselec - number of selected atoms
c pointr - is the pointer vector to the position of the selected
c          atoms
c natom  - number of atoms
c rvrs   - logical variable if .true. start the reaction from prod.
c
c output-
c grdp     - slopes of the curvlinear coordinate describing
c		the reaction coordinate at position numpth.
c		Estimated here by finite difference, NOT NORMALIZED.
c pstep    - a step vector translation from point numpth to
c		numpth+1
c x y z       - initial coordinates (minimum energy point) for
c		free energy sampling at numpth
c
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	integer ucrd,numpth,igrid,natom,nselec,level
	integer pointr(*)
	double precision grdp(3,*),pstep(3,*)
	double precision coor(3,*)
	double precision coor1(3,*)
	double precision coor2(3,*)
	logical rvrs

c
c local
c
	integer i,j,nsel3
	double precision trash

	if (numpth.gt.igrid) then
		write(stdo,*)' *** error in getpd ***'
		write(stdo,*)' numpth = path =',numpth,igrid
		level = 1
		call alert('getpd_umbr',10,'Incosistent data',16,level)
	end if
	natom = npt
	nsel3=3*nselec
	rewind ucrd
c -- reactant  --> product
c
	if (.not. rvrs) then
c
c special treatment is required for first point.
c
	if (numpth.eq.1) then
	  read(ucrd,err=999,end=999)trash,(coor(1,i),i=1,natom),
     1		(coor(2,i),i=1,natom),(coor(3,i),i=1,natom)
	  read(ucrd,err=999,end=999)trash,(coor1(1,i),i=1,natom),
     1		(coor1(2,i),i=1,natom),(coor1(3,i),i=1,natom)
C$DOIT IVDEP
	  do 1 j=1,nselec
		i=pointr(j)
		pstep(1,j)=coor1(1,i)-coor(1,i)
		pstep(2,j)=coor1(2,i)-coor(2,i)
		pstep(3,j)=coor1(3,i)-coor(3,i)
		grdp(1,j) =pstep(1,j)
		grdp(2,j) =pstep(2,j)
		grdp(3,j) =pstep(3,j)
1	  continue
	  return
	end if
	do 2 j=1,numpth-1
	  read(ucrd,err=999)trash,(coor1(1,i),i=1,natom),
     1		(coor1(2,i),i=1,natom),(coor1(3,i),i=1,natom)
2	continue
	read(ucrd,err=999)trash,(coor(1,i),i=1,natom),
     1		(coor(2,i),i=1,natom),(coor(3,i),i=1,natom)
	read(ucrd,err=999)trash,(coor2(1,i),i=1,natom),
     1		(coor2(2,i),i=1,natom),(coor2(3,i),i=1,natom)
C$DOIT IVDEP
	do 3 j=1,nselec
		i=pointr(j)
		pstep(1,j)=coor2(1,i)-coor(1,i)
		pstep(2,j)=coor2(2,i)-coor(2,i)
		pstep(3,j)=coor2(3,i)-coor(3,i)
		grdp(1,j)=coor2(1,i)-coor1(1,i)
		grdp(2,j)=coor2(2,i)-coor1(2,i)
		grdp(3,j)=coor2(3,i)-coor1(3,i)
3	continue
c product --> reactants
c
	else
c
c special treatment is required for first (last) point.
c
	if (numpth.eq.1) then
c An important reminder: igrid+1 structures are expected.
c
	do 20 j=1,igrid-1
	  read(ucrd,err=999)trash,(coor(1,i),i=1,natom),
     1		(coor(2,i),i=1,natom),(coor(3,i),i=1,natom)
20	continue
	  read(ucrd,err=999,end=999)trash,(coor1(1,i),i=1,natom),
     1	  	(coor1(2,i),i=1,natom),(coor1(3,i),i=1,natom)
	  read(ucrd,err=999,end=999)trash,(coor(1,i),i=1,natom),
     1	  	(coor(2,i),i=1,natom),(coor(3,i),i=1,natom)
C$DOIT IVDEP
	  do 21 j=1,nselec
		i=pointr(j)
		pstep(1,j)=coor1(1,i)-coor(1,i)
		pstep(2,j)=coor1(2,i)-coor(2,i)
		pstep(3,j)=coor1(3,i)-coor(3,i)
		grdp(1,j)=pstep(1,j)
		grdp(2,j)=pstep(2,j)
		grdp(3,j)=pstep(3,j)
21	  continue
	  return
	end if
	do 22 j=1,igrid-numpth+1
	  read(ucrd,err=999)trash,(coor1(1,i),i=1,natom),
     1		(coor1(2,i),i=1,natom),(coor1(3,i),i=1,natom)
22	continue
	read(ucrd,err=999)trash,(coor(1,i),i=1,natom),
     1		(coor(2,i),i=1,natom),(coor(3,i),i=1,natom)
	read(ucrd,err=999)trash,(coor2(1,i),i=1,natom),
     1		(coor2(2,i),i=1,natom),(coor2(3,i),i=1,natom)
C$DOIT IVDEP
	do 23 j=1,nselec
		i=pointr(j)
		pstep(1,j)=coor1(1,i)-coor(1,i)
		pstep(2,j)=coor1(2,i)-coor(2,i)
		pstep(3,j)=coor1(3,i)-coor(3,i)
		grdp(1,j)=coor1(1,i)-coor2(1,i)
		grdp(2,j)=coor1(2,i)-coor2(2,i)
		grdp(3,j)=coor1(3,i)-coor2(3,i)
23	continue
	end if
	return
999	continue
	write(*,*) '*** error while reading path coordinates (getpd)'
	stop
	end
