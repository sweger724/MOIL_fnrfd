           program average
c
c Xtract distance along a trajectory
c
		implicit none
	   include 'COMMON/LENGTH.BLOCK'
	   include 'COMMON/COORD.BLOCK'
	   include 'COMMON/UNITS.BLOCK'
	   include 'COMMON/CONNECT.BLOCK'
	   include 'COMMON/DEBUG.BLOCK'
	   include 'COMMON/LINE.BLOCK'
	   include 'COMMON/FREEZ.BLOCK'
	   include 'COMMON/CONVERT.BLOCK'
	   include 'COMMON/CCRD.BLOCK'
c ipick - pick subset of particles, a vector of length maxpt
c         value of zero if particle not selected one if it is
c npick  - number of picked particles
c of     - integer function for Opening a File, returned value is the
c          assigned unit number
c geti   - integer function to get integer value from a command line
c nstru  - number of structures in dynamics file
c namel  - length of program name
c inofrz - number of moving particles. used in reading the dynamics file
c level  - level of error found, level=0 only warning will be issued,
c          level=1 program stop
c urcrd,ucon - units of dynamics coord and connectivity files
c
	   integer ipick(maxpt),pointr(maxpt)
	   integer npick 
	   integer of,geti,nstru
	   integer namel,i,j,l,level,i1,i2
	   integer urcrd,uwcrd,ucon
           integer rbin
	   double precision dx,dy,dz,r

	   character*6 name
	   character*4 crd_type
	   character*5 getchar
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd,uwcrd/3*99/


	   lpstr = 1
	   norew = .false.

           stdi=5
	   stdo=6

           rbin = 1
	   totmon=0
	   npt=0
	   name='avera'
	   namel=5
	   crd_type='DYNA'
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    nstru= 1
	    pickpt=.false.

	call init_var()

1           continue
            call rline(name,namel,stdi)
	    if (find('norw')) norew=.true.
	    if (find('file')) then
	       if (find ('conn')) then  
		ucon=of()
c get connectivity 
	        call rconn(ucon)
c get coordinate file
	       else if (find ('rcrd')) then
		if (npt.eq.0) then
		 level = 1
		 call alert(name,namel,'Must read con file first',
     1			24,level)
		end if
	        urcrd=of()
	       else if (find('wcrd')) then
		uwcrd=of()
	       end if
              else 
	       nstru=geti('#str',nstru)
	       crd_type = getchar('ctyp',crd_type,4)
	       if (find('pick')) then
		call pick(ipick,npick)
		npick = 0
		do 2 i = 1,npt
		 if (ipick(i).gt.0) then
			npick = npick + 1
			pointr(npick) = i
		 end if
2		continue
               end if
	       if (find ('action')) goto  3
             end if
	     go to 1
3            continue

	do 35 i =1,npt
		do 35 j = 1,3
			coor2(j,i) = 0.d0
35	continue

          if (.not. fopen(ucon)) then
	     level=1
	     call alert(name,namel,'ucon not opened',15,level)
           else if (.not. fopen(urcrd)) then
	     level=1
	     call alert(name,namel,'urcrd not opened',16,level)
	   end if
	     
   
c read dynamics structures
	   rewind urcrd
	   j = 0
	   do 5 l=1,nstru
	    if (crd_type.eq.'DYNA') then
	    	if (.not.norew) rewind urcrd
	    	call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
	    else if (crd_type.eq.'PATH') then
		call rpath_seq(urcrd,l)
		write(*,*) ' coor(1,1) ', coor(1,1)
	    else 
		write(*,*) ' Ctyp = ',crd_type(1:4)
		call alert(name,namel,'Unrecognized ctyp',17,level)
	    end if
		do 4 i=1,npt
		 do 4 j = 1,3
		 coor2(j,i) = coor2(j,i) + coor(j,i)
4		continue
		write(*,*) ' coor2(1,1) ', coor2(1,1)
5	   continue
	   do 6 i=1,npt
		do 6 j = 1,3
			coor(j,i) = coor2(j,i)/nstru
6		continue
		write(*,*) ' coor(1,1) ', coor(1,1)
		call putcrd(uwcrd,'CHARM')

	   stop
	   end

