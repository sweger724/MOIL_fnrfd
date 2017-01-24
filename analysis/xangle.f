           program xangle
c
c Xtract angle along a trajectory
c
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
	   integer ipick(maxpt)
	   integer npick 
	   integer of,geti,nstru
	   integer namel,i,j,l,level,i1,i2,i3
	   integer urcrd,uwang,ucon
           integer rbin
	   double precision vec1x,vec1y,vec1z,vec2x,vec2y,vec2z
	   double precision cost,theta,norm1,norm2
	   double precision pi

	   character*6 name
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd,uwang/3*99/


	   i1    = 0
	   i2    = 0
	   i3    = 0
	   lpstr = 1
	   norew = .false.
	   pi    = 4.d0*datan(1.d0)

           stdi=5
	   stdo=6
           rbin = 1

	   totmon=0
	   npt=0
	   name='xangle'
	   namel=6
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    nstru= 1
	    pickpt=.false.

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
	       else if (find('wang')) then
		uwang=of()
	       end if
              else 
	       nstru=geti('#str',nstru)
	       if (find('pick')) then
		call pick(ipick,npick)
		do 11 i=1,npt
		 if (i1.eq.0 .and. ipick(i).ne.0) then
		  i1=i
		 else if (i2.eq.0 .and. ipick(i).ne.0) then
		  i2=i
		 else if (i3.eq.0 .and. ipick(i).ne.0)then
		  i3=i
	         end if
11		continue
               end if
	       if (find ('action')) goto  3
             end if
	     go to 1
3            continue


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
	   do 9 l=1,nstru
	    if (.not.norew) rewind urcrd
	    call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
            vec1x = coor(1,i1) - coor(1,i2)
            vec1y = coor(2,i1) - coor(2,i2)
            vec1z = coor(3,i1) - coor(3,i2)
            vec2x = coor(1,i3) - coor(1,i2)
            vec2y = coor(2,i3) - coor(2,i2)
            vec2z = coor(3,i3) - coor(3,i2)
	    norm1 = dsqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z)
	    norm2 = dsqrt(vec2x*vec2x + vec2y*vec2y + vec2z*vec2z)
	    cost  = (vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/(norm1*norm2)
	    if (cost.gt.  1.d0) cost = 1.d0
	    if (cost.lt. -1.d0) cost = -1.d0
	    theta = acos(cost)*180.d0/pi
	    write(*,*)' vec1 = ',vec1x,vec1y,vec1z
	    write(*,*)' vec2 = ',vec2x,vec2y,vec2z
	    write(*,*)l,theta,cost
	    write(uwang,*)l,theta
9	   continue

	   stop
	   end
