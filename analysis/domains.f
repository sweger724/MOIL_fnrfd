           program contact
c
c computing domains (structural parts with no contacts)
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
c	   number of the i-th NOT-selected  particle
c npick  - number of picked particles
c of     - integer function for Opening a File, returned value is the
c          assigned unit number
c geti   - integer function to get integer value from a command line
c getd   - double precision function to get DP value from a line
c nstru  - number of structures in dynamics file
c namel  - length of program name
c inofrz - number of moving particles. used in reading the dynamics file
c level  - level of error found, level=0 only warning will be issued,
c          level=1 program stop
c urcrd,ucon - units of dynamics coord and connectivity files
c rcut   - cutoff distance to define a collision (rcut>distance)
c          input rcut, internally modified and used as rcut^2
c xdiff ydiff zdiff - difference in x,y,z between pt (used in distance calc.)
c name - name of program (character) = contact
c find - find a charcter in line (logical function)
c fopen - check if file is open (logical function
c pickpt - true if pick instruction was found in input
c
	   parameter(maxmono1=105)
	   integer cont(maxmono1,maxmono1)
	   double precision t(maxmono1,maxmono1)
	   double precision target(maxmono1,maxmono1)
	   double precision dd(maxmono1),e(maxmono1),alpha1
	   double precision geomx(maxmono1),geomy(maxmono1)
	   double precision geomz(maxmono1)
	   double precision dx,dy,dz,dist
	   double precision sum(maxmono1)
	   integer n(maxmono1)
	   integer nres,lres,nmax
	   integer l1,l2,error
	   integer npick 
	   integer of,geti,nstru
	   integer namel,i,j,k,l,level
	   integer urcrd,ucon
           integer rbin
	   double precision rcut
	   double precision getd
	   character*7 name
	   logical find,fopen
	   logical pickpt
	   logical path,dyna
	   data ucon,urcrd/2*99/


	   lpstr = 1
	   norew = .false.

	   path = .false.
	   dyna = .false.
           stdi=5
	   stdo=6
           rbin = 1
 
	   totmon=0
	   npt=0
	   name='domains'
	   namel=7
*  open junk file for rline
*
	    jnkf=25
	    open(unit=jnkf,status='scratch')
* default parameters
	    nstru=1
	    rcut = 6.5d0

1           continue
            call rline(name,namel,stdi)
	    if (find('norw')) norew=.true.
	    if (find('file')) then
	       if (find ('conn')) then  
		ucon=of()
* get connectivity 
	        call rconn(ucon)
               end if
* get coordinate file
	       if (find ('rcrd')) then
		 if (find('DYNA')) dyna=.true.
		 if (find('PATH')) path=.true.
		 urcrd = of()
		 if (npt.eq.0) then
		  level = 1
		  call alert(name,namel,'Must read con file first',
     1			24,level)
		 end if
		end if
              else 
	       rcut = getd('rcut',rcut)
	       nstru=geti('#str',nstru)
             end if
	     if (find('acti')) go to 3
	     go to 1
3            continue

          if (.not. fopen(ucon)) then
	     level=1
	     call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd)) then
	     level=1
	     call alert(name,namel,'urcrd not opened',16,level)
	    end if
	     

	   ncount=0
* read dynamics structures
	   rewind urcrd
	   do 8 l=1,nstru
		if (dyna) then
	         call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
		else if (path) then
		 call rpath(urcrd,l)
		else
	         level = 1
		 call alert(name,namel,'File type unkn',14,level)
		end if
		nres = maxmono1
		do 9 l1=1,maxmono1
		if (moname(l1).eq.'NTG0') then 
		ncount=ncount+1
			l2 = 4
			geomx(l1)=coor(1,l2)
			geomy(l1)=coor(2,l2)
			geomz(l1)=coor(3,l2)
			go to 9
		endif
		if (moname(l1)(1:3).ne.'HEM' .or.
     1			 moname(l1).ne.'TIP3') then
		  geomx(l1) = 0.d0
		  geomy(l1) = 0.d0
		  geomz(l1) = 0.d0
		if (moname(l1)(1:3).eq.'GLY') then
       		ncount=ncount+1
 		   do 91 l2=poipt(l1-1),poipt(l1)-1
			if (ptnm(l2)(1:2).eq.'CA') then
				geomx(l1)=coor(1,l2)
				geomy(l1)=coor(2,l2)
				geomz(l1)=coor(3,l2)
				go to 9
			end if
91		    continue
		endif
		if (moname(l1)(1:3).eq.'PRO') then
       		ncount=ncount+1
		   do 92 l2=poipt(l1-1),poipt(l1)-1
			if (ptnm(l2)(1:2).eq.'CB') then
				geomx(l1)=coor(1,l2)
				geomy(l1)=coor(2,l2)
				geomz(l1)=coor(3,l2)
				go to 9
			end if
92		    continue
		endif
       		ncount=ncount+1
		  do 93 l2=poipt(l1-1),poipt(l1)-1
		    if (ptnm(l2)(1:2).ne.'N ' .or.
     1			ptnm(l2)(1:2).ne.'H ' .or.
     2			ptnm(l2)(1:2).ne.'CA' .or.
     3			ptnm(l2)(1:2).ne.'C ' .or.
     4			ptnm(l2)(1:2).ne.'O ' ) then
		   geomx(l1) = geomx(l1) + coor(1,l2)
		   geomy(l1) = geomy(l1) + coor(2,l2)
		   geomz(l1) = geomz(l1) + coor(3,l2)
		  end if
93		  continue
		  lres = poipt(l1)-poipt(l1-1)-5
C		  write(*,*) ' l1  lres = ',l1, lres
		  geomx(l1)=geomx(l1)/lres
 		  geomy(l1)=geomy(l1)/lres
		  geomz(l1)=geomz(l1)/lres
		endif
9		continue

		write(*,*) ' rcut = ',rcut
		write(*,*) ' ncount = ',ncount
c
c compute the contact matrix
c
		do 100 l1=1,maxmono1-2
		 cont(l1,l1)   = 0
		 cont(l1,l1+1) = 0
		 cont(l1,l1+2) = 0
		 cont(l1+1,l1) = 0
		 cont(l1+2,l1) = 0
		 do 100 l2=l1+3,maxmono1
		   dx = geomx(l2)-geomx(l1)
		   dy = geomy(l2)-geomy(l1)
		   dz = geomz(l2)-geomz(l1)
		   dist = dsqrt(dx*dx+dy*dy+dz*dz)
		   if (dist .le. rcut) then
			write(*,*)' l1 l2 dist ',l1,l2,dist
		   end if
		   if (dist.le.rcut) then
		   	cont(l1,l2) = 1
		   	cont(l2,l1) = 1
		   else
			cont(l1,l2) = 0
			cont(l2,l1) = 0
		   end if
100		continue

		write(*,*)' cont '
		do 10 l1=1,maxmono1
			write(*,1000)(cont(l1,l2),l2=1,maxmono1)
			write(*,*)
1000		format(20(1x,i3))
10		continue

		nmax = 0
		do 102 l1=1,maxmono1
		 n(l1) = 0
		 do 101 l2=1,maxmono1
		  if (cont(l1,l2).eq.1) then
		   n(l1) = n(l1) + 1
		  end if
101		 continue
		 if (n(l1).gt.nmax) nmax = n(l1)
102		continue

		alpha1 = 1.d0/(nmax+1)
		do 103 l1=1,maxmono1
		 t(l1,l1)=1-n(l1)*alpha1
		 do 103 l2=l1+1,maxmono1
			if (cont(l1,l2).ne.0) then		  
			 t(l1,l2)=alpha1
			 t(l2,l1)=alpha1
			else
			 t(l1,l2)=0.d0
			 t(l2,l1)=0.d0
			endif
		 target(l1,l2)=t(l1,l2)
		 target(l2,l1)=t(l2,l1)
		 target(l1,l1)=t(l1,l1)
103		 continue
		call house(target,maxmono1,maxmono1,dd,e,error)

		write(*,*)' nstru = ',l
		write(*,*) 'Eigenvalues '
		write(*,1001)dd
1001		format(1x,6(f8.3))
		do 104 l2=1,maxmono1
		write(*,1002)l2,dd(l2)
1002		format(1x,i2,1x,f8.3)
		write(*,1003) (target(l1,l2),l1=1,maxmono1)
1003		format(4(1x,f8.3))
C the header to householder routine says that
C target(i=1,n;j) is the j-th eigenvector
 104		continue

8	continue
	stop
	end
c
c compute matrix to order 10
c
