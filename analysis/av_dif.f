           program avrgdif
c
c calculate water properties from MD trajectories:
c average diffusion coefficient
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
c nstru0 - first structure read from the dynamics file
c nstru  - last structure read from the dynamics file
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
           integer ipick(maxpt),ipikpo(maxpt),npick
	   integer urcrd,ucon
	   integer of,geti,nstru,nstru0
	   integer namel,i,j,k,l,level
	   integer l1,l2,error
           integer rbin
	   double precision getd
	   double precision lbox,tau1,dens0,ddens,dens1,zero
           double precision mincoor(3),maxcoor(3),xsimbox(3)
	   double precision taucoor(3),reftaucoor(3,maxpt)
	   double precision diffcoor,sum,coeff
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
           rbin=1
	   totmon=0
	   npt=0
	   name='difdens'
	   namel=7
c	open junk file for rline
	   jnkf=25
	   open(unit=jnkf,status='scratch')
c	default parameters
	   nstru0=1
	   nstru=1
	   lbox = 0.d0
	   tau1 = 0.d0
	   dens0 = 1.5d0
	   ddens = 0.1d0
	   nrmono = 0
	   nratom = 0

1          continue
           call rline(name,namel,stdi)
	   if (find('norw')) norew=.true.
	   if (find('file')) then
		if (find ('conn')) then  
		ucon=of()
c	get connectivity 
	        call rconn(ucon)
                end if
c	get coordinate file
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
c
c read input parameters:
c lbox=edge of the cube used as histogram bin
c #str=number of MD frames
c tau1=offset time (1ps is fine)
c dens0=upper limit for density
c ddens=increments for density (both used to sort/label waters in PDB file)
c nrmono=number of solute monomers
c nratom=number of solute atoms (both used to generate PDB file)
c
 	        lbox = getd('lbox',lbox)
                nstru0=geti('lpst',nstru0)
                nstru=geti('lpen',nstru)
	        tau1=getd('tau1',tau1)
	        dens0=getd('dens0',dens0)
	        ddens=getd('ddens',ddens)
	        nrmono=geti('nrmono',nrmono)
	        nratom=geti('nratom',nratom)
	    end if

	    if (find('pick')) then
c for water analysis, pick OH2 atoms
	    call pick(ipick,i)
		npick = 0
		do 12 i=1,npt
		if (ipick(i) .ne. 0) then
			npick = npick + 1
			ipikpo(npick) = i
		endif
12		continue
	    endif

	    if (find('acti')) go to 3
	    go to 1
3           continue

            if (.not. fopen(ucon)) then
	     level=1
	     call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd)) then
	     level=1
	     call alert(name,namel,'urcrd not opened',16,level)
	    end if
	    print*,'lbox = ',lbox,'; nstru = ',nstru,'; npick = ',npick

        coeff=0.

c       find min and max of x, y, and z
c       read all structures for min/max and box calculation

        do k=1,3
        mincoor(k)=9999.
        maxcoor(k)=-9999.
        enddo

        do i1=nstru0,nstru

        rewind urcrd
        if (dyna) then
           call rdyncrd(urcrd,i1,inofrz,nofreez,rbin)
c          print*,'npt = ',npt,'; inofrz = ',inofrz
        else if (path) then
           call rpath(urcrd,1)
        else
           level = 1
           call alert(name,namel,'File type unkn',14,level)
        end if

        do i2=1,npick
        i3=ipikpo(i2)

        do k=1,3
           if (coor(k,i3) .lt. mincoor(k)) then
           mincoor(k)=coor(k,i3)
           else
              if (coor(k,i3) .gt. maxcoor(k)) then
              maxcoor(k)=coor(k,i3)
              endif
           endif
        enddo

        enddo

        enddo

c calculate box dimensions from min/max along each direction
        do k=1,3
        xsimbox(k)=maxcoor(k)-mincoor(k)
        enddo

	if (dyna) then

	do i0=nstru0,nstru-tau1
	sum=0.

	   rewind urcrd
	   call rdyncrd(urcrd,i0,inofrz,nofreez,rbin)

	   do i1=1,npick
	   i2=ipikpo(i1)
	     do k=1,3
	     reftaucoor(k,i2)=coor(k,i2)
	     enddo
	   enddo

	   rewind urcrd
	   call rdyncrd(urcrd,i0+10,inofrz,nofreez,rbin)

           do i1=1,npick
           i2=ipikpo(i1)
           diffcoor=0.
	   
          do k=1,3
	  taucoor(k)=coor(k,i2)
       if (abs(taucoor(k)-reftaucoor(k,i2)) .gt. xsimbox(k)/2.) then
                   if (taucoor(k) .lt. reftaucoor(k,i2)) then
                   taucoor(k)=taucoor(k)+xsimbox(k)
                   else
                   taucoor(k)=taucoor(k)-xsimbox(k)
                   endif
          endif
          enddo

	  do k=1,3
	  diffcoor=diffcoor+(reftaucoor(k,i2)-taucoor(k))**2/6.*10.
	  enddo

c	if (i0 .eq. 1) then
c	  write (110,101) i1,i2,diffcoor
c	  write (111,101) i1,i2,sqrt(diffcoor)
c	endif

	  sum=sum+diffcoor

	  enddo

	write (112,102) i0/10.,sum/npick*1.
        coeff=coeff+sum

	enddo

	endif

        print*,'Average diffusion coefficient = '
        print*,coeff/(nstru-nstru0-tau1)/npick*1.

101	format(2i10,f15.5)
102	format(2f15.5)

	end
