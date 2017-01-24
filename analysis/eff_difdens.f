           program difdens
           implicit none
c
c calculate water properties from MD trajectories:
c effective diffusion coefficient and density
c for the definition of effective diffusion coefficient see
c Lounnas, Pettitt, and Phillips (1994) Biophys. J. 66:601
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
c input parameters specific for diffusion and density calculation:
c lbox   - edge of the cube used as histogram bin
c nstru0 - first structure in the dynamics file
c nstru  - last structure in the dynamics file
c tau1   - offset time (usually 1ps)
c dens0  - upper value for density
c ddens  - density increment (both dens values used to sort/label waters in PDB file)
c nrmono - number of solute monomers
c nratom - number of solute atoms (both nr used to generate PDB file)

           integer ipick(maxpt),ipikpo(maxpt),npick
	   integer urcrd,ucon
	   integer of,geti,nstru0,nstru
           integer rbin
	   integer namel,i,j,k,l,level
	   integer l1,l2,error
           integer i1,i2,i3,i4,nw,itau1
	   integer nx,ny,nz,nxmax(3)
	   integer nrmono,nratom,sumcount,ndens
	   integer ncount(50,50,50),ncountframes(50,50,50)
	   integer ncwat(50,50,50),idwat(50,50,50,20)
	   integer kcounts
	   double precision dens,densnorm,getd
	   double precision effdiff(50,50,50),effdiffct,avdiff
	   double precision lbox,tau1,dens0,ddens,dens1,zero
	   double precision mincoor(3),maxcoor(3),xsimbox(3)
	   double precision taucoor(3,maxpt)
	   double precision tau1coor(3,maxpt)
	   double precision xlbox,ylbox,zlbox
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
	   name='difdens'
	   namel=7
c	open junk file for rline
	   jnkf=25
	   open(unit=jnkf,status='scratch')
c	default parameters
	   nstru0=1
	   nstru=1
           lbox = 0.d0

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
	    print*,'lbox = ',lbox,'; npick = ',npick
	    print*,'nstru0 = ',nstru0,'; nstru = ',nstru

c	find min and max of x, y, and z
c	read all structures for min/max and box calculation

        do k=1,3
        mincoor(k)=9999.
        maxcoor(k)=-9999.
        enddo

	do i1=nstru0,nstru

	rewind urcrd
	if (dyna) then
	   call rdyncrd(urcrd,i1,inofrz,nofreez,rbin)
c	   print*,'npt = ',npt,'; inofrz = ',inofrz
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

	do k=1,3
	nxmax(k)=int((maxcoor(k)-mincoor(k))/lbox+1.)
	enddo

	print*,'mincoor(k) = ',(mincoor(k),k=1,3)
	print*,'maxcoor(k) = ',(maxcoor(k),k=1,3)
	print*,'xsimbox(k) = ',(xsimbox(k),k=1,3)
	print*,'nxmax(k) = ',(nxmax(k),k=1,3)

c initialize variables
	do nx=1,nxmax(1)
	do ny=1,nxmax(2)
	do nz=1,nxmax(3)
	ncount(nx,ny,nz)=0
	ncountframes(nx,ny,nz)=0
	effdiff(nx,ny,nz)=0.
	enddo
	enddo
	enddo

c	find box indeces (nx,ny,nz) for particle i3 of structure i1
c	since the effdiff calculation is done at 1ps intervals, read
c	MD frames every 10 structures

	itau1 = int(tau1+0.5)
	do 8 i1=nstru0,nstru-itau1

	do nx=1,nxmax(1)
	do ny=1,nxmax(2)
	do nz=1,nxmax(3)
	ncwat(nx,ny,nz)=0
	   do nw=1,20
	   idwat(nx,ny,nz,nw)=0
	   enddo
	enddo
	enddo
	enddo

	rewind urcrd
	call rdyncrd(urcrd,i1,inofrz,nofreez,rbin)
	   do 9 i2=1,npick
	   i3=ipikpo(i2)

	   nx=int((coor(1,i3)-mincoor(1))/lbox+1.)
	   ny=int((coor(2,i3)-mincoor(2))/lbox+1.)
	   nz=int((coor(3,i3)-mincoor(3))/lbox+1.)

c	count waters to obtain density at histogram bin (nx,ny,nz)
	   ncount(nx,ny,nz)=ncount(nx,ny,nz)+1

c	save order, identity, and taucoor of water i3 at (nx,ny,nz)
	   ncwat(nx,ny,nz)=ncwat(nx,ny,nz)+1	
	   i4=ncwat(nx,ny,nz) 
	   idwat(nx,ny,nz,i4)=i3
	      do k=1,3
	      taucoor(k,i3)=coor(k,i3)
	      enddo

9	continue

c	once a water is localized in a histogram bin, calculate its coordinates
c	after 1ps; also check if waters escape the simulation box and
c	need to be translated by symmetry

        rewind urcrd
	call rdyncrd(urcrd,i1+10,inofrz,nofreez,rbin)
	   do 10 i2=1,npick
	   i3=ipikpo(i2)

	   do k=1,3
	   tau1coor(k,i3)=coor(k,i3)
		if (abs(tau1coor(k,i3)-taucoor(k,i3)) .gt. xsimbox(k)/2.) then
		   if (tau1coor(k,i3) .lt. taucoor(k,i3)) then
		   tau1coor(k,i3)=tau1coor(k,i3)+xsimbox(k)
		   else
		   tau1coor(k,i3)=tau1coor(k,i3)-xsimbox(k)
		   endif
		endif
	   enddo
10	   continue

c	for each water in a histogram bin cummulate terms in effdiff
	   do nx=1,nxmax(1)-1
	   do ny=1,nxmax(2)-1
	   do nz=1,nxmax(3)-1
	   i2=ncwat(nx,ny,nz)

	   if (i2 .ne. 0) then
	   ncountframes(nx,ny,nz)=ncountframes(nx,ny,nz)+1
	   do i3=1,i2
	   i4=idwat(nx,ny,nz,i3)
c	print*,'bin = ',nx,ny,nz
c	print*,'first picked water id = ',i4,'; which water at that bin = ',i2
	      do k=1,3
c	print*,'coordinate = ',k,'; effdiff b/f summation = ',effdiff(nx,ny,nz)
	      effdiff(nx,ny,nz)=effdiff(nx,ny,nz)+(tau1coor(k,i4)-
     *		taucoor(k,i4))**2/(i2*1.)
c	print*,'added to effdiff: tau1-tau = ',tau1coor(k,i4)-taucoor(k,i4)
c	print*,'effdiff a/f summation = ',effdiff(nx,ny,nz)
	      enddo
	   enddo
	   endif

	   enddo
	   enddo
	   enddo

8	continue

c	check total number of counted waters
	sumcount=0
	do nx=1,nxmax(1)
	do ny=1,nxmax(2)
	do nz=1,nxmax(3)
	sumcount=sumcount+ncount(nx,ny,nz)
	enddo
	enddo
	enddo
	print*,'CHECK: sumcount = ',sumcount*1./(nstru-itau1-nstru0+1)

c	normalization constant for bulk density
	densnorm=0.0333395*lbox**3

	avdiff=0.
	kcounts=0
	do nx=1,nxmax(1)-1
	do ny=1,nxmax(2)-1
	do nz=1,nxmax(3)-1

	xlbox=mincoor(1)+(2*nx-1)*lbox/2.
	ylbox=mincoor(2)+(2*ny-1)*lbox/2.
	zlbox=mincoor(3)+(2*nz-1)*lbox/2.

	if (ncountframes(nx,ny,nz) .gt. 0.05*(nstru-itau1-nstru0+1)) then
c	if (ncount(nx,ny,nz) .ne. 0) then

	nrmono=nrmono+1
	nratom=nratom+1

	dens=1.*ncount(nx,ny,nz)/densnorm/(nstru-itau1-nstru0+1)

c	factor of 10 to have SI units, i.e. 10**(-9) m**2/s
	effdiffct=effdiff(nx,ny,nz)/ncountframes(nx,ny,nz)/6.*10.

c	calculate average diffusion coefficient
	avdiff=avdiff+effdiffct

c	count nr of histogram bins filled with water
	kcounts=kcounts+1

	write (102,999) nratom,nrmono,xlbox,ylbox,zlbox,dens,effdiffct
	write (103,1000) nratom,nrmono,xlbox,ylbox,zlbox,dens

c	write results in PDB format to be displayed graphically
c	data at each histogram bin are entered as "TIP waters"
c	located at the center of the box; however, these waters
c	differ in according to their density: TIP1 are highly
c	probable at that location (e.g. 0.8-1.0), TIP2 are less
c	likely to be there (e.g. 0.6-0.8) etc.; moreover, to each
c	water, include "B-factors" that represent effective diffusion
c	coefficients
c	such properties coded by "TIP water" name and "B-factor" can
c	be further displayed graphically

	ndens=int(dens0/ddens)
	dens1=dens0
	do k=1,ndens
c	print*,k,dens1,dens1-ddens
	if ((dens .ge. (dens1-ddens)) .and. (dens .lt. dens1)) then
	write (104,1001) nratom,k,nrmono,xlbox,ylbox,zlbox,effdiffct
	endif
	dens1=dens1-ddens
	   if ((dens1-ddens) .lt. 0.e-16) then
	   zero=dens1-ddens
	   dens1=dens1-zero
	   endif
	enddo

        endif

	enddo
	enddo
	enddo
999	format (4hATOM,i7,2h  ,4hOH2 ,4hTIP ,i7,4h    ,3f8.3,
     *  6h  1.00,1h ,f5.2,2h  ,f5.2)
1000	format (4hATOM,i7,2h  ,4hOH2 ,4hTIP ,i7,4h    ,3f8.3,
     *  6h  1.00,1h ,f5.2)
1001	format (4hATOM,i7,2h  ,4hOH2 ,3hTIP,i1,i7,4h    ,3f8.3,
c     *  6h  1.00,1h ,f10.5)
     *  6h  1.00,1h ,f5.2)

        print*,'nr of histogram bins filled with waters = ',kcounts
        print*,'total nr of histogram bins = ',
     *  nxmax(1)*nxmax(2)*nxmax(3)

        print*,'sum diffusion coefficient = ',avdiff
        print*,'average diffusion coefficient = ',avdiff/kcounts*1.

	end
