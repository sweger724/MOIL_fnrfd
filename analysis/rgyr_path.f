           program rgyr
c
c calculate radius of gyration as a function of time from dynamic
c coordinate set.
c
C This new version manages path files -alfredo
C
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
c         value of zero if particl not selected one if it is
c ipikpo - pointer to picked particles ipikpo(i) is the particle
c          number of the i-th selected particle
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
c dd     - temporary variable storing the distance square between 2 pt
c xdiff ydiff zdiff - difference in x,y,z between pt (used in distance calc.)
c name - name of program (character) = contact
c find - find a charcter in line (logical function)
c fopen - check if file is open (logical function
c pickpt - true if pick instruction was found in input
c
	   integer ipick(maxpt), ipikpo(maxpt)
	   integer npick 
	   integer of,geti,nstru
           integer rbin
	   integer namel,i,j,k,l,level
	   integer urcrd,ucon,uwxy
	   double precision getd

	   integer ncrd
	   double precision step,cmx,cmy,cmz,cmm,tmpx,tmpy,tmpz,rg
	   double precision time
	   character*7 name
           character*4 ctype
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd,uwxy/3*99/

           ctype='DYNA'

	   lpstr = 1
	   norew = .false.

           stdi=5
	   stdo=6
           rbin = 1

	   totmon=0
	   npt=0
	   name='rgyr'
	   namel=4
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    nstru= 1
	    ncrd = 1000
	    step = 0.001d0
	    pickpt=.false.

1           continue
            call rline(name,namel,stdi)
            if (find('coor')) then
               call get4c(ctype,empty)
               write (6,*) 'file type ',ctype
            end if
	    if (find('norw')) norew=.true.
	    if (find('file')) then
	       if (find ('conn')) then  
		ucon=of()
* get connectivity 
	        call rconn(ucon)
* get coordinate file
	       else if (find ('rcrd')) then
		 if (npt.eq.0) then
		  level = 1
		  call alert(name,namel,'Must read con file first',
     1			24,level)
		 end if
	         urcrd=of()
                 else if (find('wcln')) then
                uwxy=of()
	       end if
              else 
	       nstru=geti('#str',nstru)
	       ncrd =geti('#crd',ncrd)
	       step =getd('step',step)
	       if (find('pick')) then
		call pick(ipick,npick)
		pickpt=.true.
               end if
	       if (find ('action')) goto  3
             end if
	     go to 1
3            continue

c initialized the nofreez vector
c
	     inofrz=npt
	     do 4 i=1,npt
	   	nofreez(i)=i
4            continue

          if (.not. fopen(ucon)) then
	     level=1
	     call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd)) then
	     level=1
	     call alert(name,namel,'urcrd not opened',16,level)
	    end if
	     
c           if (.not.pickpt) then
c	    level = 1
c	    call alert(name,namel,'Missing pick cmnd',17,level)
c	   end if

	   npick = 0
	   k = 0
	   do 6 i=1,npt
	    if (ipick(i) .eq. 1) then
	     npick = npick + 1
	     ipikpo(npick) = i
            end if
6          continue
           npick=npt

   
* read dynamics structures
	   rewind urcrd
	   do 9 l=1,nstru
              if (ctype.eq.'DYNA') then
		if (.not.norew) rewind urcrd
	        call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
              else if (ctype.eq.'PATH') then

                 call rpath_seq(urcrd,l)
              end if
c
c calculate center of mass cmx,cmy,cmz
c
		write(stdo,100)l
100		format(1x,'Dynamics structure # ',i5,1x,' Collisions:')
		 cmx = 0.d0
		 cmy = 0.d0
		 cmz = 0.d0
		 cmm = 0.d0
		 do 7 i=1,npick
		  j   = i
		  cmx = cmx + ptms(j)*coor(1,j)
		  cmy = cmy + ptms(j)*coor(2,j)
		  cmz = cmz + ptms(j)*coor(3,j)
		  cmm = cmm + ptms(j)
7		 continue
		 cmm = 1.d0/cmm
		 cmx = cmx*cmm
		 cmy = cmy*cmm
		 cmz = cmz*cmm
c
c calculate rgyr
c
		rg = 0.d0
		do 8 i=1,npick
		 j=i
		 tmpx = coor(1,j) - cmx
		 tmpy = coor(2,j) - cmy
		 tmpz = coor(3,j) - cmz
		 rg = rg + ptms(i)*(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz)
8	        continue
		rg = rg*cmm
		rg = dsqrt(rg)
c		time = l*ncrd*step
	   	write(uwxy,102)l, rg
102	   	format(1x,i5,1x,f15.3)
9	   continue
	   stop
	   end
