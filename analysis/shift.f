           program translate
c
c Shifting a selected fragment of a molecule along the X axis by SHIF ansgtrom
c coordinate set.
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
	   integer of,geti
	   integer namel,i,j,k,l,level
	   integer urcrd,ucon,uwcrd
	   double precision getd

	   integer ncrd
	   double precision shift
	   character*7 name
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd/2*99/


	   lpstr = 1
	   norew = .false.

           stdi=5
	   stdo=6
	   totmon=0
	   npt=0
	   name='shif'
	   namel=4
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    shift = 15.d0
	    pickpt=.false.

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
		 if (npt.eq.0) then
		  level = 1
		  call alert(name,namel,'Must read con file first',
     1			24,level)
		 end if
	         urcrd=of()
	       end if
* Open crd file for shifting.
               if (find ('shif')) then
		  write(*,*) ' find shif in file '
                  uwcrd=of()
               end if
              else 
	       shift =getd('shif',shift)
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
	     
           if (.not.pickpt) then
	    level = 1
	    call alert(name,namel,'Missing pick cmnd',17,level)
	   end if

	   npick = 0
	   k = 0
	   do 6 i=1,npt
	    if (ipick(i) .eq. 1) then
	     npick = npick + 1
	     ipikpo(npick) = i
            end if
6          continue

   
* read structures
	   rewind urcrd
	   call getcrd(urcrd,'CHARM')
* Shift the x coordinates of the selected particles
		 do 7 i=1,npick
		  j   = ipikpo(i)
		  coor(1,j) = coor(1,j) + shift
7		 continue
		write(*,*)' before writing coordinate '
	    call putcrd(uwcrd,'CHARM')
		write(*,*) ' after writing coordinates '
	    stop
	    end
