           program xcrd
	   implicit none
c
c Xtract specific coordinate sets from a dynamic file
c
	   include 'COMMON/LENGTH.BLOCK'
	   include 'COMMON/COORD.BLOCK'
	   include 'COMMON/VELOC.BLOCK'
	   include 'COMMON/ENERGY.BLOCK'
	   include 'COMMON/UNITS.BLOCK'
	   include 'COMMON/CONNECT.BLOCK'
	   include 'COMMON/DEBUG.BLOCK'
	   include 'COMMON/LINE.BLOCK'
	   include 'COMMON/FREEZ.BLOCK'
	   include 'COMMON/CONVERT.BLOCK'
	   include 'COMMON/CCRD.BLOCK'
c ipick - pick subset of particles, a vector of length maxpt
c         value of zero if particle not selected one if it is
c ipikpo - pointer to picked particles ipikpo(i) is the particle
c          number of the i-th selected particle
c npick  - number of picked particles
c of     - integer function for Opening a File, returned value is the
c          assigned unit number
c geti   - integer function to get integer value from a command line
c nstru  - number of structures to be looked in dynamics file
c stru1  - first structure to ouput
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
	   double precision rms
	   integer ipick(maxpt), ipikpo(0:maxpt)
	   integer npick 
	   integer ipic2(maxpt), ipikpo2(maxpt)
	   integer npic2
	   integer num_mon
	   integer of,geti,nstru, stru1
	   integer namel,i,i1,j,k,l,level
	   integer imon,imon1
	   integer urcrd,uwcrd,ucon
           integer rbin

	   character*7 name
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd,uwcrd/3*99/


	   lpstr = 1
	   norew = .false.
	   ipikpo(0) = -1

           stdi=5
	   stdo=6
           rbin = 1

	   totmon=0
	   num_mon = 0
	   npt=0
	   name='xcrd'
	   namel=4
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    nstru= 1
          stru1=1
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
	       else if (find('wcrd')) then
		uwcrd=of()
	       end if
              else 
	       nstru=geti('#str',nstru)
               stru1=geti('str1',stru1)
	       if (find('pick')) then
		call pick(ipick,npick)
		pickpt=.true.
               end if
	       if (find('pic2')) then
		call pick(ipic2,npic2)
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
	     
           if (.not.pickpt) then
	    level = 1
	    call alert(name,namel,'Missing pick cmnd',17,level)
	   end if

c initialized the nofreez vector
c
	     inofrz=npt
	     do 4 i=1,npt
	   	nofreez(i)=i
4            continue
	     num_mon = totmon
	     call vinit(dpot,0.0d0,3*npt)
	     call vinit(velo,0.0d0,3*npt)

	   npick = 0
	   do 6 i=1,npt
	    if (ipick(i) .eq. 1) then
	     npick = npick + 1
	     ipikpo(npick) = i
            end if
	    if (ipic2(i) .eq. 1) then
	      npic2 = npic2 + 1
	      ipikpo2(npic2) = i
	    end if
6          continue

c nstru is the final structure now
         nstru=nstru+stru1-1
c read the first structure for overlaps
	   rewind urcrd
	   call rdyncrd(urcrd,1,inofrz,nofreez,rbin)
	   call vdcopy(coor,coor2,3*npt)
	   call putcrd(uwcrd,'CHARM')
c read dynamics structures
	   rewind urcrd
	   do 9 l=1,nstru
	     j = npt
	     if (.not.norew) rewind urcrd
	     call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
	     call ovrlpck(coor2,coor,dpot,velo,ipikpo2,npic2,rms)
           if (l.ge.stru1) then
	       do 7 k=1,npick
                 i1 = ipikpo(k-1)
	         i = ipikpo(k)
                 imon1 = poimon(i1)
		 imon  = poimon(i)
		 j = j+1
		 if (imon1.ne.imon) num_mon = num_mon + 1
	         write(uwcrd,103)j,num_mon,moname(poimon(i))
     1		,ptnm(i),coor(1,i),coor(2,i),coor(3,i),BULK(1),0.d0
103	         format(i7,i7,1x,a4,1x,a4,3(f10.5),1x,a4,1x,'FREE',f10.5)
7	       continue
           endif 
9	   continue
	   stop
	   end
