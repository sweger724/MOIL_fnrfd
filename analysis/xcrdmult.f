           program xcrd
c
c Xtract distance between carboxyl oxygen of amino acid i and 
C amide hydrogen of amino acid i+4. If distance is less than
C 2.6 count it as a hydrogen bonding.
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
	   integer namel,i,j,l,level,i1,i2,i3,i4
       integer urcrd,uwtor,ucon,nh(10001),j1,j2,nhl(lgrid,maxmono)
           integer rbin

	   double precision dx1,dy1,dz1
	   double precision ux

	   character*6 name
           character*4 ctype
	   logical find,fopen
	   logical pickpt
	   data ucon,urcrd,uwtor/3*99/

           ctype='DYNA'
	   i1    = 0
	   i2    = 0
	   i3    = 0
	   i4    = 0
	   lpstr = 1
	   norew = .false.

           stdi=5
	   stdo=6
           rbin = 1

	   totmon=0
	   npt=0
	   name='xtors'
	   namel=5
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
            if (find('coor')) then
               call get4c(ctype,empty)
               write (6,*) 'file type ',ctype
            end if
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
	       else if (find('wtor')) then
		uwtor=of()
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
		 else if (i4.eq.0 .and. ipick(i).ne.0)then
		  i4=i
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
	     
   
           do j1=1,nstru
           nh(j1) =0
           do j2=1,totmon-4
           nhl(j1,j2)=0
           enddo
           enddo
           do j1=1,totmon-4
c read dynamics structures
	   rewind urcrd
           if (j1.eq.1) then
           inj=1
           else
           inj=poipt(j1-1)+1
           endif
           do j2=inj,poipt(j1)
           if (ptnm(j2).eq.'O   ') then
           i1=j2
           do j=poipt(j1+3)+1,poipt(j1+4)
           if (ptnm(j).eq.'H   ') then
           i2=j
           write (6,*) 'my pics ',i1,i2

	   do 9 l=1,nstru
              if (ctype.eq.'DYNA') then
                 if (.not.norew) rewind urcrd
                 call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
              else if (ctype.eq.'PATH') then

                 call rpath_seq(urcrd,l)
              end if

	  	dx1 = coor(1,i2) - coor(1,i1)
	  	dy1 = coor(2,i2) - coor(2,i1)
	  	dz1 = coor(3,i2) - coor(3,i1)
	

        ux=dx1*dx1+dy1*dy1+dz1*dz1
        ux = dsqrt(ux)
        if (ux.le.2.6d0) then
        nh(l)=nh(l)+1
c        write(uwtor,*) l,j1
        nhl(l,j1)=1
        endif
c	    	write(uwtor,*)ux
9	   continue
         endif
         enddo
         endif
         enddo
         enddo
         do l=1,nstru
c         write(uwtor,*) '*******',l
           do j1=1,totmon-4
c         if (nhl(l,j1).eq.1) then
         write(uwtor,*) nhl(l,j1)
c         endif
         enddo
         enddo
         do l=1,nstru
c         write(uwtor,*) l,nh(l)
         enddo

	   stop
	   end

