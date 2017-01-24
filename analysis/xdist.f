           program xdist
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
           include 'COMMON/MUTA.BLOCK'
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
	   integer namel,i,j,l,level
           integer i1(maxpt),i2(maxpt),ni1,ni2
	   integer urcrd,udist,ucon
           integer rbin
           double precision CM1(3),CM2(3),totmass1,totmass2
	   double precision dx,dy,dz,r

	   character*6 name
	   character*4 crd_type
	   character*5 getchar
	   logical find,fopen
	   logical pickpt
           double precision getd
           double precision a,b,c
           logical symm
	   data ucon,urcrd,udist/3*99/


	   ni1    = 0
	   ni2    = 0
	   lpstr = 1
	   norew = .false.

           stdi=5
	   stdo=6
           rbin = 1

	   totmon=0
	   npt=0
	   name='xdist'
	   namel=5
	   crd_type='DYNA'
c  open junk file for rline
c
	    jnkf=25
	    open(unit=jnkf,status='scratch')
c default parameters
	    nstru= 1
	    pickpt=.false.

        muta = .false.
        symm = .false.

	call init_var()

1           continue
            call rline(name,namel,stdi)
	    if (find('norw')) norew=.true.
	    if (find('file')) then
	       if (find ('conn')) then  
                if (find('muta')) muta = .true.
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
	        crd_type = getchar('ctyp',crd_type,4)
	       else if (find('wdis')) then
		udist=of()
	       end if
              else 
	       nstru=geti('#str',nstru)
		write(*,*)' nstru = ',nstru
	       if (find('pick')) then
		call pick(ipick,npick)
		do 11 i=1,npt
		 if (ipick(i).eq.1) then
		  ni1=ni1+1
                  i1(ni1) = i
                  write (*,*) 'picked type 1',i,ptnm(i),poimon(i),
     1                                        moname(poimon(i))
		 else if (ipick(i).eq.2) then
		  ni2=ni2+1
                  i2(ni2)=i
                  write (*,*) 'picked type 2',i,ptnm(i),poimon(i),
     1                                        moname(poimon(i))
	         end if
11		continue
               end if
               if (find('symm')) then
                symm = .true.
                a = getd('xtra',0.d0)
                b = getd('ytra',0.d0)
                c = getd('ztra',0.d0)
               endif
		write(*,*)' nstru = ',nstru
	       if (find ('action')) goto  3
             end if
	     go to 1
3            continue


		write(*,*)' pos 2 nstru = ',nstru
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
		write(*,*)' pos 3 nstru = ',nstru
	   do 9 l=1,nstru
            do j=1,3
             CM1(j) = 0.d0
             CM2(j) = 0.d0
            enddo
            totmass1=0.d0
            totmass2=0.d0
	    if (crd_type.eq.'DYNA') then
	    	if (.not.norew) rewind urcrd
	    	call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
	    else if (crd_type.eq.'PATH') then
		call rpath_seq(urcrd,l)
            else if (crd_type.eq.'CHAR') then
                call getcrd(urcrd,'CHARM')
	    else 
		write(*,*) ' Ctyp = ',crd_type(1:4)
		call alert(name,namel,'Unrecognized ctyp',17,level)
	    end if
            do i=1,ni1
             totmass1 = totmass1 + ptms(i1(i))
c correct for CM
             if (dabs(coor(1,i1(i))-coor(1,i1(1))).gt.0.5d0*a) then
              if (coor(1,i1(i)).gt.coor(1,i1(1))) then
               coor(1,i1(i)) = coor(1,i1(i)) - a
              else
               coor(1,i1(i)) = coor(1,i1(i)) + a
              endif
             endif
             if (dabs(coor(2,i1(i))-coor(2,i1(1))).gt.0.5d0*b) then
              if (coor(2,i1(i)).gt.coor(2,i1(1))) then
               coor(2,i1(i)) = coor(2,i1(i)) - b
              else
               coor(2,i1(i)) = coor(2,i1(i)) + b
              endif
             endif             
             if (dabs(coor(3,i1(i))-coor(3,i1(1))).gt.0.5d0*c) then
              if (coor(3,i1(i)).gt.coor(3,i1(1))) then
               coor(3,i1(i)) = coor(3,i1(i)) - c
              else
               coor(3,i1(i)) = coor(3,i1(i)) + c
              endif
             endif

             do j=1,3
              CM1(j) = CM1(j) + ptms(i1(i))*coor(j,i1(i))
             enddo
            enddo
            do i=1,ni2
             totmass2 = totmass2 + ptms(i2(i))
c correct for CM
             if (dabs(coor(1,i2(i))-coor(1,i2(1))).gt.0.5d0*a) then
              if (coor(1,i2(i)).gt.coor(1,i2(1))) then
               coor(1,i2(i)) = coor(1,i2(i)) - a
              else
               coor(1,i2(i)) = coor(1,i2(i)) + a
              endif
             endif
             if (dabs(coor(2,i2(i))-coor(2,i2(1))).gt.0.5d0*b) then
              if (coor(2,i2(i)).gt.coor(2,i2(1))) then
               coor(2,i2(i)) = coor(2,i2(i)) - b
              else
               coor(2,i2(i)) = coor(2,i2(i)) + b
              endif
             endif
             if (dabs(coor(3,i2(i))-coor(3,i2(1))).gt.0.5d0*c) then
              if (coor(3,i2(i)).gt.coor(3,i2(1))) then
               coor(3,i2(i)) = coor(3,i2(i)) - c
              else
               coor(3,i2(i)) = coor(3,i2(i)) + c
              endif
             endif

             do j=1,3
              CM2(j) = CM2(j) + ptms(i2(i))*coor(j,i2(i))
             enddo
            enddo
            do j=1,3
             CM1(j)=CM1(j)/totmass1
             CM2(j)=CM2(j)/totmass2
            enddo
		 dx = CM1(1) - CM2(1)
		 dy = CM1(2) - CM2(2)
		 dz = CM1(3) - CM2(3)
                 if (symm) then
                  dx = dx - a*anint(dx/a)
                  dy = dy - b*anint(dy/b)
                  dz = dz - c*anint(dz/c)
                 endif
		 r  = dsqrt(dx*dx + dy*dy + dz*dz)
	    	write(udist,*)l,r
9	   continue

	   stop
	   end

