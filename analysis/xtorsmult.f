           program xtors
c
c Xtract psi and phi torsion angles along a trajectory
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
	   integer ipick(maxpt),nh(lgrid,maxmono)
	   integer npick 
	   integer of,geti,nstru
	   integer namel,i,j,l,level,i1,i2,i3,i4
	   integer urcrd,uwtor,ucon
           integer rbin
	   double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
	   double precision ux,uy,uz,vx,vy,vz,uv,uu,vv,phi
           double precision phi1(lgrid),phi2(lgrid)

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
           do j1=1,lgrid
           phi1(j1)=0.0d0
           phi2(j1)=0.0d0
           enddo
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
	     
   
c read dynamics structures
           do j2=1,21
           do j1=1,nstru
           nh(j1,j2) =0
           enddo
           enddo
           ico=0
           do j1=1,npt-4
	   rewind urcrd
c           if (ptnm(j1).eq.'N   '.and.
c     &     moname(poimon(j1)).eq.'ALA ') then
           if (ptnm(j1).eq.'N   ') then
           ico=ico+1
           i1=j1
           i2=j1+2
           i3=j1+3
           jb=poimon(j1)+1
           do j2=poipt(jb-1)+1,poipt(jb)
           if (ptnm(j2).eq.'N   ') then
           i4=j2
           endif
           enddo
           
	   j = 0
           write (6,*) 'my pics ',i1,i2,i3,i4

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
	
	  	dx2 = coor(1,i3) - coor(1,i2)
	  	dy2 = coor(2,i3) - coor(2,i2)
	  	dz2 = coor(3,i3) - coor(3,i2)
	
	  	dx3 = coor(1,i4) - coor(1,i3)
	  	dy3 = coor(2,i4) - coor(2,i3)
	  	dz3 = coor(3,i4) - coor(3,i3)

	  	ux  = dy1*dz2 - dz1*dy2
	  	uy  = dz1*dx2 - dx1*dz2
	  	uz  = dx1*dy2 - dy1*dx2

	  	vx  = dy2*dz3 - dz2*dy3
	  	vy  = dz2*dx3 - dx2*dz3
	  	vz  = dx2*dy3 - dy2*dx3

	  	uu  = (ux*ux+uy*uy+uz*uz)
	  	vv  = (vx*vx+vy*vy+vz*vz)
	  	uv  = (ux*vx+uy*vy+uz*vz)/dsqrt(uu*vv)
	
	  	phi = dacos(uv)*pi180

	  	dx1 = uy*vz - uz*vy
	  	dy1 = uz*vx - ux*vz
	  	dz1 = ux*vy - uy*vx

	  	if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) phi = - phi
c	  	phi = phi - 180.d0
	  	if (phi.lt.-180) phi = phi + 360
                if (phi.ge.-77.0d0.and.phi.le.-17.0) then
                nh(l,ico)=1
                endif
c                phi1(l)=phi1(l)+phi
	    	write(uwtor,*)l,phi
c            write(uwtor,*) nh(l,ico)
9	   continue
           else
           endif
           enddo
            do j1=1,nstru
c            write(uwtor,*)j1,phi1(j1)/ico
            enddo
            write(uwtor,*)

c read dynamics structures
           do j2=1,21
           do j1=1,nstru
           nh(j1,j2) =0
           enddo
           enddo
           ico=0
           do j1=4,npt-4
	   rewind urcrd
           if (ptnm(j1).eq.'N   ') then
           ico=ico+1
           i2=j1
           i3=j1+2
           i4=j1+3
           jb=poimon(j1)-1
           do j2=poipt(jb-1)+1,poipt(jb)
           if (ptnm(j2).eq.'C   ') then
           i1=j2
           endif
           enddo
	   j = 0
           write (6,*) 'my pics ',i1,i2,i3,i4
           do 99 l=1,nstru
              if (ctype.eq.'DYNA') then
                 if (.not.norew) rewind urcrd
                 call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
              else if (ctype.eq.'PATH') then

                 call rpath_seq(urcrd,l)
              end if

                dx1 = coor(1,i2) - coor(1,i1)
                dy1 = coor(2,i2) - coor(2,i1)
                dz1 = coor(3,i2) - coor(3,i1)

                dx2 = coor(1,i3) - coor(1,i2)
                dy2 = coor(2,i3) - coor(2,i2)
                dz2 = coor(3,i3) - coor(3,i2)

                dx3 = coor(1,i4) - coor(1,i3)
                dy3 = coor(2,i4) - coor(2,i3)
                dz3 = coor(3,i4) - coor(3,i3)

                ux  = dy1*dz2 - dz1*dy2
                uy  = dz1*dx2 - dx1*dz2
                uz  = dx1*dy2 - dy1*dx2

                vx  = dy2*dz3 - dz2*dy3
                vy  = dz2*dx3 - dx2*dz3
                vz  = dx2*dy3 - dy2*dx3

                uu  = (ux*ux+uy*uy+uz*uz)
                vv  = (vx*vx+vy*vy+vz*vz)
                uv  = (ux*vx+uy*vy+uz*vz)/dsqrt(uu*vv)

                phi = dacos(uv)*pi180

                dx1 = uy*vz - uz*vy
                dy1 = uz*vx - ux*vz
                dz1 = ux*vy - uy*vx

                if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) phi = - phi
c               phi = phi - 180.d0
                if (phi.lt.-180) phi = phi + 360
                if (phi.ge.-87.0d0.and.phi.le.-27.0) then
                nh(l,ico)=1
                endif
c                phi2(l)=phi
                phi2(l)=phi2(l)+phi
c               write(uwtor,*)phi1(l),phi2(l)
               write(uwtor,*)l,phi
c            write(uwtor,*) nh(l,ico)
99          continue
           else
           endif
           enddo
            do j1=1,nstru
c            write(uwtor,*)j1,phi2(j1)/ico
c            write(6,*) ico
            enddo

	   stop
	   end

