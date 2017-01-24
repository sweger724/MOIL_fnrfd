c
c this program takes a box of water with solute and cuts off any waters
c that do not have thier oxygen inside the inputs box coordinates
c assumes cubic symmetry of the box dimensions
c
	integer i,npt,a(1000000),b(1000000),len,len2,isav,j,iwrit,iwat
	integer nptprot,nptwat
	real x(1000000),y(1000000),z(1000000),r2,box,watsav,vdw(10000)
	real xij,yij,zij,dist
	character*50 trash(2)
	character*5 c(1000000),d(1000000)
	character*40 name,name2
	character chn*1
	logical save(1000000)
	real overlap	
c
c overlap is a parameter to allow a certain amount of vdw overlap to
c get the right density that you want (given in angstroms)
c this allows waters to be a little up the repulsive vdw wall
c
	overlap=.2
	write (*,*) 'enter allowed vdw overlap'
	read (*,*) overlap
	write (*,*) 'overlap = ',overlap
c
c get crd file
c
	write (*,*) 'enter protein CRD filename (with extension)'
	read (*,*) name
	len=0
	do 1310 i=1,40
		if (name(i:i).ne.' ') then
			len=len+1
		endif
1310    continue
	open (unit=15,name=name(1:len),status='old')
c
c assumes 2 title lines!
c
	read (15,5) trash(1)
	write (*,*) trash(1)
	read (15,5) trash(2)
	write (*,*) trash(2)
	read (15,2) nptprot
	write (*,*) nptprot
	write (*,*) 'enter water CRD filename (with extension)'
	read (*,*) name
	write (*,*) name
	len=0
	do 1300 i=1,40
		if (name(i:i).ne.' ') then
			len=len+1
		endif
1300    continue
	open (unit=16,name=name(1:len),status='old')
	read (16,5) trash(1)
	write (*,*) trash(1)
	read (16,5) trash(2)
	write (*,*) trash(2)
5	format (a1)		
	read (16,2) nptwat 
	write (*,*) nptwat 
2	format (i7)
	npt=nptwat+nptprot
	write (*,*) npt, 'total particles'
	do 3 i=1,npt
	  if (i.le.nptprot) then
	 	 read (15,4) a(i),b(i),c(i),d(i),x(i),y(i),z(i)
	  else
	 	 read (16,4) a(i),b(i),c(i),d(i),x(i),y(i),z(i)
	  endif
4	  format (1x,i4,1x,i4,1x,a4,1x,a4,3(1x,f9.5))
c
c set vdwrad
c
		chn=d(i)(1:1)
		write (*,*) chn
		if (chn.eq.'H') then
			vdw(i)=0.1
		elseif (chn.eq.'O') then
			vdw(i)=1.76
		elseif (chn.eq.'N') then
			vdw(i)=1.82
		elseif (chn.eq.'C') then
			vdw(i)=2.1
		else
			vdw(i)=2.1
		endif
3	continue
	close (unit=15)
	close (unit=16)
	write (*,*) 'enter output CRD filename (with extension)'
	read (*,*) name2
	len2=0
	do 1400 i=1,40
		if (name2(i:i).ne.' ') then
			len2=len2+1
		endif
1400    continue
	open (unit=16,name=name2(1:len2),status='unk')
	write (*,*) 'enter new cubic box size desired'
	write (*,*) 'as the length of a side'
	read (*,*) box
	box=box/2.
	write (*,*) 'will save water molecules from ',-box
	write (*,*) 'to ',box, ' that do not overlap protein'
c
c calc which water should be saved
c initialize all saved
c, if out of box or overlap then discard
c
	iwat=0
	isav=0
	watsav=0.0
	do 11 i=1,npt
		save(i)=.true. 
11	continue
	do 10 i=1,npt
		if (c(i)(1:4).eq.'TIP3') then
c
c calc the coordinates, if out of box discard the whole water 
c
			  if (x(i).gt.box.or.x(i).lt.-box.or.
     &			    y(i).gt.box.or.y(i).lt.-box.or.
     &			      z(i).gt.box.or.z(i).lt.-box) then
c
c check to see if other particles are before or after it
c
c assume particle in water in order   OH2,H1,H2
c
				if (d(i)(1:1).eq.'O') then
			 write (*,*) 'discard ',i,' out of box'
			  write (*,*) 'discard ',i+1,' out of box'
			  write (*,*) 'discard ',i+2,' out of box'
			          save(i)=.false.
				  save(i+1)=.false.
				  save(i+2)=.false.
				elseif (d(i)(1:2).eq.'H1') then
			  write (*,*) 'discard ',i,' out of box'
			  write (*,*) 'discard ',i+1,' out of box'
			  write (*,*) 'discard ',i-1,' out of box'
			          save(i)=.false.
				  save(i+1)=.false.
				  save(i-1)=.false.
				else
			  write (*,*) 'discard ',i,' out of box'
			  write (*,*) 'discard ',i-1,' out of box'
			  write (*,*) 'discard ',i-2,' out of box'
			          save(i)=.false.
				  save(i-1)=.false.
				  save(i-2)=.false.
			      endif
		  	  endif		
		endif
10	continue
c
c now calc the distances for the particles that made it to here
c
c discard the whole water if any vdw overlaps with non-waters
c
	write (*,*) 'checking for vdw overlaps'
c
	do 150 i=1,npt
c
c only water on outer loop, only protein on inner
c
	  if (c(i)(1:4).ne.'TIP3') go to 150
c
c don't bother if we are discarding it already
c
	  if (.not.save(i)) go to 150

	  do 160 j=1,npt
	    if (c(j)(1:4).eq.'TIP3') go to 160
	    xij=x(i)-x(j)
	    yij=y(i)-y(j)
	    zij=z(i)-z(j)
	    r2=xij*xij+yij*yij+zij*zij
	    dist=sqrt(r2)
	    if (dist+overlap.lt.(vdw(i)+vdw(j))) then
		write (*,*) 'discard ',i, dist,vdw(i)+vdw(j)
		if (d(i)(1:1).eq.'O') then
	          save(i)=.false.
		  save(i+1)=.false.
		  save(i+2)=.false.
		elseif (d(i)(1:2).eq.'H1') then
	          save(i)=.false.
		  save(i+1)=.false.
		  save(i-1)=.false.
		else
	          save(i)=.false.
		  save(i-1)=.false.
		  save(i-2)=.false.
	        endif
		go to 150
	    endif			
160	  continue
150	continue
c	
c now count number of particles to save and number of waters
c
	do 70 i=1,npt
		if (save(i)) isav=isav+1
		if (save(i).and.c(i).eq.'TIP3') watsav=watsav+1./3.
70	continue
	write (16,30) trash(1)
	write (16,30) trash(2)
30	format (a50)
	write (16,40) isav
40	format (i7)
	iwrit=0
	do 60 i=1,npt
	  if (save(i)) then
		iwrit=iwrit+1
	  	write  (16,50) iwrit,b(i),c(i),d(i),x(i),y(i),z(i)
50		format (1x,i4,1x,i4,1x,a4,1x,a4,3(1x,f9.5))
	  endif
60	continue
	close (unit=16)
	write (*,*) 'saved ',watsav,' water molecules'
	end


