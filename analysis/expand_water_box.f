         program more_water
c special program written in order to translate water boxes and create large boxes
c we make a cluster of 8 boxes, translate in the X+/Y+/Z+ directions by the length of the current box + 0.2
c to avoid clushes.
c No need for a connectivity file. It is assumed that the water box is in CHARM format
c including only TIP3 water molecules. The title is only two lines and the number of atoms that follow
c the title is correct.
c

	character*80 filename
	character*1 char1
	character*4 c1,c2,c3,c4
	character*4 cres,cat1,cat2,cat3
	real coor(3,1000000)
	real a,b,c
	real bfactor
	integer i,j,nat,jatom, jres
	integer iatom,iwat

	write(*,*)' enter (small) box edges a, b and c '
	read(*,*)a,b,c
	write(*,*) ' enter water file name (CRD) '
	read(*,1)filename
1	format(80a)
	open (unit=1,file=filename,status='old')
	rewind 1
	write(*,*)' Enter new file name '
	read(*,1)filename
	open (unit=2,file=filename,status='new')
	rewind 2
	write(2,100)
	write(2,100)
100	format('*')
	do i=1,2
		read(1,2)char1
2		format(1a)
		if (char1.ne.'*') then
			write(*,*) ' missing title '
			stop
		end if
	end do
	read(1,*)nat
	if (8*nat.gt.1000000) then
		write(*,*) ' starting box too large '
		write(*,*)8*nat, ' must be smaller than 100000 '
		stop
	end if
	k = 0
	do i=1,nat/3
	 k = k+1
	 read(1,3)iatom,iwat,c1,c2,(coor(j,k),j=1,3)
     1	  ,c3,c4,bfactor
	 k = k+1
	 read(1,3)iatom,iwat,c1,c2,(coor(j,k),j=1,3)
     1	  ,c3,c4,bfactor
	 k = k+1
	 read(1,3)iatom,iwat,c1,c2,(coor(j,k),j=1,3)
     1	  ,c3,c4,bfactor
3	 format(i7,i7,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)

c 		add x+ water

	 k = k+1
	 coor(1,k) = coor(1,k-3) + a + 0.2
	 coor(2,k) = coor(2,k-3)
	 coor(3,k) = coor(3,k-3)
	 k = k+1
	 coor(1,k) = coor(1,k-3) + a + 0.2
	 coor(2,k) = coor(2,k-3)
	 coor(3,k) = coor(3,k-3)
	 k = k+1
	 coor(1,k) = coor(1,k-3) + a + 0.2
	 coor(2,k) = coor(2,k-3)
	 coor(3,k) = coor(3,k-3)

c		add y+ water

	 k = k+1
	 coor(2,k) = coor(2,k-6) + b + 0.2
	 coor(1,k) = coor(1,k-6)
	 coor(3,k) = coor(3,k-6)
	 k = k+1
	 coor(2,k) = coor(2,k-6) + b + 0.2
	 coor(1,k) = coor(1,k-6)
	 coor(3,k) = coor(3,k-6)
	 k = k+1
	 coor(2,k) = coor(2,k-6) + b + 0.2
	 coor(1,k) = coor(1,k-6)
	 coor(3,k) = coor(3,k-6)

c		add z+ water

	 k = k +1
	 coor(3,k) = coor(3,k-9) + c + 0.2
	 coor(1,k) = coor(1,k-9)
	 coor(2,k) = coor(2,k-9)
	 k = k +1
	 coor(3,k) = coor(3,k-9) + c + 0.2
	 coor(1,k) = coor(1,k-9)
	 coor(2,k) = coor(2,k-9)
	 k = k +1
	 coor(3,k) = coor(3,k-9) + c + 0.2
	 coor(1,k) = coor(1,k-9)
	 coor(2,k) = coor(2,k-9)

c		add xy+ water

	 k = k +1
	 coor(1,k) = coor(1,k-12) + a + 0.2
	 coor(2,k) = coor(2,k-12) + b + 0.2
	 coor(3,k) = coor(3,k-12)
	 k = k +1
	 coor(1,k) = coor(1,k-12) + a + 0.2
	 coor(2,k) = coor(2,k-12) + b + 0.2
	 coor(3,k) = coor(3,k-12)
	 k = k +1
	 coor(1,k) = coor(1,k-12) + a + 0.2
	 coor(2,k) = coor(2,k-12) + b + 0.2
	 coor(3,k) = coor(3,k-12)

c		add xz+ water

	 k = k + 1
	 coor(1,k) = coor(1,k-15) + a + 0.2
	 coor(2,k) = coor(2,k-15) 
	 coor(3,k) = coor(3,k-15) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-15) + a + 0.2
	 coor(2,k) = coor(2,k-15) 
	 coor(3,k) = coor(3,k-15) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-15) + a + 0.2
	 coor(2,k) = coor(2,k-15) 
	 coor(3,k) = coor(3,k-15) + c + 0.2

c		add yz+ water

	 k = k + 1
	 coor(1,k) = coor(1,k-18)
	 coor(2,k) = coor(2,k-18) + b + 0.2 
	 coor(3,k) = coor(3,k-18) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-18)
	 coor(2,k) = coor(2,k-18) + b + 0.2 
	 coor(3,k) = coor(3,k-18) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-18)
	 coor(2,k) = coor(2,k-18) + b + 0.2 
	 coor(3,k) = coor(3,k-18) + c + 0.2

c		add xyz+ water

	 k = k + 1
	 coor(1,k) = coor(1,k-21) + a + 0.2
	 coor(2,k) = coor(2,k-21) + b + 0.2 
	 coor(3,k) = coor(3,k-21) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-21) + a + 0.2
	 coor(2,k) = coor(2,k-21) + b + 0.2 
	 coor(3,k) = coor(3,k-21) + c + 0.2
	 k = k + 1
	 coor(1,k) = coor(1,k-21) + a + 0.2
	 coor(2,k) = coor(2,k-21) + b + 0.2 
	 coor(3,k) = coor(3,k-21) + c + 0.2

	end do

	write(2,200)k
200	format(i6)
	cres = 'TIP3'
	cat1 = 'OH2 '
	cat2 = 'H1  '
	cat3 = 'H2  '
	i = 0
4	continue
	i = i + 1
	ires = (i-1)/3 + 1
	 write(2,3)i,ires,cres,cat1,(coor(j,i),j=1,3)
     1	  ,c3,c4,bfactor
	 i = i + 1
	 write(2,3)i,ires,cres,cat2,(coor(j,i),j=1,3)
     1	  ,c3,c4,bfactor
	 i = i + 1
	 write(2,3)i,ires,cres,cat3,(coor(j,i),j=1,3)
     1	  ,c3,c4,bfactor
	if (i.ne.k) go to 4

	stop
	end
