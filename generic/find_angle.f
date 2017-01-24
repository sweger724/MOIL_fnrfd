	subroutine find_angle(iat1,iat2,iat3,iangle)
c
c find bond index (ib) given two atomic indices iat1 & iat2
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'

	integer iat1,iat2,iat3,iangle

c local
	integer sort_iat(3),sort_iangl(3)
	integer i,n,itmp,flips
	logical same

	sort_iat(1) = iat1
	sort_iat(2) = iat2
	sort_iat(3) = iat3
c@
	write(*,*)' find_angle: iat1 iat2 iat3 ',iat1,iat2,iat3

10	continue
	flips = 0
	do 1 i=1,2
	 if (sort_iat(i).gt.sort_iat(i+1)) then
	  write(*,*)' before i i+1 sort_iat(i) sort_iat(i+1) ',
     1		i,i+1,sort_iat(i),sort_iat(i+1)
	  itmp = sort_iat(i)
	  sort_iat(i) = sort_iat(i+1) 
	  sort_iat(i+1) = itmp
	  write(*,*)' after i i+1 sort_iat(i) sort_iat(i+1) ',
     1		i,i+1,sort_iat(i),sort_iat(i+1)
	  flips = 1
	 end if
1	continue
	if (flips.ne.0) go to 10
c@
	write(*,*)' find_angle: after sort iat1 iat2 iat3 '
	write(*,*) (sort_iat(i),i=1,3)

	
	do 5 n=1,nangl
		sort_iangl(1) = iangl1(n)
		sort_iangl(2) = iangl2(n)
		sort_iangl(3) = iangl3(n)
2		continue
		flips = 0
		do 3 j=1,2
		 if (sort_iangl(j).gt.sort_iangl(j+1)) then
			itmp = sort_iangl(j)
			sort_iangl(j) = sort_iangl(j+1)
			sort_iangl(j+1) = itmp
			flips =  1
		 end if
3		continue
		if (flips.ne.0) go to 2
		same = .true.
		do 4 j=1,3
		 same = same .and. (sort_iat(j).eq.sort_iangl(j))
4		continue
		if (same) then
			iangle = n
			return
		end if
5	continue
	call alert('find_angle',9,'Angle not found',15,1)
	stop
	end	
