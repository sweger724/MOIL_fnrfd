	subroutine put_mis2(success,iat)
	implicit none
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	logical success,found
	integer iat
	character*6 name
	integer namel
c
c maxb is the maxmimum number of bonds to i
c maxused in the maximum of torsion we can use to build 
c missing atoms n molecules (also the number of missing atoms)
	integer maxb,maxused
	parameter (maxb=6,maxused=100)
	integer ibond,jbond,ibond_l1,ibond_l2(maxb)
        integer ibond_l3(maxb,maxb)
        integer used(4,maxused)
	integer i,j,k,l,j_select,iused
        integer jat,kat,lat
	integer level
c
c network is for the layer of the iat atom
c no more than 4 bonds per atom are expected
c
        integer network(maxb,maxb,maxb)

	double precision rij,ril,rjk,rjl,rkl,cijk,cijl,ckjl
        double precision dij,djk,djl
	double precision cjkl,cpijkl

        logical lused

        save used,lused,iused
        data used,iused,lused/400*0,0,.false./

	ibond = 0
c generating connectivity network for iat
c
        jat = 0
        kat = 0
	lat = 0
        ibond_l1 = 0
        do i = 1,maxb
         ibond_l2(i) = 0
         do j = 1,maxb
          ibond_l3(i,j)=0
          do k = 1,maxb
             network(i,j,k)=0
           end do
          end do
         end do

	name  = 'putmis'
	namel = 6
c
c How many bonds iat has?
c
	do i=1,nb
		if (ib1(i).eq.iat) then
	         if (coor(1,ib2(i)).lt.9998.) then
		  ibond_l1 = ibond_l1 + 1
		  network(ibond_l1,1,1) = ib2(i)
                 end if
		else if (ib2(i).eq.iat) then
	         if (coor(1,ib1(i)).lt.9998.) then
		  ibond_l1 = ibond_l1 + 1
		  network(ibond_l1,1,1) = ib1(i)
		 end if
		end if
	end do

c
c check now the number of bonds of the second layer
c
	if (ibond_l1.lt.1) then
	 level = 3
	 call alert(name,namel,'No bond to undefined crd',24,level)
	end if
c
        do i=1,ibond_l1
         do j = 1,nb
         if (ib1(j).eq.network(i,1,1)) then
          if (coor(1,ib2(j)).lt.9998.) then
           if (iat.ne.ib2(j)) then
            ibond_l2(i) = ibond_l2(i) + 1
            network(i,ibond_l2(i)+1,1) = ib2(j)
           end if
          end if
         else if (ib2(j).eq.network(i,1,1)) then
          if (coor(1,ib1(j)).lt.9998.) then
            if (iat.ne.ib1(j)) then
             ibond_l2(i) = ibond_l2(i)+1
             network(i,ibond_l2(i)+1,1) = ib1(j)
            end if
           end if
          end if
         end do
        end do
c
c check if the current network can be used to construct an atom with case 1
c In case the missig atom is bonded to j which in turn is bonded to atoms
c k and l. The coordinates of j,k, andl must be defined.
c
        do i=1,ibond_l1
         if (ibond_l2(i).ge.2) then
c
c can use case 1.
c
          jat = network(i,1,1)
          kat = network(i,2,1)
          lat = network(i,3,1)
          go to 11
         end if
        end do
c If this point is reached case 1 failed. Try case 2.
        go to 100
11	continue
c
c verify we have three
c if any of the added atoms is zero, something is missing
c
	if (jat*kat*lat.eq.0) then
		level = 0
                write(*,*)' i j k l ',iat,jat,kat,lat
		call alert(name,namel,'case 1 not option',17,level)
                go to 100
        else if (kat.eq.lat .or. kat.eq.iat .or. lat.eq.iat) then
                level=0
                write(*,*)' i j k l ',iat,jat,kat,lat
                call alert(name,namel,'atoms repeat',12,level)
                go to 100
        end if
c
c get equilibrium bond paramaters
c
	do i=1,nb
         if( (iat.eq.ib1(i).and.jat.eq.ib2(i))
     1     .or. (iat.eq.ib2(i).and.jat.eq.ib1(i))) then
           dij = req(i)
         else if( (jat.eq.ib1(i).and.kat.eq.ib2(i))
     1     .or. (jat.eq.ib2(i).and.kat.eq.ib1(i))) then
           djk = req(i)
         else if( (jat.eq.ib1(i).and.lat.eq.ib2(i))
     1     .or. (jat.eq.ib2(i).and.lat.eq.ib1(i))) then
           djl = req(i)
          end if
        end do

c
c get angle parameters
c We have i-j-k i-j-l & k-j-l
c
	do i=1,nangl
         if (jat.eq.iangl2(i)) then
          if (iat.eq.iangl1(i) .and. kat.eq.iangl3(i)) then
           cijk = cos(angleq(i))
          else if (iat.eq.iangl3(i) .and. kat.eq.iangl1(i)) then
           cijk = cos(angleq(i))
          else if (kat.eq.iangl1(i) .and. lat.eq.iangl3(i)) then
           ckjl = cos(angleq(i))
          else if (kat.eq.iangl3(i) .and. lat.eq.iangl1(i)) then
           ckjl = cos(angleq(i))
          else if (iat.eq.iangl1(i) .and. lat.eq.iangl3(i)) then
           cijl = cos(angleq(i))
          else if (iat.eq.iangl3(i) .and. lat.eq.iangl1(i)) then
           cijl = cos(angleq(i))
          end if
         end if
        end do
	call case1(coor(1,iat),coor(1,jat),coor(1,kat),
     1    coor(1,lat),dij,djk,djl,ckjl,cijk,cijl)
        success = .true.
        return
100     continue
c
c the beginning of case 2
            
        do i=1,ibond_l1
	 do j=1,ibond_l2(i)
          do k=1,nb
           if (ib1(k).eq.network(i,j+1,1)) then
            if (coor(1,ib2(k)).lt.9998) then
             if (iat.ne.ib2(k)) then
              ibond_l3(i,j)=ibond_l3(i,j) + 1
              network(i,j+1,ibond_l3(i,j)+1) = ib2(k)
             end if
            end if
           end if
           if (ib2(k).eq.network(i,j+1,1)) then
            if (coor(1,ib1(k)).lt.9998) then
             if (iat.ne.ib1(k)) then
              ibond_l3(i,j)=ibond_l3(i,j) + 1
              network(i,j+1,ibond_l3(i,j)+1) = ib1(k)
             end if
            end if
           end if
          end do
         end do
        end do


c
c analyze the network to find iat-j-k-l connection
c
        lused = .false.
	do i=1,ibond_l1
         if (network(i,1,1).ne.0) then
         do j=1,ibond_l2(i)
          if (network(i,j+1,1).ne.0) then
          do 1 k=1,ibond_l3(i,j)
           if (network(i,j+1,k+1).ne.0) then
            jat = network(i,1,1)
            kat = network(i,j+1,1)
            lat = network(i,j+1,k+1)
            do l=1,iused
             if ( jat.eq.used(2,l)
     1        .and. kat.eq.used(3,l) .and. lat.eq.used(4,l)) then
               lused = .true.
               go to 1  
             end if
            end do
            iused = iused + 1
            used (1,iused) = iat
            used (2,iused) = jat
            used (3,iused) = kat
            used (4,iused) = lat
            go to 3
           end if
1	  continue
          end if
         end do 
        end if
       end do
3      continue

c
c i--j--k--l
	do i=1,nangl
         if (jat.eq.iangl2(i)) then
          if (iat.eq.iangl1(i) .and. kat.eq.iangl3(i)) then
           cijk = cos(angleq(i))
          else if (iat.eq.iangl3(i) .and. kat.eq.iangl1(i)) then
           cijk = cos(angleq(i))
          end if
         end if
         if (kat.eq.iangl2(i)) then
          if (jat.eq.iangl1(i) .and. lat.eq.iangl3(i)) then
           cjkl = cos(angleq(i))
          else if (jat.eq.iangl3(i) .and. lat.eq.iangl1(i)) then
           cjkl = cos(angleq(i))
          end if
         end if
        end do
	        do j=1,nb
                 if ((ib1(j).eq.iat .and. ib2(j).eq.jat)
     1           .or. (ib1(j).eq.jat .and. ib2(j).eq.iat))
     2           rij = req(j)
                 if ((ib1(j).eq.jat .and. ib2(j).eq.kat)
     1           .or. (ib1(j).eq.kat .and. ib2(j).eq.jat))
     2           rjk = req(j)
                 if ((ib1(j).eq.kat .and. ib2(j).eq.lat)
     1           .or. (ib1(j).eq.lat .and. ib2(j).eq.kat))
     2           rkl = req(j)
		end do

		cpijkl = -9999.
		do j=1,ntors
		 if (itor2(j).eq.jat . and. itor3(j).eq.kat) then
			if (iat.eq.itor1(j) .and. lat.eq.itor4(j)) then
			 cpijkl = phase1(j)
			end if
                 else if (itor2(j).eq.kat .and. itor3(j).eq.jat) then
                        if (iat.eq.itor4(j) .and. l.eq.itor1(j)) then
			 cpijkl = phase1(j)
                        end if
		 end if
		end do

		if (cpijkl.lt.-9998.) then
		 level = 3
c@		 call alert(name,namel,'cp not determined',17,level)
                 if (.not.lused) then
                  cpijkl = -1
                 else
                  cpijkl = 0
                 end if
		end if

	call case2(coor(1,iat),coor(1,jat),coor(1,kat),coor(1,lat)
     1     ,rij,rjk,rkl,cijk,cjkl,cpijkl)
		success = .true.
		return	
		end
