        subroutine tors_add_les(poicpy,n_cpy)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/TMP_CONNECT.BLOCK'

        integer poicpy(0:maxpt),n_cpy(maxpt)

c local
        integer i,j1,j2,j3,j4,j1beg,j2beg,j3beg,j4beg
        integer k1,k2,k3,k4,case,icase
        double precision weight

c Create temporary torsion list
	tntors = 0
        do 440 i = 1,ntors
                j1 = itor1(i)
                j2 = itor2(i)
                j3 = itor3(i)
                j4 = itor4(i)
                j1beg = poicpy(j1-1)+1
                j2beg = poicpy(j2-1)+1
                j3beg = poicpy(j3-1)+1
                j4beg = poicpy(j4-1)+1
		call test_many_ifs(case,j1,j2,j3,j4,lesid,weight,n_cpy,cplbl)
		if (case.eq.0) then
c
c torsion type: either 0-0-0-0 or 1-2-3-4 (lesid wise)
c
			do 441 k1=j1beg,poicpy(j1)
			 do 441 k2=j2beg,poicpy(j2)
			  do 441 k3=j3beg,poicpy(j3)
			   do 441 k4=j4beg,poicpy(j4)
				call fill_torsion(k1,k2,k3,k4,i,weight)
441                     continue
                else if (case/10.eq.1) then
			icase = case-10
c
c torsion type: 1-1-1-1
c
			if (icase.eq.0) then
				k2 = j2beg
				k3 = j3beg
				k4 = j4beg
				do 4415 k1=j1beg,poicpy(j1)
				 call fill_torsion(k1,k2,k3,k4,i,weight)
				 k2 = k2 + 1
				 k3 = k3 + 1
				 k4 = k4 + 1
4415				continue
			else if (icase.eq.1 .or. icase.eq.2)  then
c
c torsion type: either 1-0-0-0 or 0-1-0-0
c
			 k3 = j3beg
			 k4 = j4beg
			 if (k1.eq.k3) then
			 do 4421 k1=j1beg,poicpy(j1)
			   k3 = k3 + 1
			   k4 = k4 + 1
			   do 4421 k2=j2beg,poicpy(j2)
			    call fill_torsion(k1,k2,k3,k4,i,weight)
4421			 continue
			 else
			  do 442 k2=j2beg,poicpy(j2)
			    k3 = k3 + 1
			    k4 = k4 + 1
			    do 442 k1=j1beg,poicpy(j1)
			     call fill_torsion(k1,k2,k3,k4,i,weight)
442			 continue
			 end if
                        else if (icase.eq.3 .or. icase.eq.4) then
c
c torsion type: either 0-0-1-0 or 0-0-0-1
c
			 k1 =j1beg
			 k2 =j2beg
			 do 443 k3=j3beg,poicpy(j3)
			  do 443 k4=j4beg,poicpy(j4)
			   call fill_torsion(k1,k2,k3,k4,i,weight)
			   k1 = k1 + 1
			   k2 = k2 + 1
443			 continue
			end if
		else if (case/100.eq.2) then
			icase = case - 200
			if (icase.eq.12) then
c
c torsion type: 0-0-1-2
c 
			 k2 = j2beg
			 do 444 k1=j1beg,poicpy(j1)
			  do 444 k3=j3beg,poicpy(j3)
			   do 444 k4=j4beg,poicpy(j4)
			    call fill_torsion(k1,k2,k3,k4,i,weight)
			    k2 = k2 + 1
444			 continue
			else if (icase.eq.13 .or. icase.eq.23) then
c
c torsion type: 0-1-0-2  OR
c torsion type: 1-0-0-2
c
			 k3 = j3beg
			 do 445 k1=j1beg,poicpy(j1)
			  do 445 k2=j2beg,poicpy(j2)
			   do 445 k4=j4beg,poicpy(j4)
			    call fill_torsion(k1,k2,k3,k4,i,weight)
			    k3 = k3 + 1
445			 continue
			else if (icase.eq.14 .or. icase.eq.24
     1					.or. icase.eq.34) then
c
c torsion type: 0-1-2-0 OR
c torsion type: 1-0-2-0 OR
c torsion type: 1-2-0-0
c
			 k4 = j4beg
			 do 446 k1=j1beg,poicpy(j1)
			  do 446 k2=j2beg,poicpy(j2)
			   do 446 k3=j3beg,poicpy(j3)
			    call fill_torsion(k1,k2,k3,k4,i,weight)
			    k4 = k4 + 1
446			 continue
			end if
		else if(case/1000.eq.2) then
			icase = case - 2000
			if (icase.eq.12) then
c
c torsion type 0-0-1-1
c
			 k2 = j2beg-1
			  do 447 k1=j1beg,poicpy(j1)
			   k2 = k2 + 1
			   k4 = j4beg
			   do 447 k3=j3beg,poicpy(j3)
			    call fill_torsion(k1,k2,k3,k4,i,weight)
			    k4 = k4 + 1
447			  continue
			 else if (icase.eq.13) then
c
c torsion type: 0-1-0-1  
c
			 if (lesid(j2).ne.lesid(j3)) then
                         k3 = j3beg
                         do 448 k1=j1beg,poicpy(j1)
                          k4 = j4beg
                          k3 = k3 + 1
                          do 448 k2=j2beg,poicpy(j2)
                           call fill_torsion(k1,k2,k3,k4,i,weight)
                           k4 = k4 + 1
448                       continue

c
c torsion type: 0-1-1-0
c
			else 
			k4 = j4beg
			do 449 k1=j1beg,poicpy(j1)
			 k4 = k4 + 1
			 k3 = j3beg
			 do 449 k2=j2beg,poicpy(j2)
			   call fill_torsion(k1,k2,k3,k4,i,weight)
			   k3 = k3 + 1
449			continue
			end if
c no torsion if particles 0 in type 0-1-0-0 or in 0-0-1-0 have different cplbl

                        else if (case.eq.2300) then
                         write(*,*) 'no torsion between particles'
     &                    ,j1,j2,j3,j4,'they have different cplbl:'
     &                    ,cplbl(j1),cplbl(j2),cplbl(j3),cplbl(j4)

			 end if
		end if
440     continue

	call copy_torsions()

        return
        end

