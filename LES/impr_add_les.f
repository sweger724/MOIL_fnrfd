        subroutine impr_add_les(poicpy,n_cpy)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/TMP_CONNECT.BLOCK'

        integer poicpy(0:maxpt),n_cpy(maxpt)

c local
        integer i,j1,j2,j3,j4,j1beg,j2beg,j3beg,j4beg
        integer k1,k2,k3,k4,case,icase
        double precision weight

c Create temporary torsion list
        tnimp = 0
        do 440 i = 1,nimp
                j1 = iimp1(i)
                j2 = iimp2(i)
                j3 = iimp3(i)
                j4 = iimp4(i)
                j1beg = poicpy(j1-1)+1
                j2beg = poicpy(j2-1)+1
                j3beg = poicpy(j3-1)+1
                j4beg = poicpy(j4-1)+1
		call test_many_ifs(case,j1,j2,j3,j4,lesid,weight,n_cpy)
		if (case.eq.0) then
c
c torsion type: either 0-0-0-0 or 1-2-3-4 (lesid wise)
c
		 do 441 k1=j1beg,poicpy(j1)
		  do 441 k2=j2beg,poicpy(j2)
		   do 441 k3=j3beg,poicpy(j3)
		    do 441 k4=j4beg,poicpy(j4)
			call fill_improper(k1,k2,k3,k4,i,weight)
441              continue
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
                         call fill_improper(k1,k2,k3,k4,i,weight)
                                 k2 = k2 + 1
                                 k3 = k3 + 1
                                 k4 = k4 + 1
4415                            continue
			else if (icase.eq.1) then
c
c torsion type: 1-0-0-0 
c
                         do 4422 k1=j1beg,poicpy(j1)
                           k3 = j3beg
                           k4 = j4beg
                           do 4422 k2=j2beg,poicpy(j2)
                            call fill_improper(k1,k2,k3,k4,i,weight)
                            k3 = k3 + 1
                            k4 = k4 + 1
4422                      continue

 			else if (icase.eq.2)  then
c
c torsion type: 0-1-0-0
c
			  do 4424 k2=j2beg,poicpy(j2)
			   k3 = j3beg
			   k4 = j4beg
			   do 4424 k1=j1beg,poicpy(j1)
			    call fill_improper(k1,k2,k3,k4,i,weight)
			    k3 = k3 + 1
			    k4 = k4 + 1
4424			 continue
c
c torsion type: 0-0-1-0 
c
                        else if (icase.eq.3) then
			 do 4431 k3=j3beg,poicpy(j3)
			  k1 =j1beg
			  k2 =j2beg
			  do 4431 k4=j4beg,poicpy(j4)
			   call fill_improper(k1,k2,k3,k4,i,weight)
			   k1 = k1 + 1
			   k2 = k2 + 1
4431			 continue
c
c torsion type: 0-0-0-1 
c

			else if (icase.eq.4) then
			 k1 =j1beg-1
			 k2 =j2beg-1
			 do 443 k3=j3beg,poicpy(j3)
			  k1 = k1 + 1
			  k2 = k2 + 1
			  do 443 k4=j4beg,poicpy(j4)
			   call fill_improper(k1,k2,k3,k4,i,weight)
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
			    call fill_improper(k1,k2,k3,k4,i,weight)
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
			    call fill_improper(k1,k2,k3,k4,i,weight)
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
			    call fill_improper(k1,k2,k3,k4,i,weight)
			    k4 = k4 + 1
446			 continue
			end if
		else if(case/1000.eq.2) then
			icase = case - 2000
			if (icase.eq.12) then
c
c torsion type 0-0-1-1
c
			  k2 = j2beg
			  k4 = j4beg
			  do 447 k1=j1beg,poicpy(j1)
			   do 447 k3=j3beg,poicpy(j3)
			    call fill_improper(k1,k2,k3,k4,i,weight)
			    k2 = k2 + 1
			    k4 = k4 + 1
447			  continue
			 else if (icase.eq.13) then
c
c torsion type: 0-1-0-1  OR
c torsion type: 0-1-1-0
c
                           k3 = j3beg
                           k4 = j4beg
                           do 448 k1=j1beg,poicpy(j1)
                            do 448 k2=j2beg,poicpy(j2)
                             call fill_improper(k1,k2,k3,k4,i,weight)
                             k3 = k3 + 1
                             k4 = k4 + 1
448                        continue
			 end if
		end if
440     continue

	call copy_improper()

        return
        end

