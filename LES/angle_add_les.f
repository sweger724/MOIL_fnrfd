	subroutine angle_add_les(poicpy,n_cpy)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	integer poicpy(0:maxpt),n_cpy(maxpt)

c local
	integer i,j1,j2,j3,j1beg,j2beg,j3beg
	integer k1,k2,k3
	double precision weight

c Create temporary angle list
        tnangl = 0
        do 440 i = 1,nangl
                j1 = iangl1(i)
                j2 = iangl2(i)
                j3 = iangl3(i)
                j1beg = poicpy(j1-1)+1
                j2beg = poicpy(j2-1)+1
                j3beg = poicpy(j3-1)+1
                if (((lesid(j1).eq.0) .and. (lesid(j2).eq.0)
     *			 .and. (lesid(j3) .eq. 0))
     *                  .or. ((lesid(j1).ne.lesid(j2)) .and.
     *			 (lesid(j2).ne.lesid(j3))
     *                   .and. (lesid(j1).ne.lesid(j3))) ) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2)*n_cpy(j3))
                        do 441 k1=j1beg,poicpy(j1)
                         do 441 k2=j2beg,poicpy(j2)
                          do 441 k3=j3beg,poicpy(j3)
				call fill_angle(k1,k2,k3,i,weight)
441                     continue
                else if ((lesid(j1).eq.lesid(j2)) .and.
     *				 (lesid(j1).ne.lesid(j3))) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j3))
                        k2 = j2beg -1
                        do 442 k1=j1beg,poicpy(j1)
                          k2 = k2 + 1
                          do 443 k3=j3beg,poicpy(j3)
				call fill_angle(k1,k2,k3,i,weight)
443                       continue
442                     continue
                 else if ((lesid(j2).eq.lesid(j3)) .and.
     *				 (lesid(j1).ne.lesid(j2))) then
                        weight = 1.d0/(n_cpy(j1)*n_cpy(j2))
                        do 447 k1=j1beg,poicpy(j1)
                         k2 = j2beg
                         do 448 k3=j3beg,poicpy(j3)
				call fill_angle(k1,k2,k3,i,weight)
                         	k2 = k2 + 1
448                      continue
447                     continue
                 else if (lesid(j1).eq.lesid(j2) .and.
     1		   lesid(j2).eq.lesid(j3) .and. lesid(j1).ne.0) then
                         weight = 1.d0/n_cpy(j1)
                         k2 = j2beg
                         k3 = j3beg
                         do 449 k1=j1beg,poicpy(j1)
                                call fill_angle(k1,k2,k3,i,weight)
                         	k2 = k2 + 1
                         	k3 = k3 + 1
449                      continue
                 end if
440     continue

	call copy_angles()

	return
	end
