	subroutine shakinit(shakl,shakb,shakm)
c
c NOTE THAT THIS ROUTINE WAS MODIED BY RE TO ALLOW THE SHAKE
c OF BONDS AND ANGLES IN STO.
c *** DO NOT USE IN REGULAR DYNAMICDS
c *** DO NOT REPLACE THIS ROUTINE BY REGULAR SHAKINIT FROM DYNA
c

c output: ishak1(nshak) ishak2(nshak) nshak dissq(nshak)
c         (in SHAKE.BLOCK)
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	logical shakl,shakb,shakm
c local
	character*8 name
	integer namel,i,i1,j1,level
	integer ibond1,ibond2
        double precision deltheta
	logical shak_accpt_bnd

	 if (nshak.gt.0) then
	  level = 0
	  call alert('shakinit',8,' Shakinit called twice',22,level)
	  return
	 end if
	if (shakb) then
	 	write(stdo,100)
100	 	format(/,1x,' All bonds will be shaked ',/)
	    else if (shakl) then
		 write(stdo,101)
101     format(/,1x,' All bonds with light particles will be shaked',/)
		else
		    level = 0
		    call alert(name,namel,
     *			           'routine called but not needed',29,level)
		    return
	end if
        if (shakm) write (*,*) 'shakm is on, not shaking TIP3'
	 do 1 i=1,nb
          if(shak_accpt_bnd(shakl,shakb,shakm,i)) then
             nshak = nshak + 1
             ishak1(nshak) = ib1(i)
             ishak2(nshak) = ib2(i)
             dissq(nshak)  = req(i)*req(i)
	     shak_diff(nshak) = dissq(nshak)*epshak
          endif
1        continue

         write (*,*) 'this gives a total of ',nshak,
     &   ' actual bond constraints'


c
c Now shaking the angles 
c
c @ try first without the angles
	return
c @ try WITH angles again.
c
         do 15 i=1,nangl
          if ((.not.shakm).or.(moname(poimon(iangl1(i)))
     1		.ne.'TIP3')) then
	   nshak     = nshak + 1 
	   i1        = iangl1(i)
	   j1        = iangl3(i) 
	   ishak1(nshak) = i1
	   ishak2(nshak) = j1
	   call find_bond(iangl2(i),i1,ibond1)
	   call find_bond(iangl2(i),j1,ibond2)
c	   deltheta =  (req(ibond2)-req(ibond1)*dcos(angleq(i)))
c     1		/ (req(ibond1)*dsin(angleq(i))) 
c	   deltheta = datan(deltheta)
c	   dissq(nshak)  = req(ibond1)*req(ibond1)+req(ibond2)*req(ibond2)
c     1		+ req(ibond1)*req(ibond2)*dcos(angleq(i))
      dissq(nshak)  = req(ibond1)*req(ibond1)+req(ibond2)*req(ibond2)
     1		- req(ibond1)*req(ibond2)*dcos(angleq(i))
	   shak_diff(nshak) = dissq(nshak)*epshak
	  end if
15	 continue
	 write(stdo,*)' Number of shake constraints = ',nshak
	 if (nshak.gt.maxshak) then
	  level = 1
	  call alert(name,namel,'Max. # of constraints exceed',28,level)
	 end if
	 return
	end
	
	logical function shak_accpt_bnd(shakl,shakb,shakm,i)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	logical shakl,shakb,shakm
	integer i

c	if(zerofrz(ib1(i)).eq.0 .or. zerofrz(ib2(i)).eq.0) then
c		shak_accpt_bnd = .false.
c		return
c	end if
	if ((shakm) .and. (moname(poimon(ib1(i))).eq.'TIP3')) then
		shak_accpt_bnd = .false.
		return
	end if
	if (lesid(ib1(i)).ne.lesid(ib2(i))) then
		shak_accpt_bnd = .false.
		return
	end if
	shak_accpt_bnd = .true.
	return
	end


























