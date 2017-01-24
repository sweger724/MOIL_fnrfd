        subroutine info(unit,jstep,nofreez,inofrz,tpo,tgroup,
     1		ntemp,curtemp)
c
c integer unit
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	character*4 name
	integer unit,namel,jstep,ntemp,i
	integer inofrz,nofreez(*),tpo(*),tgroup(*)
	double precision curtemp(*)
	double precision all

	name  = 'info'
	namel = 4

        if (my_pe.eq.0) then
           write(unit,*)
           write(unit,*)
           write(unit,*)
           write(unit,*)'-----------------------------------'
        endif

	if (debug) then
	 write(unit,*)' vx  ',(velo(1,i),i=1,npt)
	 write(unit,*)' vy  ',(velo(2,i),i=1,npt)
	 write(unit,*)' vz  ',(velo(3,i),i=1,npt)
	 write(unit,*)' ptms ',(ptms(i),i=1,npt)
	end if

        if (prll_on_off.and.matshak) then
           call multemp_prll(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1		curtemp,ntemp)
        else
           call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     $          curtemp,ntemp)
        endif

	if (my_pe.eq.0) then
           write(unit,100)jstep
           write(unit,101)(curtemp(i),i=1,ntemp)
 100       format(1x,' At dynamics step ',i7)
 101       format(1x,' Current temperature(s) are  ',5(1x,f10.2))
	end if


	call wener(unit)

	all = 0.d0
	do 1 i=1,ntemp
	 all = all + 0.5d0*kboltzmann*curtemp(i)*tgroup(i)
1	continue
	all = e_total + all

	if (my_pe.eq.0) then
           write(unit,102)all
 102       format(1x,' current energy (kinetic+potential) is ',f15.3)

           write(unit,*)'-----------------------------------'
           write(unit,*)
           write(unit,*)
           write(unit,*)
        endif

	return 
	end
