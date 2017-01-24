	subroutine dump_dyna(stdo)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/DYNA.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/RESTART.BLOCK'
	include 'COMMON/EWALD.BLOCK'

	integer stdo

	write(stdo,1)
1	format(//,' *** DYNAMICS PARAMETERS ')

	write(stdo,2)start_dyna,nstep,neqstep,ninfo
	write(stdo,25)ncoor,nvelo,nlist,
     1		nscalv,ntemp,newv,nrigi
2	format(1x,' Starting step for dynamics ',i7,/
     1 ,1x,' Number of dynamics steps                    ',i7,/
     1 ,1x,' Number of equilibration steps               ',i7,/
     2 ,1x,' Num. of steps between printing information  ',i7,/)
25	format(1x,' Num. of steps between writing coordinates   ',i7,/
     4 ,1x,' Num. of steps between writing velocities    ',i7,/
     5 ,1x,' Num. of steps between updates of nbnd lists ',i7,/
     6 ,1x,' Num. of steps between velocity scaling      ',i7,/
     7 ,1x,' Num. of temperatures in the system          ',i7,/
     8 ,1x,' Num of step between velocity reassignment   ',i7,/
     9 ,1x,' period for removal of rigid body motions    ',i7)

	write(stdo,3)boltz,shakm,shakl,shakb,nori,symanneal,freeze,eqms
3	format(1x,' FLAGS FOR DYNAMICS (T/F) ',/
     1	,1x,' Sample velocities at random from Boltzmann ',l1,/
     2	,1x,' Special shake for TIP3 water molecules     ',l1,/
     3	,1x,' Shake light particles only                 ',l1,/
     4	,1x,' Shake all bonds                            ',l1,/
     5	,1x,' **  DO NOT ** fix rigid body motion        ',l1,/
     6	,1x,' Change size of period box during annealing ',l1,/
     7	,1x,' Frozen particles are present               ',l1,/
     8	,1x,' Same mass for all (for annealing)          ',l1,//)

	write(stdo,4)itershak,epshak,epshakv
4	format(1x,' SHAKE PARAMETERS: ',/
     1	,1x,' Number of iterations = ',i6,/
     2	,1x,' Converegnce criterion for coordinates ',e10.5,/
     3  ,1x,' Convergence criterion for velocities  ',e10.5,//)
c write the initial coordinates and energies if it is the first run
c	if (.not. cont) then
c insight stuff
c		if ((ustbl .ne. -1) .or. (unstbl .ne. -1))
c     *			  call init_tbl()
c	 	if (esymyes) call squeeze(a,b,c)
c	 	call nbondm()
c		if (esymyes) call syminit()
c		if (specl) call nbondm_spcl()
c	        call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
c     1		curtemp,ntemp)
c		call eforce()
c        if (.not.esymyes) call ovrlpck(coor2,coor,dpot,jpick,iorie,rms)
c	if (.not.esymyes) call prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)
c     	  call info(stdo,0,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
c	  call wdyncrd(uwcrd,nstep/ncoor+1,0,inofrz,nofreez)
cc		if ((ustbl .ne. -1) .or. (unstbl .ne. -1))
c     *			 call w_tbl(0)	
c	end if
c

	return
	end
    
