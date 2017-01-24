	subroutine frz_eval()
c
c pick freezing atoms and re-assign energy calculations
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/MSHAKE.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/DYNA.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'
	include 'COMMON/TETHER.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'

	integer i,j,istart,iend

	freeze = .true.

	call pick(ipick,i)
	inofrz = 0
c
c if parallel is on, check which of the frozen particles
c belong to the present pe. re-devide the moving particles between the pe-s
c
	if (prll_on_off) then
	 do 1 i=1,npt
	  if (ipick(i).ne.0) then
	    inofrz     = inofrz + 1
	    zerofrz(i) = 1
	  else
	    zerofrz(i) = 0
	  end if
1	 continue
	 call load_balance(inofrz,my_pe,num_pes,istart,iend)
c
c npt_par is the current number of particles that is assigned to this pe
c now MINUS the forzen particles.
c
	 npt_par = iend - istart + 1
	 j = 0
	 do 2 i=1,npt
	  if (ipick(i).ne.0) then
		j = j + 1
		if (j.lt.istart) go to 2
		if (j.gt.iend)   go to 3
		prtc_pointer(j) = i
		nofreez(j)      = i
	  end if
2	 continue
3	 continue

	else
c if NO parallel
c
	  
	do 4 i=1,npt
	   if (ipick(i).ne.0) then
	    inofrz               = inofrz + 1
	    nofreez(inofrz)      = i
	    prtc_pointer(inofrz) = i
	    zerofrz(i)           = 1
	   else
	    zerofrz(i)           = 0
	   end if
4	  continue
	  npt_par = inofrz

	end if

c Make a loop on covalent energy terms and eliminate unnnecessary
c energy elements (of particles that are frozen
c

c BONDS
	if (nb.gt.0) then
	 j = 0
5	 continue
	 j = j + 1
	 if (zerofrz(ib1(j)).eq.0 .and. zerofrz(ib2(j)) .eq.0) then
	  call rm_elemi(ib1,j,nb)
	  call rm_elemi(ib2,j,nb)
	  call rm_elemd(kbond,j,nb)
	  call rm_elemd(req,j,nb)
	  nb = nb - 1
	  j  = j - 1
	 end if
	 if (nb.gt.0 .and. (j.lt.nb)) go to 5
	end if

C MORSE & REP
	if (nmb.gt.0) then
	 j = 0
6	 continue
	 j = j + 1
	 if (zerofrz(imb1(j)).eq.0 .and. zerofrz(imb2(j)) .eq.0) then
	  call rm_elemi(imb1,j,nmb)
	  call rm_elemi(imb2,j,nmb)
	  call rm_elemd(D,j,nmb)
	  call rm_elemd(rmeq,j,nmb)
	  call rm_elemd(alpha,j,nmb)
	  call rm_elemd(beta1,j,nmb)
	  call rm_elemd(arep,j,nmb)
	  call rm_elemd(brep,j,nmb)
	  nmb = nmb - 1
	  j  = j - 1
	 end if
	 if (nmb.gt.0 .and. (j.lt.nmb)) go to 6
	end if

C ANGLES
	if (nangl.gt.0) then
	 j = 0
7	 continue
	 j = j + 1
	 if (zerofrz(iangl1(j)).eq.0 .and. zerofrz(iangl2(j)).eq.0
     1		.and. zerofrz(iangl3(j)).eq.0 ) then
	  call rm_elemi(iangl1,j,nangl)
	  call rm_elemi(iangl2,j,nangl)
	  call rm_elemi(iangl3,j,nangl)
	  call rm_elemd(kangl,j,nangl)
	  call rm_elemd(angleq,j,nangl)
	  nangl = nangl - 1
	  j  = j - 1
	 end if
	 if (nangl.gt.0 .and. (j.lt.nangl)) go to 7
	end if		 

C TORSIONS
	if (ntors.gt.0) then
	 j = 0
8	 continue
	 j = j + 1
	 if (zerofrz(itor1(j)).eq.0 .and. zerofrz(itor2(j)).eq.0
     1	.and.zerofrz(itor3(j)).eq.0.and.zerofrz(itor4(j)).eq.0) then
	  call rm_elemi(itor1,j,ntors)
	  call rm_elemi(itor2,j,ntors)
	  call rm_elemi(itor3,j,ntors)
	  call rm_elemi(itor4,j,ntors)
	  call rm_elemi(period,j,ntors)
	  call rm_elemd(ktors1,j,ntors)
	  call rm_elemd(ktors2,j,ntors)
	  call rm_elemd(ktors3,j,ntors)
	  call rm_elemd(phase,j,ntors)
	  ntors = ntors - 1
	  j  = j - 1
	 end if
	 if (ntors.gt.0 .and. (j.lt.ntors)) go to 8
	end if	


C IMPROPER TORSIONS
	if (nimp.gt.0) then
	 j = 0
9	 continue
	 j = j + 1
	 if (zerofrz(iimp1(j)).eq.0 .and. zerofrz(iimp2(j)).eq.0
     1	.and.zerofrz(iimp3(j)).eq.0.and.zerofrz(iimp4(j)).eq.0) then
	  call rm_elemi(iimp1,j,nimp)
	  call rm_elemi(iimp2,j,nimp)
	  call rm_elemi(iimp3,j,nimp)
	  call rm_elemi(iimp4,j,nimp)
	  call rm_elemd(kimp,j,nimp)
	  call rm_elemd(impeq,j,nimp)
	  nimp = nimp - 1
	  j  = j - 1
	 end if
	 if (nimp.gt.0 .and. (j.lt.nimp)) go to 9
	end if		 

C NO POINT IN COMPRESSING THE EXCLUSION LIST. IN nbond THE FIRST
C CHECK ON ANY PARTICLE PAIR IS ON zerofrz...

C 1-4 LIST
	if (totspe.gt.0) then
	 j = 0
10	 continue
	 j = j + 1
	 if (zerofrz(spec1(j)).eq.0 .and. zerofrz(spec2(j)).eq.0) then
	  call rm_elemi(spec1,j,totspe)
	  call rm_elemi(spec2,j,totspe)
	  call rm_elemd(p14,3*j-2,3*totspe)
	  call rm_elemd(p14,3*j-2,3*totspe-1)
	  call rm_elemd(p14,3*j-2,3*totspe-2)
	  totspe = totspe - 1
	  j  = j - 1
	 end if
	 if (totspe.gt.0 .and. (j.lt.totspe)) go to 10
	end if	

	return
	end	 
