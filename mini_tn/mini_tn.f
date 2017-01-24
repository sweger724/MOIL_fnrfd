	program mini_tn
	
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'
	include 'COMMON/VARBLE.BLOCK'
	include 'COMMON/CCRD.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	include 'COMMON/EWALD.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/CONVERT.BLOCK'

C SGB MODIFICATIONS
	include 'COMMON/SGB.BLOCK'
C END SGB MODIFICATIONS
	character*7 name
	character*4 cstyl,wstyl
	double precision getd
	integer namel
        integer sgbboolint
	integer of
	integer geti
	logical find

	double precision tolf,tmp
	integer ucon,ucor_in,ucor_out,umini,mistep,nlist
	integer i,istru,ig,istep,efcall,level,ncall,j
	integer ipick(maxpt)
        integer rbin


c setp-up for the truncated newton optimization
c the minimization code was taken from netlib
c author: Stephen G. Nash, George Mason University
c this part will be removed later for more efficient space
c handling
	external eforce0
	double precision dfpred,etotal
	double precision xx(3*maxpt),g(3*maxpt)
	double precision w(istore*3*maxpt)
	integer lw
	integer ier,idim

	sgbboolint=0
	stdi = 5
	stdo = 6
	stderr = 0
        rbin = 1

	totmon = 0
	npt    = 0
	nb     = 0
	nangl  = 0
	ntors  = 0
	nimp   = 0
	ncnst  = 0
	nvar   = 0
	lestyp = 0
	efcall = 0
	cstyl  = 'CHAR'
	wstyl =  'CHAR'

	jnkf = 25
	open (unit=jnkf,status='scratch')

c initialized file units and debug
	name      = 'mini_tn'
	namel     = 7
	ucon      = 101
	ucor_in   = 102
	ucor_out  = 103
	umini     = 104
	debug     = .false.

c energy default parameters
	call init_ef()
	call init_var()
	nocut  = .true.

	lcent  = .false.

	lpstr  = 1
	norew  = .false.


c minimization default parameters
c tolf  - allowed error (TOLerance) in the Force Norm
c mistep - maximum number of MInimization STEPs
	tolf   = 0.d0
	mistep = 100
	nlist  = 20

1	continue

	call rline(name,namel,stdi)

	if (find('debu')) debug = .true.

	if (find('file')) then
	 if (find('rcon') .or. find('conn')) then
	    ucon = of()
	    call rconn(ucon)
	    inofrz = npt
	    do 31 i=1,inofrz
		zerofrz(i) = 1
31	    continue
            npt_par = npt 		
	 endif
	 if (find('rcrd')) ucor_in = of()
	 if (find('wcrd')) then
           ucor_out= of()
	   wstyl = 'CHAR'
         end if
	 if (find('wpth')) then
           ucor_out= of()
	   wstyl = 'PATH'
         end if
	 if (find('wmin')) umini   = of()
	end if

	 if (find('cent')) then
	   i=1
	   lcent = .true.
	   kcenter = getd('kcnt',10.d0)
	   xeq = getd('xeqm',0.d0)
	   yeq = getd('yeqm',0.d0)
	   zeq = getd('zeqm',0.d0)
	   call pick(ipick,i)
	  icenter= 0
	  do 35 i=1,npt
	   if (ipick(i).ne.0) then
	    icenter=icenter + 1
	    center(icenter) = i
	   end if
35       continue
         endif

        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
	cutmono2 = getd('cutm',cutmono2)
	rmax     = getd('rmax',rmax)

	eps    = (getd('epsi',(eps)))
	if (.not. ctrue) ctrue  = find('cdie')
	if (find('rdie')) ctrue = .false.

	if (find('cpth')) then
		cstyl = 'PATH'
		istru = geti('istr',0)
	else if (find('cdyn')) then
		cstyl = 'DYNA'
		istru = geti('istr',0)
	end if

	if (ebyes  .and. find('nobo')) ebyes   = .false.
	if (ethyes .and. find('noan')) ethyes  = .false.
	if (etoyes .and. find('noto')) etoyes  = .false.
	if (eimyes .and. find('noim')) eimyes  = .false.
	if (evdyes .and. find('novd')) evdyes  = .false.
	if (eelyes .and. find('noel')) eelyes  = .false.
	if (nocut  .and. find('cute')) nocut   = .false.
	if (ecnyes .or.  find('cnst')) ecnyes  = .true.
	if (shift  .or.  find('shif')) shift   = .true.
	if (find('symm')) then
	 esymyes = .true.
	 a = getd('xtra',0.0d0)
	 b = getd('ytra',0.0d0)
	 c = getd('ztra',0.0d0)
	end if

        if (find('hvdw')) hvdw0 = .false.

c To be added after call rline of the relevant procedure
c
c jmjm
c
c switch for Ewald summation of long range interactions
         if (find('ewald')) then
            ewaldyes = .true.
            dtol = getd('dtol',0.0d0)
            nfft1 = geti('grdx',0)
            nfft2 = geti('grdy',0)
            nfft3 = geti('grdz',0)
            sgridx = getd('sgdx',1.0d0)
            sgridy = getd('sgdy',1.0d0)
            sgridz = getd('sgdz',1.0d0)
            intrpord = geti('iord',4)
            write (stdo,*)
     &  'PME account for long range inter. will be performed'

c if one wants to use PME for vacuum calculations
c one needs a virtual box which is sufficiently large to effectively
c separate image boxes
c in such a case stop should be commented and proper sizes
c of the virtual box  xtra,ytra,ztra should be given

            if (.not.esymyes) then
               write (stdo,*)
     &  'There is no true periodic bound. cond. - symm is missing'
               stop
c              a = getd('xtra',50.0d0)
c              b = getd('ytra',50.0d0)
c              c = getd('ztra',50.0d0)
            end if
         end if
c jmjmend
c
c virtual particles: massless, point charges displaced from the motion
c                    centers (real particles) e.g. TIP4P, 3-site CO model
c
         if (find('vprt')) then
           vp_flag = .true.
           com_flag = .true.
           gcnt_flag = .false.
           if (find('gcnt')) then
              com_flag = .false.
              gcnt_flag = .true.
           end if
         end if
c end of virtual partcles - one may actually think of moving this to con
c


        sgbboolint = geti('sgbb',sgbboolint)
	tolf   = (getd('tolf',(tolf)))
	mistep = geti('mist',mistep)
	nlist  = geti('list',nlist)

	if (find('TORS')) then
	   ncnst = ncnst + 1
	   icnst1(ncnst) = geti('atm1',0)
	   icnst2(ncnst) = geti('atm2',0)
	   icnst3(ncnst) = geti('atm3',0)
	   icnst4(ncnst) = geti('atm4',0)
	   kcns(ncnst)   = (getd('kcns',0.d0))
	   cnseq(ncnst)  = (getd('cneq',-999.d0))/pi180
	   if (debug) then
	    write(stdo,*)' kcns cnseq ',kcns(ncnst),cnseq(ncnst)*pi180
	   end if
	   if (find('loop')) then
		nvar = nvar + 1
		var1(nvar) = getd('strt',0.d0)/pi180
		var2(nvar) = getd('stop',0.d0)/pi180
		varstep(nvar) = (var2(nvar)-var1(nvar))
     1		/(getd('step',1.d0)/pi180) + 1.5
	   end if
	end if

	if (find('acti')) go to 2
	go to 1
2	continue

        if (sgbboolint.eq.0) then 
           sgbbool=.false.
        else
           sgbbool=.true.
        end if

c rmax is maintained here for old input (with a single
c cutoff to work)
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if

        if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
        if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

	write(stdo,111)mistep,nlist,ucon,ucor_in,ucor_out,
     1		cutvdw2,cutele2,tolf
111	format(1x,' Minimization parameters ',/,
     1  ' Maximum number of minimization steps ',i6,/,
     2  ' Update non-bond list each ',i6, ' steps',/,
     3  ' Connectivity file to be read from unit ',i5,/,
     4  ' Initial coordinates to be read from unit ',i5,/,
     5  ' Final coordinates will be written to unit ',i5,/,
     6  ' The cutoff distance for van der Waals is ',f10.5,/,
     7  ' The cutoff distance for electrostatic is ',f10.5,/,
     8  ' the allowed gradient error is ',f10.5)

	if (shift) write(stdo,*)' Shift dielectric will be used'

        cutvdw2  = cutvdw2*cutvdw2
        cutvbig2 = cutvbig2*cutvbig2
        cutele2  = cutele2*cutele2
        cutebig2 = cutebig2*cutebig2

	if (cutmono2.lt.0) then
		cutmono2 = cutebig2*1.44
	else
		cutmono2 = cutmono2*cutmono2
	end if

	if (nvar.ne.0) then
	 do 25 i=1,nvar
	  var2(i) = (var2(i)-var1(i))/dfloat(varstep(i)-1)
25	 continue
	end if
	tolf   = tolf*tolf


	if (cstyl.eq.'CHAR') then
	 call getcrd(ucor_in,'CHARM')
	else if (cstyl.eq.'PATH') then
	 if (istru.gt.0) then
	    
	  call rpath(ucor_in,istru)
	 else
	  level = 1
	  call alert(name,namel,'Illegal structure number',24,level)
	 end if
	else if (cstyl.eq.'DYNA') then
	 if (istru.gt.0) then
	  call rdyncrd(ucor_in,istru,inofrz,nofreez,rbin)
	 else
	  level = 1
	  call alert(name,namel,'Illegal structure number',24,level)
	 end if
	end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp



c *** HERE START A SINGLE MINIMIZATION

	if (nvar.eq.0) then

	write(umini,100)
100	format(1x,' At step 0 ',/)
	call init_wre(umini)

	ier  = 0
	idim = 3*npt
	lw   = istore*idim
	dfpred = 0.001d0*3*npt

	do 10 istep=1,mistep,nlist
C SGB MODIFICATION
	   chainno = 1
	   alphac(chainno) = 1
	   if (sgbbool) then
	      write (6,*) 'resetting alpha'
	   end if
C END SGB MODIFICATION
		call nbondm()
		if (esymyes) call syminit()
		do 3 i=1,npt
		 xx(i)       = coor(1,i)
		 xx(i+npt)   = coor(2,i)
		 xx(i+2*npt) = coor(3,i)
3		continue
		if (debug) then
		 write(stdo,*) 'idim = ',idim
		 write(stdo,*)'tolf nlist dfpred ',tolf,nlist,dfpred
		end if
		ncall = nlist
		call tn(ier,idim,xx,etotal,g,w,lw,eforce0,efcall
     1		 ,ncall,tolf,dfpred)
		write(umini,101)efcall
101		format(1x,'After ',i7,' energy calls ')
		tmp = 0.d0
		do 4 ig=1,idim
		 tmp = tmp + g(ig)*g(ig)
4		continue
		tmp = dsqrt(tmp/idim)
		write(umini,102)tmp,dsqrt(tolf)
102		format(1x,' Current gradient norm ',f10.5,
     1' requested ',f10.5)
		if (tmp*tmp.lt.tolf) then
		 write(umini,103)
103		 format(1x,' Minimization converged ')
		 efcall = efcall - 1
		 call eforce0(idim,xx,etotal,g,efcall)
		 call wener(umini)
		 go to 11
		end if
c copy "best results" to current coordinate set
		do 5 i=1,npt
		  coor(1,i)=xx(i)       
		  coor(2,i)=xx(i+npt)   
		  coor(3,i)=xx(i+2*npt) 
5		continue
6		continue
		if (ier.ne.0 .and. ier.ne.2) then
		 write(umini,*) 'Problems in minimization... ier =',ier
		end if
c get "best value" of energy and printout
		efcall = efcall - 1
		call eforce0(idim,xx,etotal,g,efcall)
		call wener(umini)
10	continue
11	continue
	if (wstyl.eq.'CHAR') then
	   call putcrd(ucor_out,'CHARM')
	else if (wstyl.eq.'PATH') then
	   call wpath(ucor_out)
	end if

	write(stdo,104)
104	format(1x,/,1x,' *** Minimization completed ')

c***	HERE END A SINGLE MINIMIZATION
c***	AND MULTI-MINIMIZATION WITH VARIABLE CONSTRAINTS STARTS

	else
	
	do 13 i=1,nvar
	 loop(i) = 1
13	continue

c 14 is the main loop on the constraints variables

	write(stdo,105)nvar,(var1(i)*pi180,i=1,nvar)
105	format(1x,//,' MINIMIZATION WITH ',i7,' CONSTRAINTS',//,
     1'Current equilibrium Positions',10(2x,f9.4),//,
     2'Below are constraints value and total energy')

14	continue
	
	do 15 i=1,nvar
	 cnseq(i) = var1(i) + (loop(i)-1)*var2(i)
15	continue

c	do minimization
	
	efcall = 0
	ier    = 0
	idim   = 3*npt
	dfpred = 0.001d0*3*npt

	do 18 istep=1,mistep,nlist
		call nbondm()
		if (esymyes) call syminit()
		do 16 i=1,npt
		 xx(i)       = coor(1,i)
		 xx(i+npt)   = coor(2,i)
		 xx(i+2*npt) = coor(3,i)
16		continue
		ncall = nlist
		call tn(ier,idim,xx,etotal,g,w,lw,eforce0,efcall
     1		 ,ncall,tolf,dfpred)
		write(umini,101)efcall
		tmp = 0.d0
		do 17 ig=1,idim
		 tmp = tmp + g(ig)*g(ig)
17		continue
		tmp = dsqrt(tmp/idim)
		write(umini,102)tmp,dsqrt(tolf)
		if (tmp*tmp.lt.tolf) then
		 write(umini,103)
		 efcall = efcall - 1
		 call eforce0(idim,xx,etotal,g,efcall)
		 call wener(umini)
		 go to 19
		end if
		if (ier.ne.0 .and. ier.ne.131) then
		 write(umini,*) 'Problems in minimization... ier =',ier
		end if
		efcall = efcall - 1
		call eforce0(idim,xx,etotal,g,efcall)
		call wener(umini)
18	continue
19	continue
	call wpath(ucor_out)

	write(stdo,106)(cnseq(i)*pi180,i=1,nvar),e_total-e_cnst,tmp
106	format(1x,'tors> ',4(1x,f10.5))
	do 21 i=1,nvar-1
	 if (loop(i).eq.varstep(i)) then
	  do 20 j=1,i
	   loop(j) = 1
20	  continue
	  loop(i+1) = loop(i+1) + 1
	  go to 22
	 else
	  loop(i) = loop(i) + 1
	  go to 22
	 end if
21	continue
22	continue
	if (nvar.eq.1) then
	 loop(1) = loop(1) + 1
	end if
	if (loop(nvar).gt.varstep(nvar)) stop
	go to 14
	end if
	stop
	end
