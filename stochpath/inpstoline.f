      subroutine inpstoline
c
c initializing as many as possible of the DYNA and STO variables
c standard input lines
c Date of creation : October 27, 1999
c Author: V. Zaloj
c Last update : November 04, 1999
c Author: V. Zaloj
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/COMM_LOCAL.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/METAL.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/PVM_LOCAL.BLOCK'
      include 'COMMON/RESTART.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/VELOC.BLOCK'
c    
      integer geti,of,ierr
      double precision getd
      logical find,fopen
c     character*80 getchar
c     
c..   new local variables.
c.....name   - name of the subroutine that calls rline
c.....namel  - lenght of name
c.....dyna   - flag to direct the program to either simple minimization 
c     (.false.) or simulated annealing (.true.).
c.....sdpyes - flag to specify the kind of minimization routine. 
c     Actually not used meanwhile.
c.....ucon   - file unit number of the connectivity file (read)
c.....urcrd  - file unit number of the coordinate input
c.....middle - index of the middle structure (actually not used here, 
c     but used in the subsequent subroutines)
c.....istart - index of the grid for the first structure used by 
c     the present processor (it is not included in any common 
c     block since it appears also in getcrd as local variable)
c.....istru  - number of grid points used in this processor
c.....udata  - unit of the standard output passed to the subroutines 
c     sto1 and sto1_dyn
c.....prep   - label to activate the preprocessing (minimization) of 
c     the structures in initial path 
c     
      character*10 name
      integer namel
c     
      integer ucon,urcrd,umini,uwmini
      integer i,ratio,k,j,ii,kk
      integer istart,istru
c     
      logical sdpyes,dyna,noinert,prep
c     
      data ucon,urcrd,uwmini,umini/4*99/
c     
c     
c.....GENERAL INITIALIZATION OF LOCAL VARS
c     
c
      sdpyes = .false.
      dyna   = .false.
      noinert= .false.
      prep   = .false.
      ratio  = 1
      name   = 'inpstoline'
      namel  = 10

c
c.....initialize shake variables
c
	shakl   = .false.
	shakb   = .true.
	shakm   = .true.
	epshak  = 1.d-5
	epshakv = 1.d-5
	itershak = 1000
c     
c.....initialization of the parallel machine
c
      call init_paral(istart,istru)
c
c open scratch file for line manipulation (check if still used).
czva  Scratch file moved here to avoid problem with shared 
czva  memory MPI (Linux) 05/26/00
c.....LINE
cout      jnkf = 25
c.....(changed to 1 for the sp/2, running poe)
      jnkf = 1
c
      open (unit=jnkf,status='scratch')
c     
c.....more parameter initialization (energy default parameters)
czva  Initialization was done already in "init_sto" (same call) 05/26/00   
c      call init_ef()
c
c.....file initialization
      if (paral) call open_inout(stdi,stdo,uwcrd,my_pe)
c     
c.....Parameter input from stdi using the line interpreter routines.
c.....Each process reads its own input.
c.....loop of the kind "do until"
c...................initializa CPU .......

cccc      call printrtc(stdo)
c     
      call inittime
c     
 1    continue
      call rline(name,namel,stdi)
      if (find('debu')) debug = .true.
      if (find('file')) then
         if (find('conn')) then
cc            if (my_pe.eq.0 ) then
c..............only the first processor reads the connectivity data.
               ucon = of()
               call rconn(ucon)
ccc            endif

c$$$c.....subroutines of TERRA for broadcast the connectivity and some 
c$$$c.....energy parameters file.
c$$$            write(stdo,*)' Elapsed time Before calling broad_conn '
c$$$            call printcpu(stdo)
c$$$            call broad_conn()
c$$$            write(stdo,*)' Elapsed time After calling broad_conn '
c$$$            call printcpu(stdo)
c 
               npts=npt    
               npt2 = 2*npt
               npt3 = 3*npt
               inofrz = npt
               iorie  = npt
c	initialize nofreez,jpick and tpo and inverse mass
c check if masses are set to be equal (for easier annealing)
               do 3 i=1,npt
                  nofreez(i) = i
                  zerofrz(i) = 1
                  jpick(i)   = i
                  tpo(i)     = 1
                  if (ptms(i).lt.1.0d-9) then
                     invms(i)=0.d0
c this is to avoid division by zero mass of vp
                  else
                     invms(i) = 1.d0/ptms(i)
                  end if              
 3          continue
	    do 7 i=1,totdmon
	       frzM(i) = .false.
 7	    continue            
            tgroup(1) = 3*npt
         end if 
         if (find('rcrd').and. my_pe.eq.0) urcrd = of()
         if (find('wpth').and. my_pe.eq.0) uwpth = of()
         if (find('wmin')) uwmini= of()
         if (find('mini')) umini = of()
      end if
c     
c     t.........parameters related to the talk() subroutine for comunication
c     t.........between workstations (tcpip protocol)  
c     t
c     t         proj_id = geti('#prj',proj_id)
c     t         max_proc = geti('#mxp', max_proc)
c     t         proc_id = geti('#pro',proc_id)
c     t         me      = getchar('here',me,lme)
c     t         upper   = getchar('uper',upper,lup)
c     t         lower   = getchar('lowr',lower,low)
c
      nstep   = geti('#ste',nstep)
      neqstep = geti('#equ',neqstep)
      ninfo   = geti('info',ninfo)
      ncoor   = geti('#crd',ncoor)
      nvelo   = geti('#vel',nvelo)
      nlist   = geti('#lis',nlist)
      nscalv  = geti('#scl',nscalv)
      ntemp   = geti('#tmp',ntemp)
      nrigi   = geti('#rig',nrigi)
      irand   = geti('rand',irand)
      newv    = geti('newv',newv)
      dt      = (getd('step',(dt)))
      hydro_scale = getd('hscl',hydro_scale)
      fmax     = getd ('fmax',fmax)
      tempi(1) = getd('tmpi',tempi(1))
      tempf(1) = getd('tmpf',tempf(1))
      if (find('hvdw')) hvdw0 = .false.
c   
      nstepopt   = geti('#sts',nstepopt)
      npri    = geti('#pri',npri)
      nwcrd   = geti('#wcr',nwcrd)
      smlp    = geti('#mlp',smlp)
      ratio   = geti('rati',ratio)
c     out         istart  = geti('strt',istart)
c     out         istru   = geti('stru',istru)
      tolg    = getd('tolg',tolg)
      estred  = getd('dfpr',estred)
      gamma   = getd('gama',gamma)
      stoscale  = getd('stos',stoscale)
      stobeta   = getd('stob',stobeta)
      dtsto     = getd('dtst',dtsto)
      if (find('ovlp')) lap   = .true.
c     pvm         if (find('nofi')) first = .false.
c     pvm         if (find('nola')) last  = .false.
      if (find('bcd1')) bc1   = .true.
      if (find('bcd2')) bc1   = .false.
      if (find('damp')) damp   = .true.
      if (find('sdpy')) sdpyes   = .true.
c     
c     
c........Energy parameters
c     
      nliststo    = geti('lsts',nliststo)
      cutvdw2  = getd('rvmx',cutvdw2)
      cutvbig2 = getd('rvbg',cutvbig2)
      cutele2  = getd('relx',cutele2)
      cutebig2 = getd('rebg',cutebig2)
      cutmono2 = getd('cutm',cutmono2)
      hydro_scale = getd('hscl',hydro_scale)
      rmax     = getd('rmax',rmax)
      eps      = getd('epsi',eps)
c     
      if (find('cdie')) ctrue = .true.
      if (find('rdie')) ctrue = .false.
c     
      if (find('cpth')) crdstyl = 'PATH'
      if (find('cchr')) crdstyl = 'CHAR'
      if (find('cini')) crdstyl = 'INIT'
      if (find('cint')) crdstyl = 'INTR'
      if (find('cext')) crdstyl = 'EXTR'
      if (find('ctta')) crdstyl = 'TETA'

      if (find('nobo')) ebyes  = .false.
      if (find('noan')) ethyes = .false.
      if (find('noto')) etoyes = .false.
      if (find('noim')) eimyes = .false.
      if (find('novd')) evdyes = .false.
      if (find('noel')) eelyes = .false.

      if (find('dyna')) dyna = .true.
      if (find('noin')) noinert  = .true.
      if (find('prep')) prep  = .true.
      if (find('mshk')) shakm  = .true.
      if (find('shkb')) shakb = .true.
      if (find('symm')) then
         esymyes = .true.
         a = getd('xtra',0.0d0)
         b = getd('ytra',0.0d0)
         c = getd('ztra',0.0d0)
c     write(stdo,*)'esymyes ',esymyes
      end if
c....................................................................     
czva.....initialize slow/fast separation for action using "nfrz"
c
c     NOTE: currently fast modes (e.g. those associated with waters
c     and/or COUNTERIONS) must come AFTER slow ones i.e.
c     'fast' particles (waters) must be located after slow part in con
c     This is not checked below!!!
c     NOTE: as one can check in getvec, it is assumed that there is no
c     covalent realtionship (i.e. bonds, angles ...) between slow and
c     fast particles 
c     if such a relation does occur (sidechains in fast!!!), loops in
c     corresponding subroutines calculating second derivative contrib.
c     (i.e. getvec, bond2_pcs, ...) must be modified, such that for
c     instance slow_i - fast_i bond gives rise to diagonal cotntribution
c     d/dslow_i * d/dslow_i
c     NOTE: as currently only waters and ions were used in fast (no torsions
c     in the fast part for instance) only loops for bonds and angles
c     have been modified in getvec (in the covalent part)
c----------------------------------------------------
c
czva  The "frozen" particles are not used for dynamics. The remaining
czva  particles are considered as fast !!! (pick_fast is droped) 
c
      if (find('nfrz')) then
         call frz_eval_sto()
         write(stdo,'(/20x,a,i6/)') 'Number of slow prts=',npts 
      end if

c     
c........Get update flag
      if (find('updt')) allupdate = .true.
c     
c........Annealing Parameters
c     
      if (find('lann')) lanneal = .true.
      dtopt     = getd('dtop',dtopt)
      tempanl   = getd('tanl',tempanl)
      irandsto  = geti('irds',irandsto)
      ianlmax   = geti('anlm',ianlmax)
      anlfac    = getd('anlf',anlfac)     
c     
      if (find('acti')) go to 2
      go to 1
 2    continue
c        
c     
      write(stdo,'(/5x,a)') 'Elapsed time after Input Action !!!'
      call printcpu(stdo)


      if (nbeta.gt.0) ehyes = .true.
      bc2=.not.bc1
c
      dt    = dt/tconv
      twodt = dt*0.5
      dt2   = dt*dt*0.5d0
      do 700 i=1,npt
         factor1(i) = dt2*invms(i)
         factor2(i) = twodt*invms(i)
 700  continue
c   
c     
c.....re-normalize dtsto,gamma and the masses
c   
      dtsto    = dtsto/tconv
      twodtsto = 2.d0*dtsto
      gamma = gamma/(2.d0*dtsto)
      stoscale=1.0d0/(stoscale**2)
      if (.not.noinert) then
         do 21 i=1,npt
            mdt2(i)      = ptms(i)/(dtsto*dtsto)
            mdt2(i+npt)  = mdt2(i)
            mdt2(i+npt2) = mdt2(i)
c            invms(i) = 1.d0
c     this is for mshake and shake
 21      continue
      else
         do 211 i=1,npt
            mdt2(i)      = 0.d0
            mdt2(i+npt)  = 0.d0
            mdt2(i+npt2) = 0.d0
c            invms(i) = 1.d0
 211     continue
      end if

c     this is for mshake and shake
      
      if (shakm) call mshakinit(debug)

c NOTE: In sto shake constrains the values of bond AND angles
c by keeping the distance of 1-2 and 1-3 fixed.
c
      if (shakb) call shakinit(shakl,shakb,shakm)

c.....initialize the cut-off distances:
c     
c.....if rmax is given (in case of working with a single cutoff distance
c.....for all long range interactions)---maintained here for old input 
c     
      if (rmax.gt.0) then
         cutvdw2 = rmax
         cutele2 = rmax
      end if
c     
      if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
      if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0
c     
      cutvdw2  = cutvdw2*cutvdw2
      cutvbig2 = cutvbig2*cutvbig2
      cutele2  = cutele2*cutele2
      cutebig2 = cutebig2*cutebig2
      if (cutmono2.lt.0) then
         cutmono2 = cutebig2*1.44
      else
         cutmono2 = cutmono2*cutmono2
      end if

c  
c.....check that required files were opened. In case of error stop.
c     
      if (my_pe.eq.0 .and. .not.fopen(ucon)) then
c     c      if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
      elseif(my_pe.eq.0.and.(.not.(fopen(urcrd) .or. urcrd.eq.5))) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)
      elseif (.not.fopen(uwcrd)) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)
      end if
c     
      if ((my_pe.eq.0).and.(urcrd .eq. 99) ) then
         write(stdo,*)' Sorry, initial coordinate file not provided'
         write(stdo,*)' missing  - urcrd number '
         write(stdo,*)' Quitting...., better luck next time'
         level = 1
         call alert(name,namel,'Problem in file',15,level)
c     else if (mod(grid,2) .eq. 0) then
c     write(stdo,*)' Sorry, number of grid points must be odd '
c     write(stdo,*)' Input is grid = ',grid
c     write(stdo,*)' Quitting..., better luck next time'
c     level = 1
c     call alert(name,namel,'Problem in grid',15,level)
      end if
      write(stdo,100)nstepopt,npri,grid,stdo,uwcrd,urcrd,crdstyl
     1     ,nwcrd,proc_id
c     
c.....check if istart and istru (defined in init_paral) were redefined.
c.....this can cause some mess if running in parallel+interpolation,
c.....or parallel+extrapolation (e.g. in multigrid mode)
      write(stdo,102) grid,pseg,istart,istru
c     
c     
c     
c.....get initial coordinates (call rchain)
c.....depending on "crdstyl", the following options exist
c.....(i)   DYNAmics files (for movies using QUANTA)
c.....(ii)  PATH files (binary, double precision code to maintain 
c.....accurate coordinates
c.....(iii) INIT. Reading formatted coordinates file for reactants and
c.....products and generating the rest of the path by
c.....linear interpolation
c.....(iv)  INTRpolate. Given a low resolution path in PATH format, add
c.....structures in between to refine the path.
c.....(v)   TETA (initialize a chain using a theta-function 
c.....interpolation,first half of the grid at reactants, second 
c.....half at products. 
c.....(vi)  EXTRaction. The amount of structures is reduced by 
c.....taking a less dense grid of structures with some skip 
c.....factor.
c.....(crdstyl is defined in init_sto.f as INIT)
c     
c     
      write(stdo,'(/5x,a)') 'Elapsed time Before Reading Chain !!!'
      call printcpu(stdo)

      call rchain(urcrd,istart,istru,ratio)

      write(stdo,'(/5x,a)') 'Elapsed time After Reading Chain !!!'
      call printcpu(stdo)
c
      write(stdo,101)tolg,gamma*2*dtsto,dtsto*tconv,
     >                           nliststo,istart,istru
      write(stdo,*) ' debug ? ',debug
c     
c.....output the energy terms used
      call init_wre(stdo)
c     
      tolg = tolg*tolg
c     
c     
c.....output number of processors used
      if (first.and.last) then
c........use only one processor (irrespective of proc)
c........yet proc needs to be defined correctly since some matrices are
c........defined using proc (!)
c     
         write(stdo,*)' *** ONLY ONE PROCESS IS USED '
         ndegf = (grid-4)*npt3
      else 
         write(stdo,*)' *** ',proc_max,' PROCESSES ARE USED '
         ndegf = pseg*npt3
         if (mod(grid-4,proc_max).ne.0) then
            level = 1
            write(stdo,'(10x,a)') " (grid-4)/proc MUST be an integer "
            call alert(name,namel,'Illegal (grid-4)/proc',21,level)
         end if
      end if
c
      write(stdo,'(10x,a,i6)') "# of degrees of freedom = ",ndegf
      write(stdo,'(10x,a,e12.5)') " ST0_scale=",stoscale
      write(stdo,'(10x,a,e12.5)') " ST0_beta=",stobeta
c
      write(stdo,'(10x,a,e12.5)') "cutvdw2  =",cutvdw2
      write(stdo,'(10x,a,e12.5)') "cutvbig2 =",cutvbig2
      write(stdo,'(10x,a,e12.5)') "cutele2  =",cutele2
      write(stdo,'(10x,a,e12.5)') "cutebig2 =",cutebig2
      write(stdo,'(10x,a,e12.5)') "cutmono2 =",cutmono2
c
ccc      do 1007 i=1,totdmon
ccc         write(stdo,*) " Input STO debug : Frozen flag =",i,frzM(i)
ccc 1007 continue
c
c-----------------------------------------------------------------     
c      if (prep) then
c     
c........In case of running large systems with INIT is necessary in
c........some cases to relax the structures in the path (a little) to 
c........avoid non-realistic deformed structures, which have a very 
c........large contribution to initial action.
c     
c      write(stdo,106)
c 106  format(1x,/,1x,'**** calling the preprocessing',
c     &     ' ( to partially relax the', 'initial structures)',/)
c     
c      call preproc(umini,uwmini)
c  
c      endif 
c----------------------------------------------------------------
c 
      write(stdo,'(/5x,a)') 'elapsed time after the preprocessing'
c
      call printcpu(stdo)
c     
      if (my_pe.eq.0.and.fopen(ucon)) close(ucon)
      if (fopen(urcrd)) close(urcrd)
c
c    
 100  format(/,1x,'PARAMETERS FOR STOCHASTIC CHAIN MINIMIZATION: ',//,
     1     1x,' number of minimization steps: ',i7,/,
     2     1x,' print each ',i7,' steps',/,
     3     1x,' number of grid points of the chain : ',i5,/,
     4     1x,' data is on unit: ',i5,' coordinates are on unit ',i5,/,
     5     1x,' initial coordinates are read from unit : ',i5,/,
     6     1x,' Coordinate are read in ',a4,' style ',/,
     7     1x,' coordinates are written at step interval : ',i5,/,
     8     1x,' Processor identification number : ',i5)
 101  format(1x,' gradient max. error ',e10.3,/,
     2     1x,' friction coefficient: ',e10.3,/,
     3     1x,' time step  - ps -(along the path) ',e10.3,/,
     4     1x,' update non-bonded list each ',i6,' steps',/,
     5     1x,' the starting structure for current processor ',i5,/,
     6     1x,' number of structures ',i5)
 102  format(1x,' before rchain grid=',i6,' pseg=',i6,
     2     1x,' istart = ',i6,' istru = ',i6)
c     
        return
        end
c

