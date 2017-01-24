      program Mfep

      implicit none

c     
c     calculate path using simulation of a chain.
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/METAL.BLOCK'
      include 'COMMON/RGYRCST.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'


c     
c     local space
c     Be careful with space allocation which is
c     largely based on "maxpt" - maximum number of particles
c     which is in LENGTH.BLOCK. "maxpt" is usually set to a large
c     number (15,000). For a long chain that will be a LOT of space.
c     You may wish to modify
c     LENGTH.BLOCK for your own purposes. Decision on proper
c     lengths of vectors (for your molecule) should not be
c     difficult, since the connectivity file must be avaliable.
c     
c     **** ALLOCATE SPACE
c     define double precision vectors of the length maxpt*3
c     temp  = desired temperature
c     dt    = time step
c     ***special parameters for free energy calculations along RC ***
c     umlst - a unit number of a file with reaction coordinate
c     urcrd - (read) file containing trajectory starting points
c     nrcrd - number of structures in urcrd
c     uwcrd - (write) file to which to write trajectories
c     uwpep - (write) file to which to write peptide-only trajectories
c     nwcrd - number of timesteps between FP traj writes
c     nsvel - number of steps between velocity rescaling
c     newv - number of steps between velocity resampling
c     mxstps - max steps; default is -1, meaning no max.
c     orth  - boolean indicating orthogonal constraint
c     *** NOTE THAT IN urcrd NPATH+1 POINTS ARE REQUIRED (INCLUDING
c     PRODUCTS) ***
c     igrid - number of path segments (active milestones)
c     num_mlsts - igrid + 1 (number of milestones)
c     debug = logical variable: true= a lot of information printed out
c     x     = coordinate vector
c     y     = coordinate vector
c     z     = coordinate vector
c     vx   = velocity vector
c     vy   = velocity vector
c     vz   = velocity vector
c     dx    = forces
c     dy    = forces
c     dz    = forces
c     dxold = old forces (previous step)
c     dyold = old forces (previous step)
c     dzold = old forces (previous step)
c     dmass = masses
c     divms = 1 over the mass
c     fact11 = a constant vector used in the Verlet integration
c     fact22 =                    "
c     the reaction coordinate at position numpth.
c     Estimated here by finite difference.
c     pstep    - a step vector translation from point numpth to
c     numpth+1
c     
c     define vectors for constraints
c     grdcm[x-z] = derivative of the center of mass with respect to [x-z]
c     grdl[x-z]  = derivatives of infitesimal rotations with respect to [x-z]
c     pointr ipick - selection of particle
c     ntest - period to test constraints
c     irand- integer seed for random number generator of velocities
c     stdo - where some information on the path is printed out
c     frep - update plane every frep steps 

      integer n
      integer igrid,num_mlsts
      double precision dst2,trceps,e0
      double precision rotat(3,3),ident(3,3)
      double precision dpold(3*maxpt)
      double precision dmass(maxpt),divms(maxpt)
      double precision fact11(maxpt),fact22(maxpt)
      double precision grdcmx(3*maxpt),grdcmy(3*maxpt)
      double precision grdcmz(3*maxpt)
      double precision grdlx(3*maxpt),grdly(3*maxpt)
      double precision grdlz(3*maxpt)
      double precision pstep(3,maxpt)
      integer nomlst
      parameter(nomlst=0)

      integer ierr
      integer status(MPI_STATUS_SIZE)


      integer pointr(maxpt),ipick(maxpt)

      double precision temp,dt,dt2,tfac, vscalar

c     pre-equilibration temperature
      double precision ptemp

      character*80 name
      integer namel

      integer ucon,umlst,ufpt,ufpp,utrc,udot,udot1,udst,nene
      integer nwcrd,uwcrd,uwpep,urcrd,nrcrd,uene,umom,nwpep
      integer upfr,upvl

c     this, previous, and next milestone numbers
      integer tmlst, pmlst, nmlst

      logical tmplog1,tmplog2,trc
      logical orth
      logical tfp
      integer npri,ntest,mxstps
      integer level,irand,npick,nlist
      integer nsvel, newv, neqstep,nmom

      integer geti,of
      double precision getd
      logical find,fopen

c     for trace
      double precision norm1, norm2, dp

c     boolean deciding between ucrb and phi/psi approaches
      logical ucrb

      double precision scalar(7),sigma(7)
      integer i,j,k,npick3,npt3
      logical massw

      logical shakl,shakb,shakm,freeze,toiyes,nori
      double precision tolshk
      integer ndegf,ii

      data tfac/4.2532764d-2/
      data uene,utrc,upvl/3*99/
      data ucon,umlst,upfr,ufpt,ufpp,urcrd,uwpep,uwcrd/8*99/

c     plp, pln (prp, prn) : left (right) fp plane
c     point and normal
      double precision plp(3,maxpt), pln(3,maxpt)
      double precision pmn(3,maxpt), pmp(3,maxpt)
      double precision pmnn(3,maxpt), pmpn(3,maxpt)
      double precision pmn0(3,maxpt), pmp0(3,maxpt),pmpt(3,maxpt)

      double precision prp(3,maxpt), prn(3,maxpt)
      double precision unitgrad(3,maxpt)
      double precision gamma,smooth

c     lnn (rnn): left (right (mid)) normal norm
      double precision lnn, rnn, mnn

      integer proc,ik,nwav,frep 
      logical paral,repul,fixend  

c     average coordinate of the current plane one before and one after  
      double precision avg_pep(3,maxpt),curr(3,maxpt)
     &,prev(3,maxpt),next(3,maxpt) 

     
c     
c     General initialization
c     
      stdi   = 5
      stdo   = 6
      stderr = 0
      totmon = 0
      npt    = 0
      nb     = 0
      nangl  = 0
      ntors  = 0
      nimp   = 0
      lestyp = 0
      nbulk  = 0

      debug  = .false.
      nocut  = .false.
      name   = 'mfep'
      namel  = 5

      shakb  = .false.
      shakl  = .false.
      tolshk = 1.d-6

      shakm  = .false.
      tolcons= 1.d-5

      lcent  = .false.
      freeze = .false. 
      toiyes = .false.
      nori   = .false.
      trc    = .false.
c     qssbp  = .false.
      ncnst  = 0

      nwcrd  = 0
      nwpep  = 0

      umom   = -1
      trceps = 0.5
      abarpsi = -2.0

      paral = .false.
      fixend = .false. 
      repul = .false. 

      proc =1
      gamma=0.00005d0
      smooth=0.3d0
      nwav = 400000
      frep  =  1
    

c ....initialization of the parallel machine
        
        call init_paral2(procID,proc)
      
        call open_inout(stdi,stdo,procID)

      
c     
c     open scratch file for line manipulation
c     
      jnkf = 25
      open (unit=jnkf,status='scratch')
c     
c     free energy default parameters
c     
      igrid   = lgrid
      nsvel   = 100
      newv    = 1000
      neqstep = 0
      nrcrd   = 1
      nene    = 100
      npri    = 1
      ntest   = 0
      mxstps  = -1
      nmom = 10
      
      dt      = 0.0005d0
      irand   = 23518284
      temp    = 300.d0
      ptemp   = 100.d0 
      tmlst   = nomlst
      pmlst   = nomlst
      nmlst   = nomlst
      orth    = .false.
      tfp     = .false.
      ucrb    = .false.
      torscstr = .false.
      alabar  = .false.
      rgcst   = .false.
c     
c     energy default parameters
c     
      call init_ef()

      nlist  = 1

c     
c     start interpret the line from stdi
c     

 1    continue


      call rline(name,namel,stdi)
      


      if (find('debu')) debug = .true.
      if (find('file')) then
         if (find('conn')) then
            ucon = of()
            call rconn(ucon)
         endif
         if (find('mlst')) umlst = of()
         if (find('upfr')) upfr = of()
         if (find('upvl')) upvl = of()
         if (find('rcrd')) urcrd = of()
         if (find('wcrd')) uwcrd = of()
         if (find('wmom')) umom = of()
         if (find('wdt1')) udot1 = of()
         if (find('wdot')) udot = of()
         if (find('wdst')) udst = of()
         if (find('wpep')) uwpep = of()
         if (find('wene')) uene  = of()

      end if

      if ( find('ucrb') ) ucrb = .true.
      if ( find('tcst') ) torscstr = .true.

      if ( find('abar') ) alabar = .true.
      abarpsi = getd('abps', abarpsi)


      trceps  = getd('trce',trceps)
      temp    = getd('temp',temp)
      ptemp   = getd('ptmp',ptemp)
      igrid   = geti('grid',igrid)
      nene    = geti('#nen',nene)
      nmom    = geti('#mom',nmom)
      nrcrd   = geti('#rcr',nrcrd)
      nwcrd   = geti('#wcr',nwcrd)
      nwpep   = geti('#pep',nwpep)
      newv    = geti('newv',newv)
      nsvel   = geti('#scl',nsvel)
      neqstep = geti('#equ',neqstep)
c      npri    = geti('#pri',npri)
      ntest   = geti('#tes',ntest)
      dt      = (getd('step',(dt)))
      irand   = geti('rand',irand)
      mxstps  = geti('#mxs', mxstps)
      gamma  = getd('gamm',gamma) 
      smooth  = getd('smth',smooth) 
      nwav  = geti('nwav',nwav) 
      frep  = geti('frep',frep) 
     
      tmlst   = procID+1
      pmlst   = tmlst-1
      nmlst   = tmlst+1
     

   

c     Energy parameters
      nlist   = geti('list',nlist)
      cutvdw2  = (getd('rvmx',(cutvdw2)))
      cutvbig2 = (getd('rvbg',cutvbig2))
      cutele2  = (getd('relx',(cutele2)))
      cutebig2 = (getd('rebg',cutebig2))
      cutmono2 = getd('cutm',cutmono2)
      rmax     = getd('rmax',rmax)
      eps     = (getd('epsi',(eps)))


      if (find('fixd')) fixend = .true.
      if (find('repl')) repul = .true.

      if (find('nocut')) nocut = .true.
      if (find('cdie')) ctrue = .true.
      if (find('rdie')) ctrue = .false.
C     deb------------
      if (find('nori')) nori = .true.

      if (find('nobo')) ebyes  = .false.
      if (find('noan')) ethyes = .false.
      if (find('noto')) etoyes = .false.
      if (find('noim')) eimyes = .false.
      if (find('novd')) evdyes = .false.
      if (find('noel')) eelyes = .false.

      if (find('cnst')) ecnyes = .true.
      if (find('TORS')) then
         write(*,*)
         write(*,*)'* Note that the keyword TORS must come after amid'
         write(*,*)'* If amid is used'
         write(*,*)
         ncnst = ncnst + 1
         if (ncnst.gt.maxcnst) then
            level = 1
            call alert(name,namel,'Maxcnst exceeded',16,level)
         endif
         icnst1(ncnst) = geti('atm1',0)
         icnst2(ncnst) = geti('atm2',0)
         icnst3(ncnst) = geti('atm3',0)
         icnst4(ncnst) = geti('atm4',0)
         kcns(ncnst)   = (getd('kcns',0.d0))
         cnseq(ncnst)  = (getd('cneq',-999.d0))/pi180
         if (debug) then
	    write(stdo,*)' kcns cnseq ',kcns(ncnst),cnseq(ncnst)*pi180
         end if
      end if
      
      if (find('mshk')) then
         shakm = .true.
         tolcons=getd('mtol',1.d-5)
         write(stdo,66)tolcons
 66      format (1x,'matrix shake on TIP3/SPCE, tolcons=',d7.1)
         call mshakinit(debug)
         if (nshakm.le.0) then
            write (stdo,*) 'nothing to shake!'
            shakm=.false.
         endif
      endif


c     added by TF, 7-Nov-2004
      if (find('shkb')) shakb = .true.


      if (find('symm')) then
         esymyes = .true.
         a = getd('xtra',0.0d0)
         b = getd('ytra',0.0d0)
         c = getd('ztra',0.0d0)
      end if


      if (find('hvdw')) hvdw0 = .false.


c     switch for Ewald summation of long range interactions
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
c     write (stdo,*) 'PME account for long range inter. will be performed'

c     if one wants to use PME for vacuum calculations
c     one needs a virtual box which is sufficiently large to effectively
c     separate image boxes
c     in such a case stop should be commented and proper sizes
c     of the virtual box  xtra,ytra,ztra should be given

         if (.not.esymyes) then
c     write (stdo,*) 'No true periodic bound. cond. Symm is missing'
            stop
c     a = getd('xtra',50.0d0)
c     b = getd('ytra',50.0d0)
c     c = getd('ztra',50.0d0)
         end if
      end if
c     jmjmend
      
c     
c     virtual particles: massless, point charges displaced from the motion
c     centers (real particles) e.g. TIP4P, 3-site CO model
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
c     end of virtual partcles - one may actually think of moving this to con
c     

c     carlos, add amide plane constraints
      if (find('amid')) call amid()

      if (find('metl')) then
         metalyes   = .true.
         a0_metal   = getd('amtl',50.d0)
         alfa_metal = getd('alfa',1.d0)
         b_wall     = getd('bwal',b)
         if (b_wall.gt.b) then
            call alert('dyna',4,'b_wall must be < b',18,1)
         end if
         v_elec     = getd('v_el',0.d0)
c     
c     compute pre-exponential  factor A0_metal such that
c     A0_metal*exp(-alfa_metal*b/2)=a0_metal
c     
c     a0_metal = a0_metal*dexp(-0.5d0*alfa_metal*b_wall)
         write(stdo,*)
         write(stdo,104)
 104     format(1x,' Metal walls perpendicular',
     1        ' to the Y axis will be set!')
         write(stdo,105)a0_metal
 105     format(1x,' Metal repulsive Wall, A0 = ',
     1        E12.6)
         write(stdo,106)b_wall/2,v_elec
 106     format(1x,' Walls are located at +/- ',f9.3,
     1        1x,' The Electrode potential is ',f10.4)
      end if
      
c     
      if (find('cent')) then
         lcent = .true.
         kcenter=getd('kcnt',10.d0)
         xeq = getd('xeqm',0.d0)
         yeq = getd('yeqm',0.d0)
         zeq = getd('zeqm',0.d0)
         call pick(ipick,i)
         icenter=0
         do 35 i=1,npt
            if (ipick(i).ne.0) then
               icenter = icenter+1
               center(icenter)=i
            endif
 35      continue
      endif
c     
      if (find('nfrz')) then
         freeze = .true.
         call pick(ipick,i)
         inofrz = 0
         do 44 i=1,npt
            if (ipick(i).ne.0) then
               inofrz = inofrz + 1
               nofreez(inofrz) = i
               zerofrz(i) = 1
            else
               zerofrz(i) = 0
            end if
 44      continue
      end if
      if (find('toiy')) toiyes = .true.
      if (find('selc')) then
         call pick(ipick,npick)
         npick=0
         do 3 i=1,npt
            npick=npick+ipick(i)
 3       continue
         j=0
         do 4 i=1,npt
            if (ipick(i).gt.0) then
               j=j+1
               pointr(j)=i
            end if
 4       continue
      endif
      

C     deb-------------------------------------------------------
C     deb+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Spherical Solvent Boundary Potential (SSBP)
C     
C     Authors:  Dmitrii Beglovd and Benoit Roux  (1994)
C     Department of Chemistry, University of Montreal
C     

C     Default is include all energy terms for the 
C     solvent boundary potential
C----------------------------------------------------------
c     if (find('ssbp')) call init_ssbp(ipick)
C     deb-------------------------------------------------------

      if (find('orth')) orth = .true.
      if (find('acti')) go to 2
      go to 1
 2    continue

c     END INPUT LOOP
c     -------------------------------------------------------------


      num_mlsts = igrid + 1

c     reaction coordinate test
c      if ( trc ) then
c         write(*,*) 'This is a TRC run, with no MD'
c         if ( .not. fopen(utrc) ) then
c            level = 1
c            call alert(name,namel,'wtrc not opened',16,level)
c         end if
         
c        The next 2 commented lines of code were the old trace,
c        which was useful for alanine dipeptide
c        write(*,*) 'trceps =', trceps
c        call wmlst(utrc,pmp,pmn,npick,trceps)

c         call getmlst(umlst,pmp,pmn,tmlst,num_mlsts,npick,pointr)
c          call getmlst(umlst,prp,prn,nmlst,num_mlsts,npick,pointr)
         
c         norm1 = 0.d0
c         norm2 = 0.d0
c         dp = 0.d0
c         do 34 j=1,npick
c            do 344 k=1,3
c               norm1 = norm1 + pmn(k,j)*pmn(k,j)
c               norm2 = norm2 + prn(k,j)*prn(k,j)
c                  dp =    dp + pmn(k,j)*prn(k,j)
c 344        continue
c 34      continue
c         dp = dp / sqrt ( norm1 * norm2 )
c         write(utrc,*) dp
c         write(*,*) 'FP STOP (TRC)'
c         stop
c      end if
      

c     if necessary, prepare the pointr for rgyrconst
      if ( rgk .ne. 0.d0 ) then
         if ( torscstr .or. ucrb ) then
            level = 1
            call alert ( name, namel, 'rgk w/ torscstr/ucrb',40,level)
         end if
         if ( rg20 .eq. 0.d0 ) then
            level = 1
            call alert ( name, namel, 'rgk w/o rg20', 40, level )
         end if
         rgcst = .true.
         rgcnpick = npick
         rgcnpt = npt
         do 835 i = 1,rgcnpt
            rgcnpointr(i) = pointr(i)
 835     continue
      end if
      write(*,*) 'rgcst: ', rgcst
      write(*,*) 'ucrb, torscstr: ', ucrb, torscstr


c     adjust pre-equilibration temperature
      if ( ptemp .lt. 0 ) ptemp = temp

      if ( orth .and. tfp ) then
         level = 1
         call alert ( name, namel, 'orth and tfp', 40, level )
      end if


      if ((tmlst .lt. 0) .or. (nmlst .lt. 0) .or. (pmlst .lt. 0)) then
         level = 1
         call alert(name,namel,'mlst number problem',20,level)
      end if


c     check that required files were opened
      if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)

      else if (.not.fopen(umlst)) then
         level = 1
         call alert(name,namel,'umlst not opened',16,level)


      else if (.not.fopen(urcrd)) then
         level = 1
         call alert(name,namel,'urcrd not opened',16,level)


      else if ( .not. fopen(uwcrd) ) then
         level = 1
         call alert(name,namel,'uwcrd not opened',16,level)

      else if ( ( nwpep .gt. 0 ) .and. .not. fopen(uwpep) ) then
         level = 1
         call alert(name,namel,'uwpep not opened',16,level)

      else if ( .not. fopen(uene) ) then
         level = 1
         call alert(name,namel,'uene not opened',16,level)

      end if



C     deb-------------------------------------------------------
c     initialze no-freeze vector
c     set the pointer to the selected particles. pointr(i) is 
c     the position of selected atom number i in the normal all 
c     atom array

      if (.not.freeze) then
         inofrz = npt
         do 21 i=1,inofrz
            nofreez(i) = i
            zerofrz(i) = 1
 21      continue
      endif


c     now that nofreez is properly defined, handle shake
      if (shakb) then
         itershak = 100
         epshak = 1.d-4
         epshakv = 1.d-4
         call shakinit(shakl, shakb, shakm, tolshk)
         if (nshak .le. 0) then
            shakb = .false.
            write (stdo,*) 'WARNING: nothing to shake!'
         end if
      end if


C     deb-------------------------------------------------------
c     rmax is maintained here for old input (with a single
c     cutoff to work)
c     
      if (rmax.gt.0) then
         cutvdw2 = rmax
         cutele2 = rmax
      end if
      
c     
c     internally the 1-4 scaling factors are defined as the inverse
      if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2*1.2d0
      if (cutebig2.lt.0.d0) cutebig2 = cutele2*1.2d0


      cutvdw2  = cutvdw2*cutvdw2
      cutvbig2 = cutvbig2*cutvbig2
      cutele2  = cutele2*cutele2
      cutebig2 = cutebig2*cutebig2

      if (cutmono2.lt.0) then
         cutmono2 = cutebig2*1.44
      else
         cutmono2 = cutmono2*cutmono2
      end if
      
      if (npick.eq.0) then
         if (debug) write(stdo,*) ' ipick = ', (ipick(j),j=1,npt)
         level = 1
         call alert(name,namel,'No selection of particles',25,level)
      end if


      write(stdo,100) temp,igrid,nmom,stdo,umlst,urcrd,uwcrd,uwpep,uene
      write(stdo,101) dt,ntest,irand,nwcrd,ptemp,neqstep,nlist,nsvel,
     1gamma,smooth,frep,nwav
   
      write(stdo,*) 'debug ? ', debug
 100  format(/,1x,' PARAMETERS FOR FREE ENERGY SIMULATION:',//,
     1     1x,' temperature: ',f5.1,/,
     2     1x,' number of milestones in reaction coordinate: ',i5,/,
     3     1x,' print distances each ',i7,' steps',/,
     4     1x,' data is on unit: ',i5,/,
     5     1x,' mlst coordinates are on unit ',i5,/,
     6     1x,' ensembles are read from unit ',i5,/,
     7     1x,' trajectory is written to unit ',i5,/,
     8     1x,' peptide traj is written to unit ',i5,/
     9     1x,' energy is written to unit ',i5,/)
 101  format(1x,' step size: ',e10.3,/,
     1     1x,' period for checking constraints ',i5,/,
     3     1x,' initial seed for random vel. selection ',i10,/,
     4     1x,' number of timesteps between traj write',i10,/,
     5     1x,' pre-equilibration temperature: ',f5.1,/,
     6     1x,' pre-equilibration steps: ',i5,/,
     7     1x,' steps between velocity re-scaling: ',i5,/,
     7     1x,' re-build non-bonded list every: ',i5,/,
     7     1x,' time step of plane updates gamma: ',e10.3,/,
     8     1x,' degree of smoothing the string, smooth : ',e10.3,/,
     8     1x,' update planes every : ',i5,/,
     8     1x,' reinitiate the average every, nwav : ',i10,/)
c            write(stdo,*) ' first and the last milestones are fixed:',fixend
c            write(stdo,*) ' repulsion btw i i+2 :',repul,smooth
    
 


      dt = dt / tfac

c     jmjm
      if ( ewaldyes ) call ewald_init()

c     end of jmjm
      if ( vp_flag ) call vp_init()
c     end of vp

      npick3 = 3 * npick
      npt3   = 3 * npt
      if ( debug ) write(stdo,*) ' npick npt  ',npick,npt
c     
c     do mass weighting 
c     
      massw=.true.
c     
c     calculate the constants fact11 and fact22 dmass and 1/m (divms)
c     
      dt2 = dt * dt / 2.d0

      do 6 i=1,npt
         dmass(i) = ptms(i)
         divms(i) = 1.d0/dmass(i)
         invms(i) = divms(i)
         fact11(i) = dt2 * divms(i)
         fact22(i) = dt * divms(i) / 2.d0
 6    continue

c     IF (QSSBP) call set_ssbp()

      if ((tmlst .gt. num_mlsts) .or. (pmlst .gt. num_mlsts) .or.
     $     (nmlst .gt. num_mlsts)) then
             if(procID.ne. (proc-1)) then 
         level = 1
         call alert(name, namel, 'mlst > num_mlsts', 15, level)
      end if
      endif 

     


c     get initial path structure,
c     path derivatives and the path step, vx vy vz & dxold dyold dzold
c     are used as temporary vectors.
C     deb-------------------------------------------


       read(umlst) e0,((pmp(k,i),i=1,npt),k=1,3) 

       call deselect(pmp,pointr,npick)

      
      if(procID . eq. (proc-1) ) then

      do i=1,npick
      do k=1,3
       next(k,i)=pmp(k,i)
      enddo
      enddo
         Call MPI_Recv(prev, 3*npick, MPI_DOUBLE_PRECISION,
     >                        proc-2,1,MPI_COMM_WORLD, status,ierr )
      endif
       if(procID.eq.(proc-2)) then
          Call MPI_Send(pmp, 3*npick, MPI_DOUBLE_PRECISION,
     >                         (proc-1),1,MPI_COMM_WORLD, ierr )
       endif


cc     first milestone does not have prev 
       if(procID . eq. 0 ) then
      do i=1,npick
      do k=1,3
      prev(k,i)=pmp(k,i)
      enddo
      enddo
         Call MPI_Recv(next, 3*npick, MPI_DOUBLE_PRECISION,
     >                        1,23,MPI_COMM_WORLD, status,ierr )
      endif
       if(procID.eq.1) then
          Call MPI_Send(pmp, 3*npick, MPI_DOUBLE_PRECISION,
     >                         0,23,MPI_COMM_WORLD, ierr )
       endif

       do i=1,proc-2
       if(procID.eq.i) then
         Call MPI_Recv(prev, 3*npick, MPI_DOUBLE_PRECISION,
     >                        i-1,13,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i-1)) then
          Call MPI_Send(pmp, 3*npick, MPI_DOUBLE_PRECISION,
     >                         i,13,MPI_COMM_WORLD, ierr )
       endif

        if(procID.eq.i) then
         Call MPI_Recv(next, 3*npick, MPI_DOUBLE_PRECISION,
     >                        i+1,24,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i+1)) then
          Call MPI_Send(pmp, 3*npick, MPI_DOUBLE_PRECISION,
     >                       i,24,MPI_COMM_WORLD, ierr )
       endif
       enddo


      call getmlst(prev,next,pmp,pmn,npick)


      call normalize ( pmn, 3*npick, unitgrad, norm1 )


    
c     get center of mass constraints and orthonormalize all
c     constraints. (dpold and velo are dummies)
      call comc(pmp,pmn,divms,npt,pointr,npick,grdcmx,
     1     grdcmy,grdcmz,grdlx,grdly,grdlz,dpold(1),
     2     dpold(1+npt),velo,sigma,ucrb,debug)


c      call dyna_mfep(procID,proc,temp,dt,dt2,pmlst,tmlst,nmlst,
      call dyna_mfep(proc,temp,dt,dt2,pmlst,tmlst,nmlst,
     1     nmom,nrcrd,mxstps,nsvel,pointr,
     2     npick,dmass,divms,grdcmx,grdcmy,grdcmz,
     3     grdlx,grdly,grdlz,fact11,fact22,stdo,urcrd,ufpt,ufpp
     $     ,uwcrd,uwpep,uene,nene,nwcrd,ntest,ptemp,neqstep,newv,irand
     $     ,debug,scalar,sigma,nlist,massw,shakm,shakb,shakl,nwpep,prp
     $     ,prn,plp,pln,pmp,pmn,unitgrad,orth,tfp
     $     ,udot,udot1,udst,umom,ucrb,umlst,upfr,upvl,num_mlsts
     $     ,gamma,smooth,nwav,frep,fixend,repul)
     

      
     
 1000 format(1x,8(f9.6))
 2000 continue
      
      write(*,*) 'MFEP STOP'
      
      stop
      end