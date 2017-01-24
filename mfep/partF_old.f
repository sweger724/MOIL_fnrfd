      program PartF 

      implicit none

c     
c     calculate force projected on the path 
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


c     urcrd - (read) file containing trajectory starting points
c     nrcrd - number of structures in urcrd
c     mxstps - max steps; default is -1, meaning no max.
c     igrid - number of path segments (active milestones)
c     num_mlsts - igrid + 1 (number of milestones)
c     x     = coordinate vector
c     y     = coordinate vector
c     z     = coordinate vector
c     dx    = forces
c     dy    = forces
c     dz    = forces
c     dmass = masses
c     divms = 1 over the mass
c     the reaction coordinate at position numpth.
c     Estimated here by finite difference.
c     pstep    - a step vector translation from point numpth to
c     numpth+1
c     
c     define vectors for constraints
c     pointr ipick - selection of particle
c     stdo - where some information on the path is printed out

      integer n,m,istr
      integer igrid,num_mlsts
      double precision dst2,e0
      double precision dmass(maxpt),divms(maxpt)
      double precision pstep(3,maxpt)
      integer nomlst
      parameter(nomlst=0)


      integer pointr(maxpt),ipick(maxpt)
      double precision temp,dt,dt2,tfac, vscalar

      character*20 name
      integer namel

      integer ucon,umlst,udot,uvfr
      integer nstru,uwcrd,uwpep,urcrd,nrcrd,uene,umom,nwpep
      integer upfr

c     this, previous, and next milestone numbers
      integer tmlst, pmlst, nmlst

      integer npick

      integer geti,of
      double precision getd
      logical find,fopen,freeze

      integer i,j,k,npick3,npt3,nlist,level  

      data ucon,umlst,upfr,urcrd,uvfr/5*99/

c     point and normal
      double precision pmn(3,maxpt), pmp(3,maxpt)
      double precision projf 
c     average coordinate of the current plane one before and one after  
      double precision avg_fp(3,maxpt),prev(3,maxpt),next(3,maxpt) 
      double precision dS, tmpx,tmpy,tmpz 
c     
c     General initialization
      stdi   = 5
      stdo   = 6
      stderr = 0
      totmon = 0
      npt    = 0
      nb     = 0
      nangl  = 0
      ntors  = 0
      nimp   = 0
      nbulk  = 0
      nstru  = 1 

      debug  = .false.
      nocut  = .false.
      name   = 'partF'
      namel  = 5

      freeze = .false. 
      umom   = -1
    

c     open scratch file for line manipulation
c     
      jnkf = 25
      open (unit=jnkf,status='scratch')
c     
c     free energy default parameters
c     
      igrid   = lgrid
      tmlst   = nomlst
      pmlst   = nomlst
      nmlst   = nomlst
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
         if (find('rcrd')) urcrd = of()
         if (find('avfr')) uvfr = of() 
      end if

      nstru   = geti('#rcr',nstru)
      igrid   = geti('grid',igrid)
      tmlst   = geti('cmls',tmlst) 
      pmlst   = geti('pmls',pmlst) 
      nmlst   = geti('nmls',nmlst) 

      write(*,*) igrid,tmlst,nmlst,pmlst,nstru 
c     Energy parameters
      cutvdw2  = (getd('rvmx',(cutvdw2)))
      cutvbig2 = (getd('rvbg',cutvbig2))
      cutele2  = (getd('relx',(cutele2)))
      cutebig2 = (getd('rebg',cutebig2))
      cutmono2 = getd('cutm',cutmono2)
      rmax     = getd('rmax',rmax)
      eps     = (getd('epsi',(eps)))

      if (find('nocut')) nocut = .true.
      if (find('cdie')) ctrue = .true.
      if (find('rdie')) ctrue = .false.
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

         if (.not.esymyes) then
c     write (stdo,*) 'No true periodic bound. cond. Symm is missing'
            stop
         end if
      end if
c     
      if (find('amid')) call amid()

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
      

      if (find('acti')) go to 2
      go to 1
 2    continue

c     END INPUT LOOP
c     -------------------------------------------------------------


      num_mlsts = igrid + 1

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


c     rmax is maintained here for old input (with a single
c     cutoff to work)
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


      write(stdo,100) igrid
      write(stdo,101) nstru 
   
      write(stdo,*) 'debug ? ', debug
 100  format(/,1x,' PARAMETERS FOR FREE ENERGY SIMULATION:',//,
     2     1x,' number of milestones in reaction coordinate: ',i5,/)
 101  format(1x,'number of configurations: ',i7,/)
    
      if ( ewaldyes ) call ewald_init()
      npick3 = 3 * npick
      npt3   = 3 * npt
      if ( debug ) write(stdo,*) ' npick npt  ',npick,npt

c     get initial path structure,
c     path derivatives and the path step, vx vy vz & dxold dyold dzold
c     are used as temporary vectors.
C     deb-------------------------------------------

      do j = 1, pmlst 
       read(umlst) e0,((prev(k,i),i=1,npt),k=1,3) 
      enddo 
       read(umlst) e0,((pmp(k,i),i=1,npt),k=1,3)  ! current 
       read(umlst) e0,((next(k,i),i=1,npt),k=1,3) 


      call deselect(prev,pointr,npick)
      call deselect(pmp,pointr,npick)
      call deselect(next,pointr,npick)
! Get a milestone point and its normal
      call getmlst(prev,next,pmp,pmn,npick)
! Get the |ds| from the path 

       do i = 1,npick
                 tmpx = next(1,i) - pmp(1,i)
                 tmpy = next(2,i) - pmp(2,i)
                 tmpz = next(3,i) - pmp(3,i)
                 dS = dS + tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
        enddo 
       dS=sqrt(dS)       

c  initial averge force array is zero  
       do k = 1,3
          do j = 1,npt 
           
               avg_fp(k,j)=0.d0 

          end do
        end do



      istr=0 
c     read the coordinate 
      do j = 1, nstru
       read(urcrd) e0,((coor(k,i),i=1,npt),k=1,3)

       call nbondm()
       if (esymyes) call syminit()
       call eforce()


c    get the force projected along the path 
      istr = istr +1 
       do k = 1,3
          do m = 1,npick
           i = pointr(m)
               avg_fp(k,i)=avg_fp(k,i)+dpot(k,i) * pmn(k,m)

          end do
        end do


c      call deselect(dpot,pointr,npick)
c      write(99,*) ((coor(k,i),i=1,npt),k=1,3)
c      write(97,*) ((dpot(k,i),i=1,npt),k=1,3)
c      write(98,*) ((pmn(k,i),i=1,npick),k=1,3)

      write(upfr,*) projf(dpot,pmn,npick,pointr)*dS 
      enddo


c   report the average force 
c    get the force projected along the path 
       
      write(uvfr) e0,(((avg_fp(i,m)*dS/istr),m=1,npt), i=1,3) 

     
 1000 format(8(f12.6))
      
 2000 continue
      
      
      stop
      end

