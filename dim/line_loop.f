       subroutine line_loop(dtopt,
     &           temper,constant_tmp,clo,interpolate,lap,
     &           select,pointr,ucon,urcrd,urint,urvel,amid_true,ufpt)
       implicit none
c
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/PATH2.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SDEL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/LD.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        double precision dtopt,temper, clo
        integer i,interpolate,j,iend,istart
        integer ucon,urcrd,urint,urvel,ufpt
        logical constant_tmp,lap,select,amid_true

C local variables        

        character*14 name
        integer namel,geti,of
        integer numremdgf, ncnsttmp
        integer*8 getL
        double precision getd
        integer ncns_start,ncns_end,tors(maxpt,4)
        logical find

        integer npick,pointr(maxpt) 

        name   = 'line_loop_sdel'
        namel  = 14

        FORCE_ENERGY = 1.d0
        nrigi = 10
        nsave = 10
        no_scaling = .false.
        ncoor = 0
        nscalv = 1.0
        LDgamma = 6.0
        ncns_start = 0
        ncns_end = 0
        Nreduced = 0    ! by default use all peptide RMSD
        select = .false.        

        cutoff=0.5
        dx = 0.1
        dx_umbrella = 0.05
        dx_approximate = 0.2
        langevin = .false.
        andersen = .false.
        coupling = 0.d0
        verbose = .false.
       
                        
c ....initialization of the parallel machine
        
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
           if (find('conn')) then
               ucon = of()
               call rconn(ucon)
               inofrz = npt
               iorie = npt
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

3              continue
               tgroup(1) = 3*npt
           endif
C ADD
           if (find('rcrd')) urcrd = of()
           if (find('rint')) urint = of()
           if (find('rvel')) urvel = of()
           if (find('wcrd')) uwcrd = of()
           if (find('wvel')) uwvel = of()
           if (find('wfpt')) ufpt = of()

C   read coarse grained model parameters (associated files)
           call Read_CG_input()

        end if

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

        dtopt=getd('dtop',dtopt)
        skpno   = geti('skno',skpno)

        if (find ('noRA'))  Random_velocities = .false.
        if (find('verb')) verbose = .true.
        
        if (find ('gbsa')) gbsabool=.true.
        gbsu=geti('gbsu',gbsu)
        nstep   = getL('#ste',nstep)
        ninfo   = geti('info',ninfo)
        ncoor   = geti('#pri',ncoor)
        igrid   = geti ('grid',igrid)
        myCell  = geti ('cell',myCell)
        myCell2 = geti ('cel2',myCell2)
        nStart = geti ('nInt',nStart)
        nEnd   = geti ('nEnd',nEnd)
        nsave   = geti('#sav',nsave)
        dt      = getd('step',dt)
        neqstep = geti('#equ',neqstep)
        
        if (find('nosc')) no_scaling = .true.
        
        if (find('mshk') .and. (.not.shakm)) then
c
c note that matrix shaking frozen particles is currently a bad idea
c and most likely will not work properly
c
          shakm  = .true.
          tolcons=getd('mtol',tolcons)
          write(stdo,1000)tolcons
1000      format (1x,'matrix shake on TIP3/SPC/SPCE, tolcons=',d7.1)
          call mshakinit(debug)
          if (nshakm.le.0) then
                write (stdo,*) 'nothing to shake!'
                shakm=.false.
          endif
          if (shakb.or.shakl) then
           write (stdo,*)'Re-calling shakinit
     1 to remove TIP3/SPC/SPCE from list'
           nshak=0
           call shakinit(shakl,shakb,shakm,epshak)
          endif
         endif

        if (find('shkl') .and. (.not.shakl)) then
          shakl  = .true.
          call shakinit(shakl,shakb,shakm,epshak)
        else if (find('shkb') .and. (.not.shakb)) then
          shakb  = .true.
          call shakinit(shakl,shakb,shakm,epshak)
        end if

c
c Set up a constraint on the geometric center of selected group
c
         if (find('cent')) then
           lcent = .true.
           kcenter=getd('kcnt',10.d0)
           xeq = getd('xeqm',0.d0)
           yeq = getd('yeqm',0.d0)
           zeq = getd('zeqm',0.d0)
           call pick(ipick,i)
           icenter=0

           do i=1,npt
             if (ipick(i).ne.0) then
               icenter = icenter+1
               center(icenter) = i
             endif
           end do

         endif

c shake parameters
        epshak   = getd('shac',epshak)
        epshakv  = getd('shav',epshakv)
        itershak = geti('itsh',itershak)

        cutoff  = getd('cutf',cutoff)
        dx = getd('delt',dx)
        dx_approximate = getd('appr',dx_approximate)
        dx_umbrella = getd('umbr',dx_umbrella)
        umbrella_K1 = getd('K1_U',umbrella_K1)
        umbrella_K2 = getd('K2_U',umbrella_K2)
        if (find('lang')) langevin=.true.
        if (find('andr')) andersen=.true.
        coupling = getd('andC',coupling)
        searchLimit  = getd('seaL',searchLimit)
        Stotal  = geti('stot',Stotal)
        ntemp = 1
        tempi(1) = getd('tmpi',tempi(1))
        tempf(1) = getd('tmpf',tempf(1))
        nwcrd   = geti('#wcr',nwcrd)
        temper  = getd('tmpr',temper)
        FORCE_ENERGY = getd('fene',FORCE_ENERGY)
        if (find ('ctmp')) constant_tmp = .true.
        LDgamma   = getd('gama',LDgamma)
        nscalv  = getd('#scl',nscalv)
        clo = getd('clog',clo)
        interpolate=(geti('itpl',(interpolate)))
        ENERGYC = (getd('pdqe',(ENERGYC)))
        if (find('ovlp')) lap = .true.
        if (find('sele')) select = .true.
  

c Energy parameters
        if (find('hvdw')) hvdw0 = .false.
        nlist   = geti('#lis',nlist)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        eps     = (getd('epsi',(eps)))

        if (find('cdie')) ctrue = .true.
        if (find('rdie')) ctrue = .false.

         irand   = geti('rand',irand)

        if (find('cpth')) crdstyl = 'PATH'
        if (find('cchr')) crdstyl = 'CHAR'
        if (find('cini')) crdstyl = 'INIT'
        if (find('cint')) crdstyl = 'INTR'

        if (find('nori')) nori = .true.
        if (find('nobo')) ebyes  = .false.
        if (find('noan')) ethyes = .false.
        if (find('noto')) etoyes = .false.
        if (find('noim')) eimyes = .false.
        if (find('novd')) evdyes = .false.
        if (find('noel')) eelyes = .false.
        if ( find('symm')) then
         esymyes = .true.
         a = getd('xtra',0.0d0)
         b = getd('ytra',0.0d0)
         c = getd('ztra',0.0d0)
        end if

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

       if (find('amid')) then
          amid_true=.true.
       endif
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
c compute pre-exponential  factor A0_metal such that
c A0_metal*exp(-alfa_metal*b/2)=a0_metal
c
c              a0_metal = a0_metal*dexp(-0.5d0*alfa_metal*b_wall)
                  write(stdo,*)
                  write(stdo,104)
104               format(1x,' Metal walls perpendicular',
     1          ' to the Y axis will be set!')
                  write(stdo,105)a0_metal
105               format(1x,' Metal repulsive Wall, A0 = ',
     1                  E12.6)
                  write(stdo,106)b_wall/2,v_elec
106               format(1x,' Walls are located at +/- ',f9.3,
     1                  1x,' The Electrode potential is ',f10.4)
       end if

c       read set of torsions to be used as reduced set of parameters       
        if (find('TORS')) then
           Nreduced = Nreduced + 1
           torsions(Nreduced,1) = geti('atm1',0)
           torsions(Nreduced,2) = geti('atm2',0)
           torsions(Nreduced,3) = geti('atm3',0)
           torsions(Nreduced,4) = geti('atm4',0)
           tor_weight(Nreduced)   = getd('weig',1.d0)
           tor_weight2(Nreduced) = tor_weight(Nreduced)**2
        end if


        if (find('acti')) go to 2
        go to 1
2       continue

        numremdgf=0
c
        if (select.and.Nreduced.eq.0) then 
        call rline(name,namel,stdi)
 
        npick=0
         call pick(ipick,i)
         do 31 i=1,npt
            npick=npick+ipick(i)
 31       continue
         j=0
         do 41 i=1,npt
            if (ipick(i).gt.0) then
               j=j+1
               pointr(j)=i
            end if
 41       continue
      endif
  
      

        if (ntemp.eq.1) then
         if ((.not.esymyes).and.(.not.qssbp).and.(.not.nori)
     >    .and.(.not.eteth_yes).and.(inofrz.eq.npt)) then
c
c we assume below that only diatomic has two rotations
c we are going to introduce error in strictly linear molecules
c which are not in our agenda anyway.
c
          if (npt.eq.2) then
                numremdgf = 5
          else
                numremdgf =6
          end if
         end if
        end if

c
c initiate GB calculations
c
      if (gbsabool) then
         call make_rborn
       endif


        call rm_degf(tpo,tempi,tempf,tgroup,ntemp,nofreez,
     1  inofrz,numremdgf)


101       format(1x,' # degrees of freedom in each temperature group '
     1,5(i7,1x))
        write(stdo,101)(tgroup(i),i=1,ntemp)
        if (prll_on_off .and. ncnst.gt.0) then
         ncnsttmp = ncns_end-ncns_start
         call load_balance(ncnsttmp,my_pe,num_pes,istart,iend)
         do 21 i=ncns_start+istart,ncns_start+iend
          j = i-istart+1
          icnst1(j) = icnst1(i)
          icnst2(j) = icnst2(i)
          icnst3(j) = icnst3(i)
          icnst4(j) = icnst4(i)
          kcns(j)   = kcns(i)
          cnseq(j)  = cnseq(i)
21       continue
        end if

        if (fmax.gt.0) sdyes = .true.

        return

        end
