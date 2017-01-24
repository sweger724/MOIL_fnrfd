        subroutine line_loop(ucon,urcrd,urvel)

c
c initializing as many as possible of the dyna variables by reading
c standard input lines
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/COUPLED_DYNA.BLOCK'

        character*6 getchar
        character*9 name
        integer namel,lc
        character*5 crdtyp
C Descriptor of input files for particular process of type fil1, fil2, ...
        character*4 fileID

        integer ucon,urcrd,urvel, ualign
        integer urtet
        integer ucon1,ucon2
        integer geti,of

        integer j,istart,iend,n,level

        double precision getd,eque
        logical find,equb,eqri

c a constant to assign new bond force constants
        double precision newb,tmp

        integer ncnsttmp

c ncns_start ncns_end - the start and the end of constraints
c need to be distributed between the processors
c
        integer ncns_start,ncns_end

        integer i,numremdgf

        name = 'line_loop'
        namel= 9
C fileID should look like fil1, fil2, etc.
        WRITE (fileID, '(A3,I1)'),'fil',procID+1 
        write(*,*) 'procID: ',procID

        lc = 5

        call init_var()

        gbsabool = .false.
        gbobcbool = .false.
        gbnpbool = .false.
        surften = 0.0005d0
        gbsu = 0
c
c open junk file for rline
c
        call open_scratch()



1       continue
        call rline(name,namel,stdi)
        if (find('debu')) then
         debug = .true.
        else if (find(fileID)) then
         if (find('alig')) then
                ualign = of()
                call ralign(ualign)
                close (ualign)
         end if
         if (find('conn')) then
c        get connectivity
                if (find('eqms')) eqms = .true.
                ucon  = of()
                call rconn(ucon)
                close (ucon)
                if (nbeta.gt.0) then
                 ehyes = .true.
                 hydro_th = 9999.d0
                 hydro_scale = 1.d0
                end if
                inofrz = npt
                iorie  = npt
c       initialize nofreez,jpick and tpo and inverse mass
c check if masses are set to be equal (for easier annealing)
                if (eqms) then
                 do 2 i=1,npt
                  ptms(i) = 10.d0
2                continue
                end if
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

3               continue
                tgroup(1) = 3*npt
         end if

         if (find('rcrd')) then
c get coordinate (currently from charmm format)
                urcrd = of()
                crdtyp = getchar('ctyp','CHARM',lc)
                if (crdtyp.eq.'PATH') then
                        call rpath(urcrd,1)
                else
                        call getcrd(urcrd,crdtyp(1:5))
                end if
                close (urcrd)
c xtwo is a copy of the original coordinate for comparison
c and overlap during the simulation. Used also in the etether
c
                do 4 i=1,npt
                        coor2(1,i) = coor(1,i)
                        coor2(2,i) = coor(2,i)
                        coor2(3,i) = coor(3,i)
4               continue
         end if
         if (find('rtet')) then
c get coordinate file for tether and overlap calculations
c they are read into the velocities and immediately transfered
c to xtwo ytwo ztwo
c
                
                urtet = of()
                crdtyp = getchar('ctyp','CHARM',lc)
                call getvel(urtet,crdtyp(1:5))
                close (urtet)
                do 5 i=1,npt
                        coor2(1,i) = velo(1,i)
                        coor2(2,i) = velo(2,i)
                        coor2(3,i) = velo(3,i)
5               continue
         end if
         if (find('wcrd')) uwcrd = of()
         if (find('wvel')) uwvel = of()
         if (find('rest')) urst  = of()
c ------------------------------------------
         if (find('con1')) ucon1 = of()
         if (find('con2')) ucon2 = of()
c ------------------------------------------
         if (find('rvel')) then
          boltz = .false.
          urvel = of()
         end if

c       rstr - A CRD file with recent coordinates for restart
         if (find('rstr')) then
          urst_crd = of()
         end if
         if (find('cont')) cont_yes = .true.

        
        else
c For restart - if true this run is a restart.
         if (find('cont')) then
                cont_yes = .true.
                call getcrd(urst_crd,'CHARM')
                close (urst_crd)
         end if

         if (find('bigb')) then
                write (stdo,*) 'changing bond force constants'
                newb=getd('newb',500.d0)
                call pick(ipick,i)
                do 6 i=1,nb
                   if (ipick(ib1(i)).ne.0.or.ipick(ib2(i)).ne.0) then
                        write (stdo,*) 'new kbond for ',ib1(i),ib2(i)
                        kbond(i)=newb
                   endif
6             continue
           end if
         nstep   = geti('#ste',nstep)
         neqstep = geti('#equ',neqstep)
         ninfo   = geti('info',ninfo)
         ncoor   = geti('#crd',ncoor)
         nvelo   = geti('#vel',nvelo)
         nlist   = geti('#lis',nlist)
         nscalv  = getd('#scl',nscalv)
         ntemp   = geti('#tmp',ntemp)
         nrigi   = geti('#rig',nrigi)
         irand   = geti('rand',irand)
         newv    = geti('newv',newv)
         dt      = (getd('step',(dt)))
c
         if (find('wfly')) then
            wfly = .true.
         endif


         if (find('gbsa')) then
                gbsabool = .true.
                gbsu = geti('gbsu',gbsu)
         end if
         if (find('gbo1')) then
                gbsabool = .true.
                gbobcbool= .true.
                gbalpha = 0.8d0
                gbbeta = 0.0d0
                gbgamma = 2.909125
                gbsu = geti('gbsu',gbsu)
         end if
         if (find('gbo2')) then
                gbsabool = .true.
                gbobcbool= .true.
                gbalpha = 1.0d0
                gbbeta = 0.8d0
                gbgamma = 4.85d0 
                gbsu = geti('gbsu',gbsu)
         end if
         if (find('npol')) then
             gbnpbool = .true.
             surften=getd('sten', surften)
             call init_gb_nonpol
         end if

         if (find('equb')) equb = .true.
         eque = getd('eque',eque)
         if (find('eqri')) eqri = .true.
         hydro_scale = getd('hscl',hydro_scale)
         fmax     = getd ('fmax',fmax)
         if (find('hvdw')) hvdw0 = .false.
         if (find('tstd')) test_d = .true.
c Brownian Dynamics dissipation
         LDgamma = getd('gama',LDgamma)
c a starting step for dynamics
c
         start_dyna = (geti('strt',start_dyna))
c
c a constraining sphere for water
c
        if (find('ball')) then
                eballyes = .true.
                fball = getd('fbal',0.d0)
                rball = getd('rbal',0.d0)
                rcball(1) = getd('xbal',0.d0)
                rcball(2) = getd('ybal',0.d0)
                rcball(3) = getd('zbal',0.d0)
                write(*,*) 'currnet values for fball rball rcball '
                write(*,*)fball,rball,rcball
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

           do 7 i=1,npt
             if (ipick(i).ne.0) then
               icenter = icenter+1
               center(icenter) = i
             endif
7         continue
c
c if parallel, spread the "center" terms between the pe-s
c
          if (prll_on_off) then
            call load_balance(icenter,my_pe,num_pes,istart,iend)
            do 71 i=istart,iend
                j = i - istart + 1
               center(j)=center(i)
71          continue
          end if
                
         endif
c
c Note that the selection is done for particles that WILL NOT
c               be freezed
c
         if (find('nfrz')) then
          call frz_eval()
        end if

c shake parameters
        epshak  = getd('shac',epshak)
        epshakv = getd('shav',epshakv)
        itershak= geti('itsh',itershak)
        cg_shak_it = geti('cgpt',cg_shak_it)
        cg_shakv_it = geti('cgvl',cg_shakv_it)
        if (find('cgsk')) matshak = .true.

        if (find('mshk') .and. (.not.shakm)) then
c
c note that matrix shaking frozen particles is currently a bad idea
c and most likely will not work properly
c
          shakm  = .true.
          tolcons=getd('mtol',tolcons)
          write(stdo,1000)tolcons
1000      format (1x,'matrix shake on TIP3, tolcons=',d7.1)
          call mshakinit(debug)
          if (nshakm.le.0) then
                write (stdo,*) 'nothing to shake!'
                shakm=.false.
          endif
          if (shakb.or.shakl) then
           write (stdo,*)'Re-calling shakinit to remove TIP3 from list'
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
         if (find('nori')) nori = .true.
         if (find('nosc')) no_scaling = .true.
         if (find('orie')) then
          call pick(ipick,i)
          iorie = 0
          do 9 i=1,npt
           if (ipick(i).ne.0) then
            iorie = iorie + 1
            jpick(iorie) = i
           end if
9         continue
         end if

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


c carlos, add amide plane constraints
         if (find('amid')) call amid()
                
         if (ecnyes .or.  find('cnst')) ecnyes  = .true.
c Note that the keyword TORS must come after amid
         ncns_start = ncnst
         if (find('TORS')) then
           write(*,*)
           write(*,*)'* Note that the keyword TORS must come after amid'
           write(*,*)'* If amid is used'
           write(*,*)
           ncnst = ncnst + 1
           ncns_end = ncnst
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


c energy parameters
          cutvdw2  = (getd('rvmx',(cutvdw2)))
          cutvbig2 = (getd('rvbg',cutvbig2))
          cutele2  = (getd('relx',(cutele2)))
          cutebig2 = (getd('rebg',cutebig2))
          cutmono2 = getd('cutm',cutmono2)
          rmax     = getd('rmax',rmax)
          hydro_scale = getd('hscl',hydro_scale)
          hydro_th    = getd('hthr',hydro_th)
  
         eps    = (getd('epsi',(eps)))
         if (find('cdie')) ctrue = .true.
         if (find('rdie')) ctrue = .false. 

c shake parameters
         epshak   = getd('shac',epshak)
         epshakv  = getd('shav',epshakv)
         itershak = geti('itsh',itershak)

c Morse parameters
c
        if (nmb.gt. 0 ) then
         if (nmb.gt.maxmorsb) then
          level = 1
          call alert(name,namel,'Maxmorse exceeded',17,level)
         end if
         emyes0    = .true.

         do 13 n=1,nmb
          emyes(n) = .true.
          D(n)     = getd('Dmor',D(n))
          alpha(n) = getd('alph',alpha(n))
13       continue
        end if

C SPEC IS NOT SET UP FOR PARALLELIZATION (YET)
        if (find('spec')) then
         specl = .true.
         call rcon_specl1(ucon1)
         call rcon_specl2(ucon2)
         do 14 n=1,nmb
         rcut(n)   = getd('rcut',rcut(n)) 
         lamda(n)  = getd('lmda',lamda(n))
14       continue
        end if

C REPL IS NOT SET UP FOR PARALLELIZATION (YET)
      if (find('repl')) then
       repyes0 = .true.
         do 15 n=1,nmb
          repyes(n) = .true.
          Arep(n)   = getd ('Arep',Arep(n))
          Brep(n)   = getd ('Brep',Brep(n))
          beta1(n)  = getd('beta',beta1(n))
15       continue
        endif
C SWITCH IS NOT SET UP FOR PARALLELIZATION (YET)        
        if ( find('swit')) switch = .true.
        delt  = getd ('delt',delt )
        Rcros = getd ('Rcro',Rcros)
        dRcros = getd ('dRcr',dRcros)
        Force = getd ('Forc',Force)
        Rcros1 = Rcros - dRcros/2.d0
        Rcros2 = Rcros + dRcros/2.d0
        if(switch) then
           
           if(Rcros.le.0. .or. dRcros.le.0.) then
             level=1
             call alert(name,namel,'crosing point Rcro or
     $   mixing region is 0',16,level)
           end if

        end if

c---------------------------------------------


         if (find('nocut')) nocut = .true.
         if (shift  .or.  find('shif')) shift   = .true.
         if (ebyes  .and. find('nobo')) ebyes   = .false.
         if (ethyes .and. find('noan')) ethyes  = .false.
         if (etoyes .and. find('noto')) etoyes  = .false.
         if (eimyes .and. find('noim')) eimyes  = .false.
         if (evdyes .and. find('novd')) evdyes  = .false.
         if (eelyes .and. find('noel')) eelyes  = .false.
         if (ehyes  .and. find('nohy')) ehyes   = .false.
         if (find('nbfi')) then
                efyes = .true.
                write(stdo,*)
                write(stdo,*) 'FINITE VAN DER WAALS AT ZERO DISTANCE'
                write(stdo,*)
         end if
         if (find('symm')) then
          esymyes = .true.
          a = getd('xtra',0.0d0)
          b = getd('ytra',0.0d0)
          c = getd('ztra',0.0d0)
         end if
       if (find('sym2')) then
          symanneal=.true.
          syma2 = getd('xtr2',0.0d0)
          symb2 = getd('ytr2',0.0d0)
          symc2 = getd('ztr2',0.0d0)
          write (stdo,*)
     &  'Symmetry box size equilibration will be performed'
        end if


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
            if (.not.esymyes) then
               write (stdo,*)
     &  'There is no true periodic bound. cond. - symm is missing'
               a = getd('xtra',50.0d0)
               b = getd('ytra',50.0d0)
               c = getd('ztra',50.0d0)
            end if
         end if
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




         if (find('tthr')) then
                eteth_yes = .true.
                tmp       = getd('frcc',1.d0)
                call pick(ipick,i)
                do 16 i=1,npt
                 if (ipick(i).ne.0) then
                        n_tether = n_tether + 1
                        tether(n_tether) = i
                        tf(n_tether) = tmp
                 end if
16        continue
          write(stdo,*)' total tether # = ',n_tether
          if (prll_on_off) then
           call load_balance(n_tether,my_pe,num_pes,istart,iend)
           do 17 i=1,n_tether
                j = i+istart-1
                tether(i) = j
                tf(i)     = tmp
17         continue 
          end if
         end if

         if (find('ssbp')) then
                write(*,*)' Calling init_ssbp '
                call init_ssbp(ipick)   
         end if
 
 
         if (find('mult')) then
c
c If 'mult' then multiple temperatures are present
c tpo is an integer vector with a length equal to number
c of particles which points to the tempi/tempf arrays
c Note that if les particles are present but mult was not found,
c no les conscious scaling will be performed. I.e. all
c particles will be scaled with the same factor.
c
          call picktemp(tpo,tempi,tempf,tgroup,ntemp,nofreez,inofrz)
          write(stdo,100)ntemp
          write(stdo,102)(tempi(i),i=1,ntemp)
          write(stdo,103)(tempf(i),i=1,ntemp)
100       format(1x,' Number of different temperatures:',i5)
101       format(1x,' # degrees of freedom in each temperature group '
     1          ,5(i7,1x))
102       format(1x,' Initial temperatures: ',5(f10.4,1x))
103       format(1x,' Final   temperatures: ',5(f10.4,1x))

         else
          tempi(1) = getd('tmpi',tempi(1))
          tempf(1) = getd('tmpf',tempf(1))
         end if
        end if
        if (find('acti')) go to 20
        go to 1

20      continue
        
        numremdgf=0
c
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

c
c if parallel is on divide the particles between the processors
c
c               if (prll_on_off.and.inofrz.eq.npt) then
c                npt_par = npt
c                call load_balance(inofrz,my_pe,num_pes,istart,iend)
c                do 31 i=istart,iend
c                 j          = i - istart + 1
c                 nofreez(j) = i
c31              continue
c               end if
        
        return
        end
