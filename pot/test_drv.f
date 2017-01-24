       program test_drv
c
c a program to test current derivatives of the energy
c function useful for debugging
c
       
       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/CONSPECL1.BLOCK'
       include 'COMMON/CONSPECL2.BLOCK'
       include 'COMMON/UNITS.BLOCK'
       include 'COMMON/DEBUG.BLOCK'
       include 'COMMON/LINE.BLOCK'
       include 'COMMON/NBLIST.BLOCK'
       include 'COMMON/SYMM.BLOCK'
       include 'COMMON/COORD.BLOCK'
       include 'COMMON/ENERGY.BLOCK'
       include 'COMMON/SPECL.BLOCK'
       include 'COMMON/CONSTRAN.BLOCK'
       include 'COMMON/CONVERT.BLOCK'
       include 'COMMON/FREEZ.BLOCK'
       include 'COMMON/CCRD.BLOCK'
       include 'COMMON/VELOC.BLOCK'
       include 'COMMON/EWALD.BLOCK'
       include 'COMMON/METAL.BLOCK'
       include 'COMMON/SGB.BLOCK'
       include 'COMMON/EBALL.BLOCK'


       character*4 name,coortyp
       double precision getd
       integer namel,level,nstru
       integer of,geti
       integer i,j,k,n
       logical find

       integer ucon,ucor,uwene
       integer ipick(maxpt)
       integer ucon1,ucon2

       double precision diff1,diff2

       stdi = 5
       stdo = 6
       stderr = 0
       totmon = 0
       npt    = 0
       nb     = 0
       nmb    = 0
       nangl  = 0
       ntors  = 0
       nimp   = 0
       lestyp = 0
       nstru  = 1
       lpstr  = 1

       coortyp = 'CHAR'

       jnkf = 25
       open (unit=jnkf,status='scratch')

        lcent = .false.
       ctrue = .true. 
       shift = .false.
       name  = 'pote'
       namel = 4
       ucon  = -1
       ucor  = -1
       debug = .false.
       hydro_scale = 1.d0
       call init_ef()

1       continue

       call rline(name,namel,stdi)

       if (find('debu')) debug = .true.

       if (find('file')) then
         if (find('rcon').or.find('conn'))  then
         ucon  = of()
         call rconn(ucon)
c check for hydrophobic potential
         if (nbeta.gt.0) then
          ehyes = .true.
          hydro_th = 9999.d0
          hydro_scale = 1.d0
         end if
c initialze no freez vector
         inofrz = npt
         do 31 i=1,inofrz
       	zerofrz(i) = 1
31         continue
        end if
        if (find('con1')) ucon1 = of()
        if (find('con2')) ucon2 = of()
        if (find('rcrd')) then
         if (find('PATH')) then
          coortyp = 'PATH'
          ucor    = of()
         else if (find('DYNA')) then
          coortyp = 'DYNA'
          ucor    = of()
         else
          coortyp = 'CHAR'
          ucor  = of()
          call getcrd(ucor,'CHARM')
         end if
        end if
        if (find('wene')) uwene = of()
       end if

c---------------------------------------------
       if (find('mors')) then
        if (nmb.gt.maxmorsb) then
         level = 1
         call alert(name,namel,'Maxmorse exceeded',17,level)
        end if
       emyes0 = .true.
       do 55 n=1,nmb
        emyes(n)      = .true.
        D(n)     = getd('Dmor',D(n))
        alpha(n) = getd('alph',alpha(n))
55     continue
        end if

        if (find('spec')) then
          specl = .true.
         do 44 n=1,nmb
          rcut(n)   = getd('rcut',rcut(n))
          lamda(n) = getd('lmda',lamda(n))
44       continue
          call rcon_specl1(ucon1)
          call rcon_specl2(ucon2)
        end if
       if (find('repl')) then
        do 56 n=1,nmb
         repyes(n) = .true.
        Arep(n)  = getd('Arep',Arep(n))
        Brep(n)  = getd('Brep',Brep(n))
        beta1(n) = getd('beta',beta1(n))
56        continue
       end if
c---------------------------------------------

        if (find('gbsa')) then
                gbsabool=.true.
        end if
        if (find('gbo1')) then
           gbsabool=.true.
           gbobcbool=.true.
           gbalpha = 0.8d0
           gbbeta = 0.0d0
           gbgamma = 2.909125
        end if
        if (find('gbo2')) then
           gbsabool=.true.
           gbobcbool=.true.
           gbalpha = 1.0d0
           gbbeta = 0.8d0
           gbgamma = 4.85d0
        end if
        if (find('npol')) then
           gbnpbool=.true.
           surften=getd('sten',surften)
           call init_gb_nonpol
        end if
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

        gbsu    = geti('gbsu',gbsu)
        nstrub = geti('#stb',nstrub)

       nstru    = geti('#str',nstru)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
       cutmono2 = getd('cutm',cutmono2)
       rmax     = getd('rmax',rmax)
       eps    = (getd('epsi',(eps)))
       if (.not. ctrue) ctrue  = find('cdie')
       if (find('rdie')) ctrue = .false.
       hydro_scale = getd('hscl',hydro_scale)

       if (shift  .or.  find('shif')) shift   = .true.
       if (ebyes  .and. find('nobo')) ebyes   = .false.
       if (ethyes .and. find('noan')) ethyes  = .false.
       if (etoyes .and. find('noto')) etoyes  = .false.
       if (eimyes .and. find('noim')) eimyes  = .false.
       if (evdyes .and. find('novd')) evdyes  = .false.
       if (eelyes .and. find('noel')) eelyes  = .false.
       if (ecnyes .or.  find('cnst')) ecnyes  = .true.
       if ( find('symm')) then
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

c carlos, add amide plane constraints
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

       if (find('acti')) go to 2
       go to 1
2       continue

c rmax is maintained here for old input (with a single
c cutoff to work)
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if


        if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
        if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

        cutvdw2  = cutvdw2*cutvdw2
        cutvbig2 = cutvbig2*cutvbig2
        cutele2  = cutele2*cutele2
        cutebig2 = cutebig2*cutebig2

       if (cutmono2.lt.0) then
       	cutmono2 = cutebig2*1.44d0
       else
       	cutmono2 = cutmono2*cutmono2
       end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp



       if (debug) then
       	write(stdo,*)' after getcrd '
       	do 21 i=1,npt
       	 write(stdo,*)i,coor(1,i),coor(2,i),coor(3,i)
21       	continue
       end if
       if(nmb.gt.0 .and.(D(nmb).le.0. .or. alpha(nmb).le.0.)) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
        call init_wre(uwene)
        if (coortyp.eq.'CHAR') then
         call nbondm()
         if (esymyes) call syminit()
c =====================================================
        if(specl) call nbondm_spcl()
c ====================================================
         call num_drv()
         call eforce()
         write(stdo,*)' i ANALY dx dy dz NUM dx dy dx diff '
         do 22 i=1,npt
          diff1 = 0.d0
          diff2 = 0.d0
          do 23 j=1,3
           diff1 = diff1 + dpot(j,i)*dpot(j,i)
           diff2 = diff2 + coor2(j,i)*coor2(j,i)
23          continue
          diff1 = diff1-diff2
          if (abs(diff1).gt.1.d-4) then
           write(stdo,*)i,(dpot(k,i),coor2(k,i),k=1,3),diff1
          end if
22       continue
       end if

       stop
       end
