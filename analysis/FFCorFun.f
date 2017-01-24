       program FFCorFun

       implicit none

c common blocks

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
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
        include 'COMMON/LD.BLOCK'
        include 'COMMON/PT.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PRESS.BLOCK'

c declare the variables you need
c loops
        integer i,i1,i2,istr,nstru
c I/O input
        integer ucon,urcrd,urvel,of,namel,rbin
        integer unitFF
        integer MAXnGAMMA,MAXnstru,MAXNCOR,NCOR
        parameter(MAXnGAMMA = 20)
        parameter(MAXnstru = 100000)
        parameter(MAXNCOR = 1000)
        integer nGAMMA,iGAMMA(MAXnGAMMA),iNOEF(maxpt)
        logical NOEFBO
        double precision FX(MAXnGAMMA,MAXnstru),
     1                   FY(MAXnGAMMA,MAXnstru),
     2                   FZ(MAXnGAMMA,MAXnstru)


        double precision FXi1(MAXnstru),FYi1(MAXnstru),FZi1(MAXnstru)
        double precision FXi2(MAXnstru),FYi2(MAXnstru),FZi2(MAXnstru)

        double precision CFXX(MAXNCOR),CFXY(MAXNCOR),CFXZ(MAXNCOR),
     1                   CFYX(MAXNCOR),CFYY(MAXNCOR),CFYZ(MAXNCOR),
     2                   CFZX(MAXNCOR),CFZY(MAXNCOR),CFZZ(MAXNCOR)
        double precision intCFXX,intCFXY,intCFXZ,
     1                   intCFYX,intCFYY,intCFYZ,
     2                   intCFZX,intCFZY,intCFZZ
        double precision CSI(3*MAXnGAMMA,3*MAXnGAMMA)

        logical find
        character*8 name

c initialization

        stdi = 5
        stdo = 6
        rbin = 1
 
        name = 'FFCorFun'
        namel = 8
        NOEFBO = .FALSE.

        call readinput(ucon,urcrd,urvel,name,namel,
     1  nstru,NCOR,dt,iNOEF,NOEFBO)

c remove all the interactions with a selected subset of the system
        DO i=1,npt
         if (iNOEF(i).eq.1) then
          epsgm12(i) = 0.d0
          epsgm6(i)  = 0.d0
          ptchg(i)   = 0.d0
         endif
        ENDDO
        if (NOEFBO) write (*,*) 'BONDED INTERACTIONS WILL NOT BE USED'
        DO i=1,nb
         if (iNOEF(ib1(i)).eq.1 .or. iNOEF(ib2(i)).eq.1 .or. 
     1        NOEFBO) then
          kbond(i) = 0.d0
         endif
        ENDDO
        DO i=1,nangl
         if (iNOEF(iangl1(i)).eq.1 .or. iNOEF(iangl2(i)).eq.1 .or.
     1       iNOEF(iangl3(i)).eq.1 .or.
     2        NOEFBO) then
          kangl(i) = 0.d0
         endif
        ENDDO
        DO i=1,ntors
         if (iNOEF(itor1(i)).eq.1 .or. iNOEF(itor2(i)).eq.1 .or.
     1       iNOEF(itor3(i)).eq.1 .or. iNOEF(itor4(i)).eq.1 .or.
     2        NOEFBO) then
          ktors1(i) = 0.d0
          ktors2(i) = 0.d0
          ktors3(i) = 0.d0
         endif
        ENDDO
        DO i=1,nimp
         if (iNOEF(iimp1(i)).eq.1 .or. iNOEF(iimp2(i)).eq.1 .or.
     1       iNOEF(iimp3(i)).eq.1 .or. iNOEF(iimp4(i)).eq.1 .or.
     2        NOEFBO) then
          kimp(i) = 0.d0
         endif
        ENDDO
        DO i=1,totspe
         if (iNOEF(spec1(i)).eq.1 .or. iNOEF(spec2(i)).eq.1) then
          p14(1,i) = 0.d0
          p14(2,i) = 0.d0
          p14(3,i) = 0.d0
         endif
        ENDDO

        inofrz = npt

        DO i=1,inofrz
         nofreez(i) = i
        ENDDO

        nGAMMA = 0

        DO i=1,npt
         if (ptnm(i).eq.'CA') then
          nGAMMA = nGAMMA + 1
          iGAMMA(nGAMMA) = i
          write (*,*) nGAMMA,iGAMMA(nGAMMA)
         endif
        ENDDO

        DO 1 istr=1,nstru

         rewind urcrd
         call rdyncrd(urcrd,istr,inofrz,nofreez,rbin)

         if (esymyes) call squeeze(a,b,c)
         call nbondm()
         if (esymyes) call syminit()
         if (specl) call nbondm_spcl()

         call eforce()
         call wener(stdo)

         DO i=1,nGAMMA
          write (*,*) 'CA',i,iGAMMA(i),coor(1,iGAMMA(i)),
     1                       coor(2,iGAMMA(i)),coor(3,iGAMMA(i))
          FX(i,istr) = -dpot(1,iGAMMA(i))
          FY(i,istr) = -dpot(2,iGAMMA(i))
          FZ(i,istr) = -dpot(3,iGAMMA(i))
         ENDDO

1       CONTINUE

        DO i1=1,nGAMMA
         DO i2=1,nGAMMA
          DO istr=1,nstru
           FXi1(istr) = FX(i1,istr)
           FYi1(istr) = FY(i1,istr)
           FZi1(istr) = FZ(i1,istr)
           FXi2(istr) = FX(i2,istr)
           FYi2(istr) = FY(i2,istr)
           FZi2(istr) = FZ(i2,istr)
          ENDDO
c          DO i=1,nstru
c           FXi2(i) = FXi1(i)
c          ENDDO 
          call CORFUN(NCOR,nstru,dt,FXi1,FXi2,CFXX,intCFXX)
          call CORFUN(NCOR,nstru,dt,FXi1,FYi2,CFXY,intCFXY)
          call CORFUN(NCOR,nstru,dt,FXi1,FZi2,CFXZ,intCFXZ)
          call CORFUN(NCOR,nstru,dt,FYi1,FXi2,CFYX,intCFYX)
          call CORFUN(NCOR,nstru,dt,FYi1,FYi2,CFYY,intCFYY)
          call CORFUN(NCOR,nstru,dt,FYi1,FZi2,CFYZ,intCFYZ)
          call CORFUN(NCOR,nstru,dt,FZi1,FXi2,CFZX,intCFZX)
          call CORFUN(NCOR,nstru,dt,FZi1,FYi2,CFZY,intCFZY)
          call CORFUN(NCOR,nstru,dt,FZi1,FZi2,CFZZ,intCFZZ)
          unitFF = 100+(i1-1)*nGAMMA+i2
          write (unitFF,'(1i5,9f10.5)') -1,
     1                                  intCFXX,intCFXY,intCFXZ,
     2                                  intCFYX,intCFYY,intCFYZ,
     3                                  intCFZX,intCFZY,intCFZZ
          DO i=1,NCOR
           write (unitFF,'(1i5,9f10.5)') 
     1                                i,CFXX(i),CFXY(i),CFXZ(i),
     2                                  CFYX(i),CFYY(i),CFYZ(i),
     3                                  CFZX(i),CFZY(i),CFZZ(i)
          ENDDO
          CSI((i1-1)*3+1,(i2-1)*3+1) = intCFXX
          CSI((i1-1)*3+1,(i2-1)*3+2) = intCFXY
          CSI((i1-1)*3+1,(i2-1)*3+3) = intCFXZ
          CSI((i1-1)*3+2,(i2-1)*3+1) = intCFYX
          CSI((i1-1)*3+2,(i2-1)*3+2) = intCFYY
          CSI((i1-1)*3+2,(i2-1)*3+3) = intCFYZ
          CSI((i1-1)*3+3,(i2-1)*3+1) = intCFZX
          CSI((i1-1)*3+3,(i2-1)*3+2) = intCFZY
          CSI((i1-1)*3+3,(i2-1)*3+3) = intCFZZ
         ENDDO
        ENDDO

      stop
      end

        subroutine readinput(ucon,urcrd,urvel,name,namel,
     1  nstru,NCOR,dt,iNOEF,NOEFBO)

        implicit none

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
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'

        character*8 name
        character*4 coortyp
        double precision getd,dt
        integer namel,level,nstru
        integer of,geti
        integer i,n,nstrub
        integer rbin
        logical find

        integer ucon,urcrd,urvel
        integer ipick(maxpt)
        integer ucon1,ucon2

        integer NCOR,iNOEF(maxpt)
        logical NOEFBO
 
c        integer MAXnGAMMA
c        integer nGAMMA,iGAMMA(MAXnGAMMA)

        stdi = 5
        stdo = 6
        stderr = 0
        rbin = 1

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

        my_pe = 0
        num_pes = 1

        coortyp = 'CHAR'
        jnkf = 25
        open (unit=jnkf,status='scratch')

        lcent = .false.
        ctrue = .true.
        shift = .false.
        ucon  = -1
        urcrd  = -1
        urvel  = -1
        debug = .false.
        hydro_scale = 1.d0
        surften=0.005d0

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
          npt_par= npt
          do 31 i=1,inofrz
                zerofrz(i) = 1
                prtc_pointer(i) = i
31        continue
         end if
         if (find('rcrd')) urcrd = of()
         if (find('rvel')) urvel = of()

C   read coarse grained model parameters (associated files)
         call Read_CG_input()

        end if

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

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
56       continue
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

        dt = getd('step',0.001)
        NCOR = geti('NCOR',100)

        gbsu    = geti('gbsu',gbsu)
C let nstru use old style #str and new style
        nstru    = geti('#str',nstru)
        nstru = geti('#ste',nstru)
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

c repulsive wall
         if (find('rwal')) then
          rwall = .TRUE.
          nwalls = geti('#wal',nwalls)
          w0(1) = getd('wpo1',w0(1))
          normw(1) = geti('nwa1',normw(1))
          w0(2) = getd('wpo2',w0(2))
          normw(2) = geti('nwa2',normw(2))
          w0(3) = getd('wpo3',w0(3))
          normw(3) = geti('nwa3',normw(3))
          w0(4) = getd('wpo4',w0(4))
          normw(4) = geti('nwa4',normw(4))
          w0(5) = getd('wpo5',w0(5))
          normw(5) = geti('nwa5',normw(5))
          w0(6) = getd('wpo6',w0(6))
          normw(6) = geti('nwa6',normw(6))
          weps = getd('weps',weps)
                  write(stdo,*) 'rwal = ',rwall
          do i=1,nwalls
           if (normw(i).eq.1) then
            write (stdo,1045) i
1045        format(1x,' Repulsive wall',1i5,' perpedicular to x axis')
           else if (normw(i).eq.2) then
            write (stdo,1046) i
1046        format(1x,' Repulsive wall',1i5,' perpedicular to y axis')
           else if (normw(i).eq.3) then
            write (stdo,1047) i
1047        format(1x,' Repulsive wall',1i5,' perpedicular to z axis')
           else
            call alert(name,namel,'wall out of 3 dimensions!',25,1)
           endif
          enddo
          do i=1,nwalls
           write (stdo,1048) i,normw(i),w0(i)
1048       format(1x, 'wall #',1i5,'is perpendicular to axis',1i5,
     1    'and located at ',f9.3)
           if (normw(i).eq.0)
     1      call alert(name,namel,'no norm for a wall!',19,1)
          enddo
          write (stdo,1049) weps
1049      format(1x, 'rep wall is weps/r^6 --- weps= ',f9.3)
          endif

c electric field
          if (find('efie')) then
           efield_yes = .TRUE.
           DV = getd('DelV',DV)
           efield = getd('elec',efield)
           nefield = geti('nefi',nefield)
           DVtype = geti('DVty',DVtype)
           if (efield.ne.0.d0 .and. DV .ne. 0.d0) then
            call alert(name,namel,'both E and DV defined!',22,1)
           else if (efield.eq.0.d0 .and. DV.eq.0.d0) then
            call alert(name,namel,'neither E nor DV defined!',25,1)
           endif
          endif

c take some particles out from the computation of energy
        if (find('noef')) then
         if (find('nobo')) NOEFBO = .true.
         call pick(ipick,i)
         do i=1,npt
          iNOEF(i) = ipick(i)
         enddo
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

        if (debug) then
                write(stdo,*)' after getcrd '
                do 21 i=1,npt
                 write(stdo,*)i,coor(1,i),coor(2,i),coor(3,i)
21              continue
        end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp

         if (eCGyes) then
           call CGinit()
         endif

        if (gbsabool) call make_rborn
        if(nmb.gt.0) then
         if (D(nmb).le.0.d0 .or. alpha(nmb).le.0.d0) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
        end if
 
        return
        end 

        subroutine checkANDmove(fix,move,coor,a,b,c)
        implicit none
        integer fix,move
        double precision coor(3,*)
        double precision a,b,c
c check X
        if ((coor(1,fix)-coor(1,move)).gt.0.5d0*a)
     1                   coor(1,move) = coor(1,move) + a
        if ((coor(1,fix)-coor(1,move)).lt.-0.5d0*a)
     1                   coor(1,move) = coor(1,move) - a
c check Y
        if ((coor(2,fix)-coor(2,move)).gt.0.5d0*b)
     1                   coor(2,move) = coor(2,move) + b
        if ((coor(2,fix)-coor(2,move)).lt.-0.5d0*b)
     1                   coor(2,move) = coor(2,move) - b
c check Z
        if ((coor(3,fix)-coor(3,move)).gt.0.5d0*c)
     1                   coor(3,move) = coor(3,move) + c
        if ((coor(3,fix)-coor(3,move)).lt.-0.5d0*c)
     1                   coor(3,move) = coor(3,move) - c
        return
        end
