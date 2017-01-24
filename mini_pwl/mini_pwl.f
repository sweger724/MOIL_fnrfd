        program mini_pwl
        implicit none   
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
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/TETHER.BLOCK'

        character*8 name
        character*4 cstyl
        double precision getd
        integer namel
        integer of
        integer ucon1,ucon2
        integer geti
        integer rbin,wbin
        logical find

        double precision tolf,tmp
        integer ucon,ucor_in,ucor_out,umini,mistep,nlist
        integer i,istru,ig,istep,efcall,npt3,level,ncall,j
        integer ipick(maxpt)
        integer lpst,lpend,idyna,n
        integer istart,iend

        logical write_path


c setp-up for the powell routine.
c this part will be removed later for more efficient space
c handling
        double precision estred

        stdi = 5
        stdo = 6
        stderr = 0
        rbin = 1
        wbin = 1

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
        write_path = .false.


        jnkf = 25
        open (unit=jnkf,status='scratch')

        surften = 0.005d0

c initialized file units and debug
        name      = 'mini_pwl'
        namel     = 8
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
        norew  = .true.


c minimization default parameters
c tolf  - allowed error (TOLerance) in the Force Norm
c mistep - maximum number of MInimization STEPs
        tolf   = 0.d0
        mistep = 100
        nlist  = 20
        lpend = -1

1       continue

        call rline(name,namel,stdi)

        if (find('debu')) debug = .true.

        if (find('file')) then

         if (find('rcon').or.find('conn'))  then
          ucon    = of()
          call rconn(ucon)
          if (nbeta.gt.0) then
              ehyes = .true.
              hydro_th = 9999.d0
              hydro_scale = 1.d0
          end if
          inofrz = npt
c       initialize nofreez, zerofrz
          do 2 i=1,npt
                nofreez(i) = i
                zerofrz(i) = 1
2         continue

          if (prll_on_off) then
c
c divide the particles between the pe-s
c
           npt_par = npt
           call load_balance(npt_par,my_pe,num_pes,istart,iend)
           do 3 i=istart,iend
                        j               = i - istart + 1
                        prtc_pointer(j) = i
3          continue

          else

           npt_par = npt
           istart  = 1
           iend    = npt
           do 4 i=1,npt
            prtc_pointer(i) = i
4          continue

          end if

         endif
         if (find('con1')) ucon1   = of()
         if (find('con2')) ucon2   = of()
         if (find('rcrd')) then

        if (find('cpth')) then
                cstyl = 'PATH'
                istru = geti('istr',0)
        else if (find('cdyn')) then
                cstyl = 'DYNA'
                istru = geti('istr',0)
        else if (find('DYNA')) then
                cstyl = 'DYNA'
                lpst  = geti('lpst',1)
                lpend = geti('lpen',-1)
                istru = 0
        else if (find('PATH')) then
                cstyl = 'PATH'
                lpst  = geti('lpst',1)
                lpend = geti('lpen',-1)
                istru = 0
        end if

            ucor_in = of()
            if (cstyl.eq.'CHAR') then
                 call getcrd(ucor_in,'CHARM')
            end if
            rewind ucor_in
         end if
         if (find('wcrd')) ucor_out= of()
         if (find('wpth')) then
            ucor_out= of()
            write_path = .true.
         end if
         if (find('wmin')) umini   = of()

C   read coarse grained model parameters (associated files)
         call Read_CG_input()

        end if

c---------------------------------------------

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

        if (find('mors')) then
         if (nmb.gt.maxmorsb) then
          level = 1
          call alert(name,namel,'Maxmorse exceeded',17,level)
         end if
         emyes0     = .true.
        do 5 n=1,nmb
         D(n)     = getd('Dmor',D(n))
         alpha(n) = getd('alph',alpha(n))
         emyes(n) = .true.
5       continue
        end if
          if (find('repl')) then
           repyes0 = .true.
           do 6 n=1,nmb
           repyes(n)= .true.
           Arep(n)  = getd('Arep',Arep(n))
           Brep(n)  = getd('Brep',Brep(n))
           beta1(n) = getd('beta',beta1(n))
6          continue
          end if
c---------------------------------------------
         if (find('cent')) then
           i=1
           lcent = .true.
           kcenter = getd('kcnt',10.d0)
           xeq = getd('xeqm',0.d0)
           yeq = getd('yeqm',0.d0)
           zeq = getd('zeqm',0.d0)
           call pick(ipick,i)
          icenter= 0
          if (.not. prll_on_off) then
                istart = 1
                iend   = npt
          end if
          do 7 i=istart,iend
           if (ipick(i).ne.0) then
            icenter=icenter + 1
            center(icenter) = i
           end if
7       continue
         endif

         if (find('tthr')) then
                eteth_yes = .true.
                tmp       = getd('frcc',1.d0)
                call pick(ipick,i)
                do 161 i=1,npt
                 if (ipick(i).ne.0) then
                        n_tether = n_tether + 1
                        tether(n_tether) = i
                        tf(n_tether) = tmp
                 end if
161        continue
          write(stdo,*)' total tether # = ',n_tether
          if (prll_on_off) then
           call load_balance(n_tether,my_pe,num_pes,istart,iend)
           do 171 i=1,n_tether
                j = i+istart-1
                tether(i) = j
                tf(i)     = tmp
171         continue
          end if
                do 41 i=1,npt
                        coor2(1,i) = coor(1,i)
                        coor2(2,i) = coor(2,i)
                        coor2(3,i) = coor(3,i)
41               continue
         end if


        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        if (find('gbsa')) gbsabool = .true.
c YS
        if (find('gbo1')) then
            gbsabool=.true.
            gbobcbool=.true.
            gbalpha = 0.8D0
            gbbeta = 0.0D0
            gbgamma = 2.909125
        end if
        if (find('gbo2')) then
            gbsabool=.true.
            gbobcbool=.true.
            gbalpha = 1.0D0
            gbbeta = 0.8D0
            gbgamma = 4.85D0 
        end if
        if (find('npol')) then
           gbnpbool=.true.
           surften=getd('sten',surften)
           call init_gb_nonpol
        end if

        gbsu     = geti('gbsu',gbsu)
        eps    = (getd('epsi',(eps)))
        if (.not. ctrue) ctrue  = find('cdie')
        if (find('nbfi')) efyes = .true.
        if (find('rdie')) ctrue = .false.
        hydro_scale = getd('hscl',hydro_scale)
        if (ebyes  .and. find('nobo')) ebyes   = .false.
        if (ethyes .and. find('noan')) ethyes  = .false.
        if (etoyes .and. find('noto')) etoyes  = .false.
        if (eimyes .and. find('noim')) eimyes  = .false.
        if (evdyes .and. find('novd')) evdyes  = .false.
        if (eelyes .and. find('noel')) eelyes  = .false.
        if (ehyes  .and. find('nohy')) ehyes   = .false.

        if (find('hvdw')) hvdw0 = .false.

        if (nocut  .and. find('cute')) nocut   = .false.
        if (ecnyes .or.  find('cnst')) ecnyes  = .true.
        if (shift  .or.  find('shif')) shift   = .true.
        if (find('symm')) then
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
                  write(stdo,114)
114               format(1x,' Metal walls perpendicular',
     1          ' to the Y axis will be set!')
                  write(stdo,115)a0_metal
115               format(1x,' Metal repulsive Wall, A0 = ',
     1                  E12.6)
                  write(stdo,116)b_wall/2,v_elec
116               format(1x,' Walls are located at +/- ',f9.3,
     1                  1x,' The Electrode potential is ',f10.4)
       end if

        tolf   = getd('tolf',tolf)
        mistep = geti('mist',mistep)
        mistep = geti('#ste',mistep)
        nlist  = geti('list',nlist)
        nlist  = geti('#lis',nlist)

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
     1          /(getd('step',1.d0)/pi180) + 0.5d0
           end if
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


        if (find('acti')) go to 8
        go to 1
8       continue

c If number of constraints is greater than zero and parallel processing
c is on distribute the constraints between the processors
c
        if (ncnst.gt.0 .and. prll_on_off) then
         call load_balance(ncnst,my_pe,num_pes,istart,iend)
          if (my_pe.ne.0) then
            do 9 i=1,ncnst
                j = i - istart + 1
                icnst1(i) = icnst1(j)
                icnst2(i) = icnst2(j)
                icnst3(i) = icnst3(j)
                icnst4(i) = icnst4(j)
                kcns(i)   = kcns(j)
                cnseq(i)  = cnseq(j)
9           continue
           end if
        end if  

c rmax is maintained here for old input (with a single
c cutoff) to work
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if

        if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
        if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

        if (gbsabool) call make_rborn

        if(mistep .lt. nlist) nlist = mistep

        write(stdo,100)mistep,nlist,ucon,ucor_in,ucor_out,
     1          (cutvdw2),(cutele2),tolf
100     format(1x,' Minimization parameters ',/,
     1  ' Maximum number of minimization steps ',i6,/,
     2  ' Update non-bond list each ',i6, ' steps',/,
     3  ' Connectivity file to be read from unit ',i5,/,
     4  ' Initial coordinates to be read from unit ',i5,/,
     5  ' Final coordinates will be written to unit ',i5,/,
     6  ' The cutoff distance for van der Waal is ',f10.5,/,
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
         do 10 i=1,nvar
          var2(i) = (var2(i)-var1(i))/dfloat(varstep(i))
10       continue
        end if
        tolf   = tolf*tolf


        if (cstyl.eq.'PATH') then
         if (istru.gt.0) then
c@
          write(*,*)' istru = ',istru
          call rpath(ucor_in,istru)
         else if (lpend.eq.-1) then
          level = 1
          call alert(name,namel,'Illegal structure number',24,level)
         end if
        else if (cstyl.eq.'DYNA') then
         if (istru.gt.0) then
          norew  = .false.
          call rdyncrd(ucor_in,istru,inofrz,nofreez,rbin)
         else if (lpend.eq.-1) then
          level = 1
          call alert(name,namel,'Illegal structure number',24,level)
         end if
        end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp
        if (eCGyes) then
          call CGinit()
        endif


c *** HERE START A SINGLE MINIMIZATION

        npt3   = 3*npt
        estred = 0.001d0*npt3


        if (nvar.eq.0 .and. lpend.eq.-1) then

        write(umini,101)0
101     format(1x,' At step ',i5,/)
        call init_wre(umini)
        if (emyes0) then
         do 11 n=1,nmb
          if(alpha(n) .le.0. .or. D(n).le.0.) then
           level=1
           call alert(name,namel,'Dmor or alph is 0',16,level)
          end if
          write (stdo,*)'D=',D(n)
          write (stdo,*)' alpha= ',alpha(n)
11       continue
        end if

        do 13 istep=1,mistep,nlist
                write (*,*) 'mini step=',istep
                call nbondm()
                if (esymyes) call syminit()
                ncall = nlist
                if (ncall + istep - 1 .gt. mistep)
     >              ncall = mistep - istep + 1
                call powel(tolf,ncall,estred,efcall)
                if (prll_on_off) call reduce_energies()
                if (my_pe.eq.0) then
                 write(umini,102)istep+efcall-1
102              format(1x,'After ',i7,' minimization steps ')
                end if
                tmp = 0.d0
                do 12 i=1,npt_par
                  ig = prtc_pointer(i)
                  tmp=tmp+dpot(1,ig)*dpot(1,ig)+dpot(2,ig)*dpot(2,ig)
     6                  +dpot(3,ig)*dpot(3,ig)
12               continue
                 if (prll_on_off) call reduce_1(tmp)
                 tmp = dsqrt(tmp/npt3)
                 if (my_pe.eq.0) then
                  write(umini,103)tmp,dsqrt(tolf)
103               format(1x,' Current gradient norm ',f10.5,
     1' requested ',f10.5)
                 end if
                 if (tmp*tmp.lt.tolf) then
                  call eforce()
                  if (prll_on_off) call reduce_energies()
                  if (my_pe.eq.0) then
                   write(umini,104)
104                format(1x,' Minimization converged ')
                   call wener(umini)
                  end if
                  go to 14
                 end if
c get "best value" of energy and printout
                 call eforce()
                 if (prll_on_off) call reduce_energies()
                 if (my_pe.eq.0) call wener(umini)
13      continue
14      continue
        if (write_path .and. (my_pe.eq.0)) then
                call wpath(ucor_out)
        else if (my_pe.eq.0) then
                call putcrd(ucor_out,'CHARM')
        end if

        write(stdo,105)
105     format(1x,/,1x,' *** Minimization completed ')

c *** HERE END A SINGLE MINIMIZATION AND DYNAMICS LOOP MINIMIZATION STARTS
c ***

        else if (nvar.eq.0 .and. lpend.ne.-1) then
        write(umini,106)
106     format(1x,' Structure number   total energy   gradient norm ')

        if (emyes0) then
         if(alpha(nmb) .le.0. .or. D(nmb).le.0.) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
         write (stdo,*)'D=',D(nmb)
         write (stdo,*)'alpha=',alpha(nmb)
        end if

        do 18 idyna=lpst,lpend
        if (cstyl.eq.'DYNA') then
                call rdyncrd(ucor_in,idyna,inofrz,nofreez,rbin)
        else if (cstyl.eq.'PATH') then
                call rpath(ucor_in,idyna)
        end if
        do 16 istep=1,mistep,nlist
                call nbondm()
                if (esymyes) call syminit()
                ncall = nlist
                call powel(tolf,ncall,estred,efcall)
                if (prll_on_off) call reduce_energies()
                tmp = 0.d0
                do 15 ig=1,npt_par
                        tmp=tmp+dpot(1,ig)*dpot(1,ig)
     1                          +dpot(2,ig)*dpot(2,ig)
     2                          +dpot(3,ig)*dpot(3,ig)
15                      continue
                if (prll_on_off) call reduce_1(tmp)
                tmp = dsqrt(tmp/npt3)
                if (tmp*tmp.lt.tolf) then
                 call eforce()
                 if (prll_on_off) call reduce_energies()
                  if (my_pe.eq.0) then
                   write(stdo,*)' structure number ',idyna
                   call wener(stdo)
                   write(umini,107)idyna,e_total,tmp
107                format(1x,i10,10x,f10.5,4x,f10.5)
                  end if
                 go to 17
                end if
16      continue
        call eforce()
        if (prll_on_off) call reduce_energies()
        if ( my_pe.eq.0) then
         write(umini,107)idyna,e_total,tmp
         call wener(stdo)
        end if
17      continue
        if (cstyl.eq.'DYNA') then
         if (my_pe.eq.0) 
     1  call wdyncrd(ucor_out,lpend-lpst+1,idyna,inofrz,nofreez,wbin)
        else if (cstyl.eq.'PATH') then
         if (my_pe.eq.0) call wpath(ucor_out)
        end if

18      continue

        write(stdo,108)
108     format(1x,/,1x,' *** Minimization completed ')

c***    HERE ENDS DYNAMICS MINIMIZATION
c***    AND MULTI-MINIMIZATION WITH VARIABLE CONSTRAINTS STARTS

        else
        
        do 19 i=1,nvar
         loop(i) = 1
19      continue

c 20 is the main loop on the constraints variables

        write(stdo,109)nvar,(var1(i)*pi180,i=1,nvar)
109     format(1x,//,' MINIMIZATION WITH ',i7,' CONSTRAINTS',//,
     1'Current equilibrium Positions',2(2x,f9.4),//,
     2'Below are constraints value and total energy')

20      continue
        
        do 21 i=1,nvar
         cnseq(i) = var1(i) + (loop(i)-1)*var2(i)
21      continue

c       do minimization
        
        efcall = 0

        do 23 istep=1,mistep,nlist
                write (*,*) 'mini step=',istep
                call nbondm()
                if (esymyes) call syminit()
                ncall = nlist
                call powel(tolf,ncall,estred,efcall)
                if (my_pe.eq.0) write(umini,102)istep+efcall-1
                tmp = 0.d0
                do 22 ig=1,npt_par
                 do 22 i=1,3
                 tmp=tmp+dpot(i,ig)*dpot(i,ig)
22              continue
                if (prll_on_off) call reduce_1(tmp)
                tmp = dsqrt(tmp/npt3)
                if (my_pe.eq.0) write(umini,103)tmp,dsqrt(tolf)
                if (tmp*tmp.lt.tolf) then
                 if (my_pe.eq.0) write(umini,104)
                 efcall = efcall - 1
                 call eforce()
                 if (prll_on_off) call reduce_energies()
                 if (my_pe.eq.0) call wener(umini)
                 go to 24
                end if
                efcall = efcall - 1
                call eforce()
                if (prll_on_off) call reduce_energies()
                if (my_pe.eq.0) call wener(umini)
23      continue
24      continue
        if (my_pe.eq.0) then
         call wpath(ucor_out)
         write(stdo,110)(cnseq(i)*pi180,i=1,nvar),e_total-e_cnst,tmp
110      format(1x,'tors> ',4(1x,f10.5))
        end if
        do 26 i=1,nvar
           loop(i) = loop(i) + 1
           if (loop(i).le.varstep(i)) then
              go to 27
           else
             if (i.lt.nvar) then
                 loop(i) = 1
              else
                 stop
              end if
           end if
26         continue
27         continue
cold    do 26 i=1,nvar-1
cold     if (loop(i).eq.varstep(i)) then
cold      do 25 j=1,i
cold       loop(j) = 1
cold25    continue
cold      loop(i+1) = loop(i+1) + 1
cold      go to 27
cold     else
cold      loop(i) = loop(i) + 1
cold      go to 26
cold     end if
cold26  continue
cold27  continue
cold    if (nvar.eq.1) then
cold     loop(1) = loop(1) + 1
cold    end if
cold    if (loop(nvar).gt.varstep(nvar)) stop
        go to 20
        end if
        stop
        end
