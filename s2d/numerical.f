        program numer0
        
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
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
        include 'COMMON/SCNDRV.BLOCK'
        include 'COMMON/SC2.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'

        character*5 name
        double precision getd
        integer namel
        integer of,i
        integer geti
        logical find

        integer ucon,ucor,uwene
        integer ipick(maxpt)
        integer ucon1,ucon2

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

        jnkf = 25
        open (unit=jnkf,status='scratch')

        lcent = .false.
        ctrue = .true. 
        shift = .false.
        name  = 'numer'
        namel = 5
        ucon  = -1
        ucor  = -1
        debug = .false.


c energy default flags and parameters
        call init_ef()


1       continue

        call rline(name,namel,stdi)

        if (find('debu')) debug = .true.

        if (find('file')) then
         if (find('rcon').or.find('conn'))  then
          ucon  = of()
          call rconn(ucon)
          if (maxpt2d.lt.npt)  then
                write(*,*) ' maxpt2d is too small  '
                write(*,*) '-- check second derivative arrays '
                   call alert('numerical',9,'array size problem ',19,1)
                   stop
          end if
c initialze no freez vector
          inofrz = npt
          do 31 i=1,inofrz
                zerofrz(i) = 1
31        continue
         end if
         if (find('con1')) ucon1 = of()
         if (find('con2')) ucon2 = of()
         if (find('rcrd')) then
          ucor  = of()
          call getcrd(ucor,'CHARM')
         end if
         if (find('wene')) uwene = of()
        end if

        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        eps    = (getd('epsi',(eps)))
        if (.not. ctrue) ctrue  = find('cdie')
        if (find('rdie')) ctrue = .false.

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
c@ INSERT gbsa here
        if (find('acti')) go to 2
        go to 1
2       continue

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp



        call numer1(d2vec)

        stop
        end

        subroutine numer1(d2v)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
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
        include 'COMMON/SCNDRV.BLOCK'


        character*5 name
        integer namel,level
        integer i,ii,j,kx,ky,kz,jj

        double precision delta,factor
c vectors neededfor thedigonalization of the matrix
c
c       double precision eigenv(3*maxpt),work(3*maxpt)
        double precision d2v(3*npt,3*npt)


c rmax is maintained here for old input (with a single
c cutoff) to work
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
                cutmono2 = cutebig2*1.44
         else
                cutmono2 = cutmono2*cutmono2
         end if


        if (debug) then
                write(stdo,*)' after getcrd '
                do 21 i=1,npt
                 write(stdo,*)i,coor(1,i),coor(2,i),coor(3,i)
21              continue
        end if
        if(nmb.gt.0 .and.(D(nmb).le.0. .or. alpha(nmb).le.0.)) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
        call nbondm()
        if (esymyes) call syminit()
        do 20 i=1,3*npt
         do 20 j=1,3*npt
          d2v(i,j) = 0.d0
20      continue
        
        delta = 1.d-3
        factor = 1.d0/(2.d0*delta)
        do 22 i=1,npt
        coor(1,i) = coor(1,i) - delta
        call eforce()
        do 23 j=1,npt
                coor2(1,j) = dpot(1,j)
                coor2(2,j) = dpot(2,j)
                coor2(3,j) = dpot(3,j)
23      continue
        coor(1,i) = coor(1,i) + 2.d0*delta
        call eforce()
        coor(1,i) = coor(1,i) - delta
        kx = 3*(i-1)+1
        do 24 j=1,npt
                jj = 3*(j-1)+1
                d2v(jj,kx)   = factor*(dpot(1,j)-coor2(1,j))
                d2v(jj+1,kx) = factor*(dpot(2,j)-coor2(2,j))
                d2v(jj+2,kx) = factor*(dpot(3,j)-coor2(3,j))
24      continue
        coor(2,i) = coor(2,i) - delta
        call eforce()
        do 25 j=1,npt
                coor2(1,j) = dpot(1,j)
                coor2(2,j) = dpot(2,j)
                coor2(3,j) = dpot(3,j)
25      continue
        coor(2,i) = coor(2,i) + 2.d0*delta
        call eforce()
        coor(2,i) = coor(2,i) - delta
        ky = 3*(i-1)+2
        do 26 j=1,npt
                jj = 3*(j-1)+1
                d2v(jj,ky)   = factor*(dpot(1,j)-coor2(1,j))
                d2v(jj+1,ky) = factor*(dpot(2,j)-coor2(2,j))
                d2v(jj+2,ky) = factor*(dpot(3,j)-coor2(3,j))
26      continue
        coor(3,i) = coor(3,i) - delta
        call eforce()
        do 27 j=1,npt
                coor2(1,j) = dpot(1,j)
                coor2(2,j) = dpot(2,j)
                coor2(3,j) = dpot(3,j)
27      continue
        coor(3,i) = coor(3,i) + 2.d0*delta
        call eforce()
        coor(3,i) = coor(3,i) - delta
        kz = i*3
        do 28 j=1,npt
                jj = 3*(j-1)+1
                d2v(jj,kz)   = factor*(dpot(1,j)-coor2(1,j))
                d2v(jj+1,kz) = factor*(dpot(2,j)-coor2(2,j))
                d2v(jj+2,kz) = factor*(dpot(3,j)-coor2(3,j))
28      continue
22      continue
        write(stdo,*)
        write(stdo,*)' Second derivative matrix '
c       do 29 i=1,3*npt
c        write(stdo,100)(d2v(i,j),j=1,3*npt)
c100     format(12(f5.1,1x))
c29     continue
        do 291 i=1,npt
         do 291 j=1,npt
          ii = 3*(i-1) + 1
          jj = 3*(j-1) + 1
          write(stdo,*)i,j
          write(*,100)d2v(ii,jj),d2v(ii+1,jj),d2v(ii+2,jj)
          write(*,100)d2v(ii,jj+1),d2v(ii+1,jj+1),d2v(ii+2,jj+1)
          write(*,100)d2v(ii,jj+2),d2v(ii+1,jj+2),d2v(ii+2,jj+2)
100       format(3(f9.2,1x))
291     continue

c       call house(d2v,3*npt,3*npt,eigenv,work,i)
c       write(stdo,*)' Eigenvalues '
c       write(stdo,100)(eigenv(j),j=1,3*npt)
c       write(stdo,*)' Error = ',i

        return
        end
