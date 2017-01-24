       subroutine line_loop(uwcrd,dtopt,npri,
     &           nstep,temper,constant_tmp,clo,interpolate,lap,
     &           select,nlist,ucon,urcrd,amid_true)
       implicit none
c
c calculate path using SDEL algorithm
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PATH.BLOCK'
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
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        double precision dtopt,temper,clo
        integer uwcrd,nstep,npri,interpolate
        integer ucon,urcrd,nlist,of
        logical constant_tmp,lap,select,amid_true

C local variables        

        character*10 name
        integer namel,geti
        double precision getd
        logical find

        name   = 'line_loop'
        namel  = 10
c ....initialization of the parallel machine
       
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
           if (find('conn')) then
               ucon = of()
               call rconn(ucon)
           endif

           if (find('rcrd').and.my_pe.eq.0) urcrd = of()
           if (find('wcrd').and.my_pe.eq.0) uwcrd = of()

C   read coarse grained model parameters (associated files)
           call Read_CG_input()

        end if

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

        dtopt=getd('dtop',dtopt)
        igrid   = geti('grid',igrid)
        skpno   = geti('skno',skpno)

        if (find ('noRA')) Random_velocities = .false.

        if (find ('gbsa')) gbsabool=.true.
        gbsu=geti('gbsu',gbsu)
        nstep   = geti('#ste',nstep)
        npri    = geti('#pri',npri)
        nwcrd   = geti('#wcr',nwcrd)
        temper  = getd('tmpr',temper)
        if (find ('ctmp')) constant_tmp = .true.

        numprocs    = geti('proc',numprocs)
        gamma   = (getd('gama',(gamma)))
        clo = getd('clog',clo)
        interpolate=(geti('itpl',(interpolate)))

        if (find('ovlp')) lap = .true.

        if (find('sele')) select = .true.
c Energy parameters
        if (find('hvdw')) hvdw0 = .false.
        nlist   = geti('list',nlist)
        cutvdw2  = (getd('rvmx',(cutvdw2)))
        cutvbig2 = (getd('rvbg',cutvbig2))
        cutele2  = (getd('relx',(cutele2)))
        cutebig2 = (getd('rebg',cutebig2))
        cutmono2 = getd('cutm',cutmono2)
        rmax     = getd('rmax',rmax)
        eps     = (getd('epsi',(eps)))

        if (find('cdie')) ctrue = .true.
        if (find('rdie')) ctrue = .false.


        if (find('cpth')) crdstyl = 'PATH'
        if (find('cchr')) crdstyl = 'CHAR'
        if (find('cini')) crdstyl = 'INIT'
        if (find('cint')) crdstyl = 'INTR'

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

        Hamilt  = getd('hami',Hamilt)
        
        irand   = geti('rand',irand)

        if (find('cpth')) crdstyl = 'PATH'
        if (find('cchr')) crdstyl = 'CHAR'
        if (find('cini')) crdstyl = 'INIT'
        if (find('cint')) crdstyl = 'INTR'
        if (find('acti')) go to 2       

        go to 1
2       continue

        return

        end
