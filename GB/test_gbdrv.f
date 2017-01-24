c test derivative, copied from energy.f 
	program test_gbdrv 
	
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
	include 'COMMON/SGB2DRV.BLOCK'
	include 'COMMON/EBALL.BLOCK'
	character*4 name,coortyp
	double precision getd
	integer namel,level,nstru
	integer of,geti
	integer i,n,nstrub
        integer rbin
	logical find

	integer ucon,ucor,uwene
	integer ipick(maxpt)
	integer ucon1,ucon2

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
	surften=0.005d0
	call init_ef()



1	continue

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
31	  continue
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
56	 continue
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

	if (find('acti')) go to 2
	go to 1
2	continue

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
21		continue
	end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp

	if (gbsabool) call make_rborn
	if(nmb.gt.0) then
	 if (D(nmb).le.0.d0 .or. alpha(nmb).le.0.d0) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
	end if
	 if (my_pe.eq.0) call init_wre(uwene)
	 if (coortyp.eq.'CHAR') then
	  call nbondm()
	  if (esymyes) call syminit()
c =====================================================
        if(specl) call nbondm_spcl()
c ====================================================

C@
C	write(*,*)' before calling eforce '
	  call eforce()
	    if (my_pe.eq.0) call wener(uwene)
	 else if (coortyp.eq.'PATH') then

C AVI MODIFICATION :: fix rpath reads so as not to rewind
	   rewind(ucor)
C END MODIFICATION
	   do 30 i=1,nstru

C CHANGED BELOW from rpath(ucor,i)
	    call rpath_seq(ucor,1)
	    call nbondm()
	    if (esymyes) call syminit()
		call eforce()
	    if (my_pe.eq.0) call wener(uwene)
30	   continue
	 else if (coortyp.eq.'DYNA') then
	   rewind ucor
	   norew =  .true.
	   do 32 i=1,nstru
	    call rdyncrd(ucor,i,inofrz,nofreez,rbin)
	    call nbondm()
	    if (esymyes) call syminit()
	    call eforce()
	    if (my_pe.eq.0) call wener(uwene)
32	   continue
	  else
	   level = 1
	   call alert(name,namel,'Unknown coordinate type',23,level)
	  end if
c ---------------------------------------------------------------------  
          
         call testderi2nd(npt)

	stop
	end

c ---------------------------------------------------------------------
c ** the following subroutines are included for debugging purpose **
c    better to be  deleted from production code
c --------------------------------------------------------------------- 

c --------------------------------------------------------------------- 
      SUBROUTINE egb_all_pair(natom,qmoil,diel_ext,
     $     ifsndbool,moilvector,moilvectorout)
c ---------------------------------------------------------------------

      IMPLICIT NONE

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'

c - - - arg types - - -
      integer natom,ifsndbool,moildebug
      double precision qmoil(*),moilvector(*),moilvectorout(*),
     $     diel_ext, multfact

c Local variables
      double precision q(MAXPT), epol
      integer i, j

c Begin execution

      do i = 1, natom
          sumdeijda(i) = 0.0
	  dpot(1,i) = 0.0
	  dpot(2,i) = 0.0
	  dpot(3,i) = 0.0
      end do

      multfact = 18.220867158288597763891252209377D0
      do i = 1,natom
         q(i) = qmoil(i)*multfact
      end do 

      epol = 0.0
      egb_inter_ene = 0.0

      do i = 1,natom
         do j= i+1, natom
            if ( (q(i) .ne.0D0) .and. ( q(j) .ne. 0D0)) THEN
                if (lesid(i).eq.0 .or. lesid(j).eq.0) then
                   call egb_calc_pair_f(coor(1,i), coor(1,j),
     $                 q(i), q(j), reff(i), reff(j), epol,
     $                 dpot(1,i),dpot(1,j),
     $                 moilvectorout(3*(i-1)+1),
     $                 moilvectorout(3*(j-1)+1), 
     $                 moilvector(3*(i-1)+1),
     $                 moilvector(3*(j-1)+1),
     $                 massfac(i), massfac(j),
     $                 diel_ext, ifsndbool,i,j)
                end if
            end if 
         end do
      end do 

      egb_inter_ene = epol
c      write(*,*) "E_pol = ", egb_inter_ene 
      return
      END

c --------------------------------------------------------------------- 
      SUBROUTINE egb_radii(natom)
c ---------------------------------------------------------------------

      IMPLICIT NONE

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

c Arguments
      integer natom

c     - - - local declarations - - -

      integer count,i,j,k

      double precision 
     &        xi,yi,zi,ri,ri1i,sumi,xij,yij,zij,r2,dij1i,dij,sj,sj2,
     &        temp,uij,dij2i,tmpsd,dumbo,theta


c YS variables for Onufriev's GB model 
      double precision BOFFSET, KSCALE
      parameter ( BOFFSET = 0.09D0, KSCALE = 0.73D0 )

c YS LCPO stuff follows
      double precision rgbmax2,rgbmax1i,rgbmax2i,tmpcs, rgbmaxpsmax2 

c YS FGB taylor coefficients follow from A to H;
c    1/3, 2/5, 3/7, 4/9, 5/11, 4/3, 12/5, 24/7, 40/9, 60/11
      double precision TA, TB, TC, TD, TDD, TE, TF, TG, TH, THH
      parameter (TA =0.33333333333333333333)
      parameter (TB =0.4D0)
      parameter (TC =0.42857142857142857143)
      parameter (TD =0.44444444444444444444)
      parameter (TDD=0.45454545454545454545)
      parameter (TE =1.33333333333333333333)
      parameter (TF =2.4D0)
      parameter (TG =3.42857142857142857143)
      parameter (TH =4.44444444444444444444)
      parameter (THH=5.45454545454545454545)

c Begin

c -----------------------------------------------------------
c Algorithm Step 1
c Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275
c -----------------------------------------------------------

c YS  Get the "effective" Born radii via the approximate pairwise
c     method Use Eqs 9-11 of Hawkins, Cramer, Truhlar; J. Phys. 
c     Chem. 100:19824(1996). Also added a modified method based on
c     Onufriev, Bashford, Case, Proteins, 55:383(2004). eqs 6
c
      rgbmax2 = rgbmax*rgbmax
      rgbmax1i = 1.0/rgbmax
      rgbmax2i = rgbmax1i*rgbmax1i
      rgbmaxpsmax2 = ( rgbmax + Fsmax ) * ( rgbmax + Fsmax )

      count = 1      
      do i = 1,natom
          xi = coor(1, i)
          yi = coor(2, i)
          zi = coor(3, i)
          ri = rborn( i ) - BOFFSET 
          ri1i = 1.0D0/ri
          sumi = 0.0D0 

          do j = 1,natom
              if ( i .eq. j ) go to 100 
              if (lesid(i).eq.0 .or. lesid(j).eq.0) then

                  xij = xi - coor(1, j)
                  yij = yi - coor(2, j)
                  zij = zi - coor(3, j)

                  r2 = xij*xij + yij*yij + zij*zij
                  if ( r2 .gt. rgbmaxpsmax2 ) go to 100
                  
                  dij1i = 1.0D0/DSQRT(r2)
                  dij = r2*dij1i
                  sj = fs( j ) * ( rborn( j ) - BOFFSET  )
                  sj2 = sj * sj
 
c YS: from NAB code: 
c     followings are from the Appendix of Schaefer and Froemmel
c     JMB 216:1045-1066, 1990; Taylor series expansion for d >> s
c     is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
c     Andreas Svrcek-Seiler and Alexy Onufriev  

                  if ( dij .gt. rgbmax+sj ) go to 100
                  if ( dij .gt. rgbmax-sj ) then 
                     uij = 1.0 / ( dij - sj )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uij+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - sj2)+
     &                      2.0 * DLOG((dij - sj) * rgbmax1i))
                     sumi = sumi - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = sj2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sumi = sumi - sj * tmpsd * dij2i * dumbo                
                  else if (dij .gt. ri + sj ) then 
		     temp =  0.5 *( sj/(r2-sj2) +
     $                    0.5 * dij1i * DLOG((dij-sj)/(dij+sj)))
                     sumi = sumi - temp
                  else if (dij .gt. DABS(ri-sj) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + ri*ri -sj2)
                     uij = 1.0 / (dij+sj)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uij + dij1i*DLOG(ri*uij))
                     sumi = sumi - temp
                  else if ( ri .lt. sj ) then
		     temp = 0.5 * (sj / (r2 - sj2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((sj - dij)/(sj + dij)))
                     sumi = sumi - temp

                  end if 
              end if

              if( gbnpbool ) then 
                  if( (P0(ptsaid(i)) + P0(ptsaid(j))) .gt. dij ) then
                      if ( (P0(ptsaid(i)) .gt. 2.5 ) .and. 
     &                     (P0(ptsaid(j)) .gt. 2.5 ) ) then
                           ineighbor(count) = j
                           count = count + 1
                      end if
                  end if
              end if 

100       continue  
          end do 

          if(gbobcbool) then 
C YS  Onufriev effective radii
C     note sumi is negative here
               psi(i) = -ri*sumi
               reff(i) = 1.0 / (ri1i-DTANH(( gbalpha - gbbeta*psi(i)+
     &                  gbgamma*psi(i)*psi(i))*psi(i))/rborn(i))

          else 
C YS "standard" HCT effective radii
               reff( i ) = 1.0 /(ri1i+ sumi)
          end if 
          if ( gbnpbool ) then 
                ineighbor(count) = 0
                count = count + 1
          end if
c print out born radii for debug
c          write(*,*) "BORN RADII> i=", i, " reff=", reff(i)

      end do 
      return
      END

c --------------------------------------------------------------------- 
      SUBROUTINE egb_self(natom,qmoil,diel_ext)

c YS     Computing  the "diagonal" term
c        Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

c Parameters
      integer natom
      double precision qmoil(*), diel_ext

c Local variables
      integer i
      double precision KSCALE, BOFFSET
      parameter (KSCALE=0.73D0, BOFFSET=0.09D0)
      double precision  dax(MAXPT), day(MAXPT), daz(MAXPT)      
      double precision rgbmaxpsmax2, r2, dij1i, dij2i,dij3i, dij,
     &                 temp1, temp7, sj,sj2,datmp, rgbmax2, rgbmax1i,
     &                 rgbmax2i, dumbo, tmpsd
      double precision sumda, thi, multfac
      double precision q(MAXPT)


      
c Begin

      rgbmax2 = rgbmax*rgbmax
      rgbmax1i = 1.0 / rgbmax
      rgbmax2i = rgbmax1i*rgbmax1i
      rgbmaxpsmax2 = ( rgbmax + Fsmax ) * ( rgbmax + Fsmax )
    
      multfact = 18.220867158288597763891252209377D0
      do i = 1,natom
         q(i) = qmoil(i)*multfact
      end do 
c ---------------------------------------------------------------------
c Compute the "diagnal" energy that is a function of only the
c effective radius Ri. Also compute the contribution of 
c the diagnal energy term to the sum by which the derivative of
c Ri will be multiplied
c --------------------------------------------------------------------- 

      
      egb_self_ene = 0.0
      do i =1, natom 
	 expmkf = DEXP(-KSCALE * kappa * reff(i) )/ diel_ext
	 dielfac = 1.0D0 - expmkf
	 qi2h = 0.5 * q(i)*q(i)
	 qid2h = qi2h * dielfac
         egb_self_ene = egb_self_ene - qid2h / reff(i)
c	 sumdeijda(i) = sumdeijda(i) + qid2h - KSCALE*kappa*qi2h*
c     &               expmkf*reff(i)
      end do
c      write(*,*) "EGB> Self-energy: ", egb_self_ene

      return
      END


c ---------------------------------------------------------------------
      SUBROUTINE deriv_rc(i,j,datmpx,datmpy,datmpz)
c i,j - atom i and j
c datmpx - derivative of radius of atom i w.r.t to x coord of atom j
      
      implicit none
c args
      integer i, j
      double precision datmpx, datmpy, datmpz

      return
      END


c ---------------------------------------------------------------------
      SUBROUTINE check_derivative(natom)

      IMPLICIT none
      INCLUDE 'COMMON/LENGTH.BLOCK'
      INCLUDE 'COMMON/SGB.BLOCK'
      INCLUDE 'COMMON/CONNECT.BLOCK'
      INCLUDE 'COMMON/COORD.BLOCK'
      INCLUDE 'COMMON/ENERGY.BLOCK'

c args type
      INTEGER natom, i, j

c local variables
      DOUBLE PRECISION deltad,temp,deri, U1,U2, diff,sumdiff
      PARAMETER (deltad= 0.0000001)
      DOUBLE PRECISION analydx, analydy, analydz

c begin
c      write (*,*) "natom=", natom

c ---------------------------------------------------------------------
c Checking derivatives of total energies w.r.t Cartesian coordinates
c --------------------------------------------------------------------- 
      write(*,*) " "
      write(*,*) "CHK_DERI> =======================================" 
      write(*,*) "CHK_DERI>  Checking derivatives of total "
      write(*,*) "CHK_DERI>  energies w.r.t coords"
      write(*,*) "CHK_DERI> =======================================" 
      sumdiff = 0
      do i = 1, natom 
          write(*,*) "CHK_DERI> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"  
          write(*,*) "CHK_DERI> Atom:", i  
          write(*,*) "CHK_DERI> ** obtaining analytical derivatives"
          call eforce()
          analydx = dpot(1,i)
          analydy = dpot(2,i)
          analydz = dpot(3,i)


c          write(*,*) "CHK_DERI> *** compare dE/dX"  
c          write(*,*) "CHK_DERI>  calculate numerical dE/dX"
          temp = coor(1,i)
          coor(1,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(1,i) = temp - deltad
          call eforce()
          U2 = e_total
          deri = (U1-U2)/(2*deltad)

          coor(1,i) = temp
	  diff = DABS(deri - analydx)

          write(*,*) "CHK_DERI> X  Num Anal: ", deri,analydx  
          write(*,*) "CHK_DERI> atom diff: ",i, diff  
	  
          sumdiff = sumdiff + diff

c          write(*,*) "CHK_DERI> *** compare dE/dY"  
c          write(*,*) "CHK_DERI>  calculate numerical dE/dY"
          temp = coor(2,i)
          coor(2,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(2,i) = temp - deltad
          call eforce()
          U2 = e_total
          deri = (U1-U2)/(2*deltad)

          coor(2,i) = temp
          
	  diff = DABS(deri - analydy)

          write(*,*) "CHK_DERI> Y  Num Anal: ", deri,analydy  
          write(*,*) "CHK_DERI> atom diff: ",i, diff  
          sumdiff = sumdiff + diff

c          write(*,*) "CHK_DERI> *** compare dE/dZ"  
      
c calculate numerical dE/dZ
          temp = coor(3,i)
          coor(3,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(3,i) = temp - deltad
          call eforce()
          U2 = e_total
          deri = (U1-U2)/(2*deltad)

          coor(3,i) = temp
          
	  diff = DABS(deri - analydz)

          write(*,*) "CHK_DERI> Z  Num Anal: ", deri,analydz  
          write(*,*) "CHK_DERI> atom diff: ",i, diff  
          sumdiff = sumdiff + diff
      end do

      write(*,*) "CHK_DERI> ***********************************" 
      write(*,*) "CHK_DERI> Average diff : ", sumdiff/(3*natom) 

c skip other checkings for derivatives
      go to 9001 

c ---------------------------------------------------------------------
c Checking derivatives of GB polarized solvation energy
c w.r.t effective radii, assuming coordiates  hold constant
c --------------------------------------------------------------------- 
      sumdiff = 0.0

      write(*,*) " "
      write(*,*) "DERI_RADI> =======================================" 
      write(*,*) "DERI_RADI>  Checking derivatives of GB interaction "
      write(*,*) "DERI_RADI>  energy w.r.t effective radii "
      write(*,*) "DERI_RADI> =======================================" 
      
      call egb_radii(npt)
      call egb_all_pair(npt, ptchg, dielectric, 0, 0.d0,0.d0)

      do i = 1, natom
          analydx = DABS(sumdeijda(i) / ( reff(i) * reff(i) ))
          write(*,*) "DERI_RADI> ******* i = " , i
          write(*,*) "DERI_RADI> Analytical dE/dRi = ", analydx  
          write(*,*) "DERI_RADI> Now computing numerical ... " 
	  temp = reff(i)
	  reff(i) = temp + deltad
          call egb_all_pair(npt, ptchg, dielectric, 0, 0.d0,0.d0)
	  U1 = egb_inter_ene
	  reff(i) = temp - deltad
          call egb_all_pair(npt, ptchg, dielectric, 0, 0.d0,0.d0)
	  U2 = egb_inter_ene
	  reff(i) = temp
          deri = DABS((U1-U2)/(2*deltad))
          if ( deri .ne. 0 ) then
	      diff = DABS(deri - analydx)/deri
          else
	      diff = 0
          end if
          write(*,*) "DERI_RADI> Numerical  dE/dRi: ", deri  
          write(*,*) "DERI_RADI> Analytical dE/dRi: ", analydx 
          write(*,*) "DERI_RADI> diff (%)         : ", diff 
          sumdiff = sumdiff + diff
      end do
      write(*,*) "DERI_RADI> ***********************************" 
      write(*,*) "DERI_RADI> Average diff : ", sumdiff/natom 


c ---------------------------------------------------------------------
c Checking derivatives of GB interaction energy w.r.t Cartesian coord
c assuming effective radii hold constant
c --------------------------------------------------------------------- 

      write(*,*) " "
      write(*,*) "DERI_INTE> =======================================" 
      write(*,*) "DERI_INTE>  Checking derivatives of GB interaction "
      write(*,*) "DERI_INTE>  energies w.r.t coords"
      write(*,*) "DERI_INTE> =======================================" 

      sumdiff = 0.0
      do i = 1, natom 

          write(*,*) "DERI_INTE> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"  
          write(*,*) "DERI_INTE> Atom:", i  

          write(*,*) "DERI_INTE> ** obtaining analytical derivatives"
	  call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
	  analydx = DABS(dpot(1,i))
	  analydy = DABS(dpot(2,i))
	  analydz = DABS(dpot(3,i))

          write(*,*) "DERI_INTE> *** compare dE/dX"  
          write(*,*) "DERI_INTE>  calculate numerical dE/dX"
          temp = coor(1,i)
          coor(1,i) = temp + deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U1 = egb_inter_ene
          coor(1,i) = temp - deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U2 = egb_inter_ene
          deri = DABS((U1-U2)/(2*deltad))
          coor(1,i) = temp
          if ( deri .ne. 0 ) then
		  diff = DABS(deri - analydx)/deri
          else
                  diff = 0
          end if
          write(*,*) "DERI_INIE> Numerical  dE/dX: ", deri  
          write(*,*) "DERI_INIE> Analytical dE/dX: ", analydx 
          write(*,*) "DERI_INIE> diff %          : ", diff 
          sumdiff = sumdiff + diff

          write(*,*) "DERI_INTE> *** compare dE/dY"  
          write(*,*) "DERI_INTE>  calculate numerical dE/dY"
          temp = coor(2,i)
          coor(2,i) = temp + deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U1 = egb_inter_ene
          coor(2,i) = temp - deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U2 = egb_inter_ene
          deri = DABS((U1-U2)/(2*deltad))
          coor(2,i) = temp
          if ( deri .ne. 0 ) then
	      diff = DABS(deri - analydy)/deri
          else
	      diff = 0
          end if
          write(*,*) "DERI_INIE> Numerical  dE/dY: ", deri  
          write(*,*) "DERI_INIE> Analytical dE/dY: ", analydy 
          write(*,*) "DERI_INIE> diff %          : ", diff 
          sumdiff = sumdiff + diff

          write(*,*) "DERI_INTE> *** compare dE/dZ"  
          write(*,*) "DERI_INTE>  calculate numerical dE/dZ"
          temp = coor(3,i)
          coor(3,i) = temp + deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U1 = egb_inter_ene
          coor(3,i) = temp - deltad
          call egb_all_pair(npt,ptchg,dielectric,0,0.d0,0.d0)
          U2 = egb_inter_ene
          deri = DABS((U1-U2)/(2*deltad))
          coor(3,i) = temp
          if ( deri .ne. 0 ) then
	      diff = DABS(deri - analydz)/deri
          else
	      diff = 0
          end if
          write(*,*) "DERI_INIE> Numerical  dE/dZ: ", deri  
          write(*,*) "DERI_INIE> Analytical dE/dZ: ", analydz 
          write(*,*) "DERI_INIE> diff %          : ", diff 
          sumdiff = sumdiff + diff
      end do

      write(*,*) "DERI_INTE> ***********************************" 
      write(*,*) "DERI_INTE> Average diff : ", sumdiff/(3*natom) 

c ---------------------------------------------------------------------
c Checking derivatives of effective Born radius of atom i w.r.t 
c Cartesian coord of atom j
c --------------------------------------------------------------------- 

      write(*,*) " "
      write(*,*) "DERI_RADI2> =======================================" 
      write(*,*) "DERI_RADI2>  Checking derivatives of radii  "
      write(*,*) "DERI_RADI2>  w.r.t coords (not done)"
      write(*,*) "DERI_READ2> =======================================" 

c      do i = 1, natom 
c          do j = 1, natom
c	      write(*,*) "DERI_RADI2> Radius of atom i, 
c     & coord of atom j", i,j
	     
c	  end do 
c      end do

9001  return
      end
