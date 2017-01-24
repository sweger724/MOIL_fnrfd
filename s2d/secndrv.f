      program secndrv
       
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
      include 'COMMON/GETVEC.BLOCK'
c      include 'COMMON/SC2.BLOCK'

      character*7 name
      double precision getd
      integer namel,level
      integer of,i,j,ii,jj
      logical find,numer,analyt,compar,eigen
      
      integer ucon,ucor,uwene
      integer ipick(maxpt)
      integer ucon1,ucon2

      double precision d2vn(3*maxpt,3*maxpt),d2va(3*maxpt,3*maxpt)
      double precision diff, toler
            
      name  = 'secndrv'
      namel = 7
      ucon  = -1
      ucor  = -1
      toler = 1.d-6
      numer=.false.
      analyt=.false.
      compar=.false. 
      eigen=.false.

c default flags and parameters
      call init_d2v()
c
c energy  flags and parameters
      call init_ef()

c.....Parameter input from stdi using the line interpreter routines.
c.....Each process reads its own input.
c.....loop of the kind "do until"
c
 1    continue

         call rline(name,namel,stdi)

         if (find('debu')) debug = .true.

         if (find('file')) then
            if (find('rcon').or.find('conn'))  then
               ucon  = of()
               call rconn(ucon)
c              initialize no freez vector
               inofrz = npt
               do 31 i=1,inofrz
                  zerofrz(i) = 1
 31            continue
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
	 cutmono2 = getd('cutm'cutmono2)
         rmax     = getd('rmax',rmax)
         eps      = (getd('epsi',(eps)))
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
 35         continue
         endif
         
         if (find('anlt')) analyt=.true.
         if (find('numr')) numer =.true.
         if (find('cmpr')) compar=.true.


         if (find('acti')) go to 2

      go to 1
 2    continue
c
c
c rmax is maintained here for old input (with a single cutoff) to work
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

      if (nmb.gt.0 .and.(D(nmb).le.0. .or. alpha(nmb).le.0.)) then
         level=1
         call alert(name,namel,'Dmor or alph is 0',16,level)
      end if
c
      if (debug) then
         write(stdo,*)' after getcrd '
         do 21 i=1,npt
            write(stdo,*)i,coor(1,i),coor(2,i),coor(3,i)
 21      continue
      end if
c
c initialize the second derivative matrix
      do 110 i=1,npt*3
         do 111 j=1,npt*3
            d2vn(i,j)=0.d0
            d2va(i,j)=0.d0
 111     continue
 110  continue
c 
      if (analyt) then
         call d2vanlyt(d2va,eigen)
c
         if (debug)  then
            write(stdo,*)
            write(stdo,*)' Second derivative matrix (analytic)'
c            do 29 i=1,3*npt
c               write(stdo,100)(d2v(i,j),j=1,3*npt)
c100            format(12(f5.1,1x))
c29          continue
            do 291 i=1,npt
               do 292 j=1,npt
                  ii = 3*(i-1) + 1
                  jj = 3*(j-1) + 1
                  write(stdo,*)i,j
                  write(*,101)d2va(ii,jj),d2va(ii+1,jj),d2va(ii+2,jj)
                  write(*,102)d2va(ii,jj+1),d2va(ii+1,jj+1),
     1                 d2va(ii+2,jj+1)
                  write(*,103)d2va(ii,jj+2),d2va(ii+1,jj+2),
     1                 d2va(ii+2,jj+2)
 101              format(1x,'xx,xy,xz',3(f9.2,1x))
 102              format(1x,'yx,yy,yz',3(f9.2,1x))
 103              format(1x,'zx,zy,zz',3(f9.2,1x))
 292           continue
 291        continue
         end if
c
      end if
c
      if (numer) then
         call d2vnumer(d2vn)
c
         if (debug)  then
            write(stdo,*)
            write(stdo,*)' Second derivative matrix (numeric)'
            do 391 i=1,npt
               do 392 j=1,npt
                  ii = 3*(i-1) + 1
                  jj = 3*(j-1) + 1
                  write(stdo,*)i,j
                  write(*,101)d2vn(ii,jj),d2vn(ii+1,jj),d2vn(ii+2,jj)
                  write(*,102)d2vn(ii,jj+1),d2vn(ii+1,jj+1),
     1                 d2vn(ii+2,jj+1)
                  write(*,103)d2vn(ii,jj+2),d2vn(ii+1,jj+2),
     1                 d2vn(ii+2,jj+2)
 392           continue
 391        continue
         end if
c
      end if
c 
      if (compar) then
         do 120 i=1,npt*3
            do 121 j=1,npt*3
               diff=d2va(i,j)-d2vn(i,j)
               if (abs(diff).gt.toler) then
                  ii=mod(i-1,3)+1
                  jj=mod(j-1,3)+1
                  write(stdo,*)'  large discrepancy in the derivative',
     &                 '  with respect to particles'
                  write(stdo,*) (i-1)/3+1,
     &              '(',ii,') and',(j-1)/3+1, '(',jj,'). diff=', diff 
               end if
 121        continue
 120     continue
      end if
c
c
      end
