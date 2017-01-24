      Program SolvateCRD_gen
c====================================================================
c This program is a modification of SolvateCRD. The database file
c that uses may contain any type of molecule, not only water. The
c solute is inserted with the same non overlap rule as before
c====================================================================

ccccccc======================================================================
ccccccc
ccccccc     This program performs the solvatation of the protein molecule with
ccccccc     the water molecules stored in another database file. The output file
ccccccc     contains the initial protein molecule and water molecules that are
ccccccc     note overlap with the structure.
ccccccc
ccccccc     Version 1.00
ccccccc     Authors:  V. Zaloj, J. Meller
ccccccc
ccccccc=====================================================================

      implicit none
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'

c I/O input
      character*4 styl1,styl2
      character*4 name
      integer namel
      integer ucon,urcrd,uwcrd,ursbx,uwpol,ucons
c
      integer geti,of,ierr
      double precision getd
      logical find,fopen
c
      character*6 getchar
      integer level,lc
      character*5 ctyp
      character*4 mononame

c indexes
      integer i,j,n

c solvent/solute special arrays

      integer npt_solvent, npt_solute
      integer totmon_solvent,totmon_solute
      integer poimon_solvent(maxpt), poimon_solute(maxpt)
      integer poipt_solvent(maxmono),poipt_solute(maxmono)
      double precision coor_solvent(3,maxpt),coor_solute(3,maxpt)
      double precision ptms_solvent(maxpt),ptms_solute(maxpt)
      double precision epsgm12_solvent(maxpt),epsgm12_solute(maxpt)
      double precision epsgm6_solvent(maxpt),epsgm6_solute(maxpt)
      character*4 ptnm_solvent(maxpt),ptnm_solute(maxpt)
      character*4 moname_solvent(maxmono),moname_solute(maxmono)

c monomers skipping
      integer skip(maxmono)
      double precision sigma,cut,r2

c box and centering
      double precision xb,yb,zb,xb2,yb2,zb2
      double precision totmass
      double precision crdcent(3)

      character*1 junk

      data xb,yb,zb/3*1.0d0/

c adjusting concentration
      integer inc,inw,seed,idnc,idnw
      integer poi_water(maxpt),poi_cosolvent(maxpt)
      integer nions
      real random(1)
      double precision nc,nw,xc,xw,dnc,dnw
      character*4 cosolvent
      logical CORRECT

c
c
c     General initialization
c

      seed = 1
      CORRECT = .false.
      nions = 0
      lc     = 5
      stdi   = 5
      stdo   = 6
      stderr = 0
      totmon = 0
      npt    = 0
      nb     = 0
      nangl  = 0
      ntors  = 0
      nimp   = 0
      lestyp = 0
      nbulk  = 0
      debug  = .false.
      name   = 'solv'
      namel  = 4

      jnkf = 25
      open (unit=jnkf,status='scratch')

 1    continue

      call rline(name,namel,stdi)
      if (find('acti')) goto 2
      xb     = getd('xsbx',  xb  )
      yb     = getd('ysbx',  yb  )
      zb     = getd('zsbx',  zb  )
      if (find('CORR')) CORRECT=.TRUE.
      nions = geti('nion',nions)
c
      if (find('debu')) debug = .true.
      if (find('file')) then
	 if (find('wpol')) then
            uwpol = of()
	 end if
	 if (find('rsbx')) then
            ursbx = of()
            call getcrd(ursbx,'CHARM')
            do i=1,npt
             coor_solvent(1,i) = coor(1,i)
             coor_solvent(2,i) = coor(2,i)
             coor_solvent(3,i) = coor(3,i)
            enddo
	 end if
c connectivity file for the solute
	 if (find('conn')) then
            ucon = of()
            call rconn(ucon)
cccccc collect the information on the solute particle
            write (*,*) 'solute',npt,totmon
            npt_solute = npt
            totmon_solute = totmon
            do i=1,npt
             epsgm12_solute(i) = epsgm12(i)
             epsgm6_solute(i)  = epsgm6(i)
             ptms_solute(i)    = ptms(i)
             ptnm_solute(i)    = ptnm(i)
             poimon_solute(i)  = poimon(i)
             moname_solute(poimon(i)) = moname(poimon(i))
            enddo
            do i=1,totmon
             poipt_solute(i) = poipt(i)
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 end if
c connectivity file for the solvent
         if (find('cons')) then
            ucons = of()
            call rconn(ucons)
            write (*,*) 'solvent',npt,totmon
cccccc collect the information on the solvent particle
            npt_solvent = npt
            totmon_solvent = totmon
            do i=1,npt
             epsgm12_solvent(i) = epsgm12(i)
             epsgm6_solvent(i)  = epsgm6(i)
             ptms_solvent(i)    = ptms(i)
             ptnm_solvent(i)     = ptnm(i)
             poimon_solvent(i)  = poimon(i)
             moname_solvent(poimon(i)) = moname(poimon(i))
            enddo
            do i=1,totmon
             poipt_solvent(i) = poipt(i)
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         endif
	 if (find('rcrd')) then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','CHARM',lc)
            urcrd = of()
            call getcrd(urcrd,ctyp)
            do i=1,npt
             coor_solute(1,i) = coor(1,i)
             coor_solute(2,i) = coor(2,i)
             coor_solute(3,i) = coor(3,i)
            enddo
	 end if
	 if (find('wcrd'))then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','CHARM',lc)
            uwcrd = of()
	 end if
      end if

      goto 1

 2    continue

c     use getcrd for the solute
c      write (*,*) 'get coordinates for the solute'
c      npt = npt_solute
c      totmon = totmon_solute
c      do i=1,npt_solute
c       ptnm(i) = ptnm_solute(i)
c       poimon(i) = poimon_solute(i)
c       moname(poimon(i)) = moname_solute(poimon_solute(i))
c      enddo
c      do i=1,totmon_solute
c       poipt(i) = poipt_solute(i)
c      enddo
c      !call getcrd(urcrd,ctyp)

c place the center of mass of the solute at zero
      write (*,*) coor_solute(1,1)
      do  n=1,3
         crdcent(n)=0.d0
      enddo
      totmass = 0.d0
      do i=1,npt_solute
         totmass = totmass + ptms_solute(i)
         do n=1,3
            crdcent(n)= crdcent(n) + coor_solute(n,i) * ptms_solute(i)
         enddo
      enddo

      do n=1,3
         crdcent(n)= crdcent(n)/totmass
      enddo

      do i=1,npt_solute
         do n=1,3
            coor_solute(n,i) = coor_solute(n,i) - crdcent(n)
         enddo
      enddo

      write (*,*) coor_solute(1,1),totmass,npt_solute

c     now get the box sizes and prepare waters
      xb=abs(xb)
      yb=abs(yb)
      zb=abs(zb)
      xb2=xb/2.d0
      yb2=yb/2.d0
      zb2=zb/2.d0
c     the cutoff used to decide whether a certain molecule
c     will be kept in the solvated box is going to be 
c     2**(1/6)*sigma, which is the point at which the vdw
c     interaction is equal to zero.

      cut = 2.d0**(1.d0/6.d0)

c
      write(stdo,'(/10x,a/)') "Reading solvent box ..."

c      write (*,*) 'get coordinates for the solvent'
c      npt = npt_solvent
c      totmon = totmon_solvent
c      do i=1,npt_solvent
c       ptnm(i) = ptnm_solvent(i)
c       poimon(i) = poimon_solvent(i)
c       moname(poimon(i)) = moname_solvent(poimon_solvent(i))
c      enddo
c      do i=1,totmon_solvent
c       poipt(i) = poipt_solvent(i)
c      enddo
c      !call getcrd(ursbx,ctyp)

c      do i=1,npt_solvent
c       coor_solvent(1,i) = coor(1,i)
c       coor_solvent(2,i) = coor(2,i)
c       coor_solvent(3,i) = coor(3,i)
c      enddo

      do i=1,totmon_solvent
       skip(i) = 0
      enddo

      do i = 1,npt_solvent

         if (.not.((abs(coor_solvent(1,i)).lt.xb2).and.
     2             (abs(coor_solvent(2,i)).lt.yb2).and.
     3             (abs(coor_solvent(3,i)).lt.zb2))) then
          skip(poimon_solvent(i)) = 1
         endif

         if (skip(poimon_solvent(i)).eq.1) goto 100

         do j=1,npt_solute
            r2 = (coor_solvent(1,i) - coor_solute(1,j))**2 +
     2           (coor_solvent(2,i) - coor_solute(2,j))**2 +
     3           (coor_solvent(3,i) - coor_solute(3,j))**2
            sigma = (epsgm12_solute(j)*epsgm12_solvent(i))
            sigma = sigma/(epsgm6_solute(j)*epsgm6_solvent(i))
            sigma = sigma**(1.d0/6.d0)
            if (r2 .lt. (cut*sigma)**2) then 
             skip(poimon_solvent(i)) = 1
             goto 100
            endif
          enddo
100   continue   
      enddo

c copy all the coordinates of the solute in the new file
      do i=1,npt_solute
       coor(1,i) = coor_solute(1,i)
       coor(2,i) = coor_solute(2,i)
       coor(3,i) = coor_solute(3,i)
       ptnm(i)   = ptnm_solute(i)
       poimon(i) = poimon_solute(i)
       moname(poimon(i)) = moname_solute(poimon_solute(i))
      enddo
      npt = npt_solute
      totmon = totmon_solute


      if (CORRECT) then
c ONLY FOR A BINARY MIXTURE
c ensure that the ratio n(cosolvent)/(n(water)+n(cosolvent)) is correct
      cosolvent = 'XXXX'
      nw = 0.d0
      nc = 0.d0
c read without skipping
      do 111 i=1,totmon_solvent
       if (moname_solvent(i).eq.'TIP3'.or.
     1     moname_solvent(i).eq.'SPCE') then
        nw = nw + 1.d0
        inw = inw + 1
        poi_water(inw) = i
       else if (cosolvent.eq.'XXXX') then
        cosolvent = moname_solvent(i)
        nc = nc + 1.d0
        inc = inc + 1
        poi_cosolvent(inc) = i
       else if (moname_solvent(i).eq.cosolvent) then
        nc = nc + 1.d0
        inc = inc + 1
        poi_cosolvent(inc) = i
       else
        write (*,*) 'not a binary solvent...'
        write (*,*) 'the code will not fix the concentrations'
        write (*,*) 'of the different components'
        goto 110
       endif
111   continue
      xc = nc/(nc+nw)
      xw = nw/(nc+nw)
c read with skipping
      nw = 0.d0 - dble(nions)
      nc = 0.d0
      do 113 i=1,totmon_solvent
       if (skip(i).eq.1) goto 113
       if (moname_solvent(i).eq.cosolvent) then
        nc = nc + 1.d0
       else
        nw = nw + 1.d0
       endif
113   continue

      call RLUXGO(223,seed,0,0)  

c case 1 compute dnc
      idnc = anint((xc*(nw+nc)-nc)/(1.d0-xc))
      if (idnc.lt.0) then
       write (*,*) 'remove',idnc,' ',cosolvent

114   continue 
      if (idnc.ge.0 ) goto 112 
      call RANLUX(random,1)
      j = int(real(inc)*random(1)) 
      if (skip(poi_cosolvent(j)).eq.1) then
       goto 114
      else
       skip(poi_cosolvent(j))=1
       idnc = idnc + 1
      endif

      goto 114

      else

c case 2 compute dnw
      idnw = anint((xw*(nw+nc)-nw)/(1.d0-xw))
      if (idnw.lt.0) then
       write (*,*) 'remove',idnw,'waters'
      else
       write (*,*) 'dnw gt zero... error'
      endif    
115   continue 
      if (idnw.ge.0 ) goto 112
      call RANLUX(random,1)  
      j = int(real(inw)*random(1))            
      if (skip(poi_water(j)).eq.1) then
       goto 115
      else
       skip(poi_water(j))=1
       idnw = idnw + 1
      endif

      goto 115

      endif

112   continue

      nw = 0.d0 - dble(nions)
      nc = 0.d0
      do 116 i=1,totmon_solvent
       if (skip(i).eq.1) goto 116
       if (moname_solvent(i).eq.cosolvent) then
        nc = nc + 1.d0
       else
        nw = nw + 1.d0
       endif
116   continue

      write (*,*) 'xc expected',xc
      write (*,*) 'xc measured',nc/(nc+nw)

110   continue
      ENDIF



c copy the coordinates of the solvent only if not skipped
      do i=1,npt_solvent
       !write (*,*) poimon_solvent(i),skip(poimon_solvent(i))
       if (skip(poimon_solvent(i)).ne.1) then
        npt = npt + 1
         if (i.eq.1) then 
          totmon = totmon + 1
         else if (poimon_solvent(i).ne.poimon_solvent(i-1)) then
          totmon = totmon + 1
         endif
        coor(1,npt) = coor_solvent(1,i)
        coor(2,npt) = coor_solvent(2,i)
        coor(3,npt) = coor_solvent(3,i)
        ptnm(npt)   = ptnm_solvent(i)
        poimon(npt) = totmon 
        moname(totmon) = moname_solvent(poimon_solvent(i))
       endif
       enddo

      call putcrd(uwcrd,ctyp)

      close(uwcrd)
c
      mononame=BULK(1)
      write(uwpol,'(a/a,a4,a,i5.5/a)')
     $     "~","MOLC=(",mononame,")   #mon=",totmon,"~"
      write(uwpol,'(10a5)') (moname(i), i=1,totmon)
      write(uwpol,'(a4)') "*EOD"
      close(uwpol)
c
      write(stdo,'(/3(3x,a,f6.1)/)') ' xb=',xb,' yb=',yb,
     >     ' zb=',zb
      write(stdo,'(//3x,a,i7//3x,a,i7/)') " natom=",npt,
     >     " nmono=",totmon
      stop
      end
c--------------------------------------------------------






