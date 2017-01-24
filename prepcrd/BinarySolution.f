      Program BinarySolution

C this program is based upon SolvateCRD
C it creates a binary solution of water + other solute
C it needs:
C A) a crd and wcon of a single molecule of the solute
C B) the usual MOIL file with the water box
C C) the relative number of particles of solute in the mixture, ie xs = ns/(ns + nw)
C D) the volume of the box
C it works as following
C call NRemWat = n water removed to add a solute/ns (average number of water removed)
C so
C xs = ns/(ns + nw(0) - NRemWat*ns)
C START LOOP
C 1) introduce a solute
C 2) remove waters and update NRemWat
C 3) is xs within boundaries from expected? 
C    YES -> GOTO END LOOP
C    NO -> is it smaller? YES -> GOTO START LOOP
C                         NO  -> some problem occurred: end the program
C END LOOP

c======================================================================
c
c     This program performs the solvatation of the protein molecule with
c     the water molecules stored in another database file. The output file
c     contains the initial protein molecule and water molecules that are
c     note overlap with the structure.
c
c     Version 1.00
c     Authors:  V. Zaloj, J. Meller
c
c=====================================================================
      implicit none
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/MUTA.BLOCK'
c
      character*14 name
      integer namel
      integer ucon,urcrd,uwcrd,urwbx,uwpol
      character*1 junk
c
      integer geti,of
      double precision getd
      logical find,fopen
c
      integer level

      integer i,j,k,l,loop
      double precision e,ovdw,hvdw,omass,hmass

      integer nrec,NW0,remove(maxpt),id(maxpt)
      double precision water(3,maxpt),iWatRemoved(maxpt)
      double precision solute(3,maxpt),reference(3,maxpt)
      double precision O(3,maxpt),H1(3,maxpt),H2(3,maxpt)
      double precision ptvdw(maxpt)

      integer nptS
      character*4 ptnmS(maxpt),SoluteName

      integer NS
      double precision XS,XSt,NRemWat,tol

      double precision xb,yb,zb
      double precision cm(3)
      double precision r,dx,dy,dz
      double precision xx,yy,zz,xx1,yy1,zz1,xx2,yy2,zz2

      integer kinit,kend
      double precision totmass,massd,ratio

      integer irand
      real random(3)

c
c
c     General initialization
c
      stdi   = 5
      stdo   = 6
      stderr = 0
      totmon = 0
      npt    = 0
      name   = 'BinarySolution'
      namel  = 14
      XS = 0.d0
      xb = 99.d0
      yb = 99.d0
      zb = 99.d0

      massd = -1.d0

      loop = 0

      jnkf = 25
      open (unit=jnkf,status='scratch')

 1    continue

      call rline(name,namel,stdi)
      if (find('acti')) goto 2
      xb     = getd('xwbx',  xb  )
      yb     = getd('ywbx',  yb  )
      zb     = getd('zwbx',  zb  )
c
      if (find('debu')) debug = .true.
      if (find('file')) then
	 if (find('wpol')) then
            uwpol = of()
	 end if
	 if (find('rwbx')) then
            urwbx = of()
	 end if
	 if (find('wcon') .or. find('conn')) then
            ucon = of()
            call rconn(ucon)
	 end if
	 if (find('rcrd')) then
            urcrd = of()
	 end if
	 if (find('wcrd'))then
            uwcrd = of()
	 end if
      end if
      XS = getd('conc',XS)
      irand = geti('rand',irand)
      massd = getd('masd',massd)
      goto 1

 2    continue

      if (massd.lt.0.d0) then
       write (*,*) 'INSERT THE PROPER MASS DENSITY FOR THE SOLUTION!'
       call alert(name,namel,'no mass density found!',22,1)
      endif

      if (XS .le. 0.d0) then
       write (*,*) 'mole fraction =',XS
       call alert(name,namel,'incorrect mole fraction',23,1)
      endif

      call RLUXGO(223,irand,0,0)

c     use getcrd for the solute
      call getcrd(urcrd,'CHARM')

      nptS = npt

      cm(1) = 0.d0
      cm(2) = 0.d0
      cm(3) = 0.d0

      do i=1,npt
       cm(1) = cm(1) + coor(1,i)
       cm(2) = cm(2) + coor(2,i)
       cm(3) = cm(3) + coor(3,i)
      enddo

      cm(1) = cm(1)/dble(npt)
      cm(2) = cm(2)/dble(npt)
      cm(3) = cm(3)/dble(npt)

      SoluteName = moname(poimon(1))

      do i=1,npt
       do l=1,3
        reference(l,i) = coor(l,i) - cm(l)
       enddo
       ptnmS(i) = ptnm(i)
      enddo

      xb=abs(xb)
      yb=abs(yb)
      zb=abs(zb)
c     here adjust the vdrw radius
      omass=15.999d0
      hmass=1.008d0
      hvdw=dsqrt(4.d0*0.0498)*0.30d0**(6)
      ovdw=dsqrt(4.d0*0.1521)*3.15d0**(6)
c
c
      write(stdo,'(/10x,a/)') "Reading water box ..."
      rewind(urwbx)
      read(urwbx,1000)junk
      read(urwbx,1000)junk
      read(urwbx,1000)junk
1000  format(a1)
      read(urwbx,'(i7)') nrec
      nrec=nrec/3
      write(stdo,'(/10x,a,i6,a/)') "Water box file has n=",nrec,
     $         " water molecules !!!"
c
      NW0 = 0
      do 625 i = 1,nrec
c         read(urwbx,'(20x,3f10.5)',end=630) xx,yy,zz
c         read(urwbx,'(20x,3f10.5)') xx1,yy1,zz1
c         read(urwbx,'(20x,3f10.5)') xx2,yy2,zz2
         read(urwbx,'(24x,3f10.5)') xx,yy,zz
         read(urwbx,'(24x,3f10.5)') xx1,yy1,zz1
         read(urwbx,'(24x,3f10.5)') xx2,yy2,zz2
         if (.not.((abs(xx).lt.xb).and.(abs(yy).lt.yb)
     &        .and.(abs(zz).lt.zb))) go to 625
         NW0 = NW0 + 1
         O(1,NW0) = xx - 0.5d0*xb
         O(2,NW0) = yy - 0.5d0*yb
         O(3,NW0) = zz - 0.5d0*zb
         H1(1,NW0) = xx1 - 0.5d0*xb
         H1(2,NW0) = yy1 - 0.5d0*yb
         H1(3,NW0) = zz1 - 0.5d0*zb  
         H2(1,NW0) = xx2 - 0.5d0*xb
         H2(2,NW0) = yy2 - 0.5d0*yb
         H2(3,NW0) = zz2 - 0.5d0*zb         
         remove(NW0) = 0
625    continue


       write (*,*) 'the water box has now',NW0,'molecules'
       
       NS = 0
       NRemWat = 0.d0
       npt = 0
       totmon = 0

       tol = XS*0.001d0

20     CONTINUE


       XSt = dble(NS)/(dble(NS) + NW0 - NRemWat)
       write (*,*) 'concentration',XSt
       write (*,*) 'number of solutes',NS
       write (*,*) 'number of water removed',NRemWat
       if (XSt.gt.(XS-tol) .and. XSt.lt.(XS+tol)) goto 21
       if (XSt.gt.XS) goto 21

       call RANLUX(random,3)
       cm(1) = (random(1)-0.5e0)*xb
       cm(2) = (random(2)-0.5e0)*yb
       cm(3) = (random(3)-0.5e0)*zb


       do i=1,nptS     
        do l=1,3 
         solute(l,i) = reference(l,i) + cm(l)
        enddo
       enddo

c check if it clashes with other solutes
       if (NS.gt.0) then
        e = 0.d0
        do i=1,totmon
         e = 0.d0
         do k=1,npts
          do j=1,npts
           dx = coor(1,(i-1)*npts+k) - solute(1,j)
           dy = coor(2,(i-1)*npts+k) - solute(2,j)
           dz = coor(3,(i-1)*npts+k) - solute(3,j)
           dx = dx - xb*anint(dx/xb)
           dy = dy - yb*anint(dy/yb)
           dz = dz - zb*anint(dz/zb)
           r = dsqrt( dx**2 + dy**2 + dz**2) 
           e = e + epsgm12(j)*
     1             ptvdw((i-1)*npts+k)/r**(12)
          enddo
         enddo
         if (e .gt. 25.d0) then
          write (*,*) 'solute rejected'
          goto 20
         endif
        enddo
       endif


       NS = NS + 1
       totmon = totmon + 1
       do i = 1, nptS
        coor(1,npt + i) = solute(1,i)
        coor(2,npt + i) = solute(2,i)
        coor(3,npt + i) = solute(3,i)
        ptms(npt + i) = ptms(i)
        ptnm(npt + i) = ptnmS(i)
        poimon(npt + i) = totmon
        id(npt + i) = i
        ptvdw(npt + i) = epsgm12(i)
        !write (111,*) npt,NS,i,id(npt+i)
       enddo 
       moname(totmon) = SoluteName
       poipt(totmon) = npt + nptS
       npt = npt + nptS


c now create the list of water to remove
       do 111 i=1,NW0
        if (remove(i).eq.1) goto 111
        e = 0.d0
        do j=1,nptS
         dx = O(1,i) - solute(1,j)
         dy = O(2,i) - solute(2,j)
         dz = O(3,i) - solute(3,j)
         dx = dx - xb*anint(dx/xb)
         dy = dy - yb*anint(dy/yb)
         dz = dz - zb*anint(dz/zb)
         r = dsqrt( dx**2 + dy**2 + dz**2 )
         e = e + epsgm12(j)*ovdw/r**(12)

         dx = H1(1,i) - solute(1,j)
         dy = H1(2,i) - solute(2,j)
         dz = H1(3,i) - solute(3,j)
         dx = dx - xb*anint(dx/xb)
         dy = dy - yb*anint(dy/yb)
         dz = dz - zb*anint(dz/zb)
         r = dsqrt( dx**2 + dy**2 + dz**2)
         e = e + epsgm12(j)*hvdw/r**(12)

         dx = H2(1,i) - solute(1,j)
         dy = H2(2,i) - solute(2,j)
         dz = H2(3,i) - solute(3,j)
         dx = dx - xb*anint(dx/xb)
         dy = dy - yb*anint(dy/yb)
         dz = dz - zb*anint(dz/zb)
         r = dsqrt( dx**2 + dy**2 + dz**2 )
         e = e + epsgm12(j)*hvdw/r**(12)     

        enddo
        if (e.gt.50.d0) then
         remove(i) = 1  
         NRemWat = NRemWat + 1.d0
         write (*,*) 'number of removed waters =',NRemWat
        endif       
111    continue      
       loop = loop + 1
       if (loop.gt.10000) stop
       goto 20

21     CONTINUE

c create the crd file
       do 112 i=1,NW0

        if (remove(i).eq.1) goto 112
        npt = npt + 1
        totmon = totmon + 1
        coor(1,npt) = O(1,i)
        coor(2,npt) = O(2,i)
        coor(3,npt) = O(3,i)
        poimon(npt) = totmon
        moname(totmon) = 'TIP3'
        ptnm(npt) = 'OH2 '
        ptms(npt) = omass
        ptvdw(npt) = ovdw

        npt = npt + 1
        coor(1,npt) = H1(1,i)
        coor(2,npt) = H1(2,i)
        coor(3,npt) = H1(3,i)
        poimon(npt) = totmon
        ptnm(npt) = 'H1  '
        ptms(npt) = hmass
        ptvdw(npt) = hvdw

        npt = npt + 1
        coor(1,npt) = H2(1,i)
        coor(2,npt) = H2(2,i)
        coor(3,npt) = H2(3,i)
        poimon(npt) = totmon
        ptnm(npt) = 'H2  '
        ptms(npt) = hmass
        ptvdw(npt) = hvdw
                   
        poipt(totmon) = npt
  
112    continue

      close(urwbx)

      write(stdo,'(/10x,a,i6/)') "Water box file was closed !!!"

      rewind(uwcrd)
c
c     The file was formed by putcrd so there are only 2 coment lines
c

c compute energy
      e = 0.d0
      do i=1,npt-1
       do j=i+1,npt
        if (poimon(i).eq.poimon(j)) goto 123
        dx = coor(1,i) - coor(1,j)
        dy = coor(2,i) - coor(2,j)
        dz = coor(3,i) - coor(3,j)
        dx = dx - xb*anint(dx/xb)
        dy = dy - yb*anint(dy/yb)
        dz = dz - zb*anint(dz/zb)
        r  = dsqrt(dx**2+dy**2+dz**2)
        e  = e + ptvdw(i)*ptvdw(j)/r**(12)
        if ((ptvdw(i)*ptvdw(j)/r**(12)).gt.100.d0) then
         write (*,*) e,r,i,j,ptnm(i),ptnm(j),ptvdw(i),ptvdw(j)
        endif
123    continue
       enddo
      enddo

      write (*,*) 'TOTAL REPULSIVE ENERGY = ',e

C change volume
      totmass = 0.d0
      do i=1,npt
       write (555,*) ptnm(i),moname(poimon(i)),ptms(i)
       totmass = totmass + ptms(i)
      enddo
      write (*,*) 'expected mass density',massd
      write (*,*) 'total mass',totmass
      write (*,*) 'measured mass density',totmass/(xb*yb*zb)*1.66d0
      ratio = massd/(totmass/(xb*yb*zb)*1.66d0)
      write (*,*) 'scaling volume by',ratio
      ratio = ratio**(1.d0/3.d0)
      xb = xb/ratio
      yb = yb/ratio
      zb = zb/ratio
      write (*,*) 'measured mass after scaling',
     1      totmass/(xb*yb*zb)*1.66d0
      write (*,*) 'SIZE OF THE BOX:'
      write (*,*) 'XB = ',xb
      write (*,*) 'YB = ',yb
      write (*,*) 'ZB = ',zb
c      do i=1,totmon
c       if (i.eq.1) kinit = 1
c       if (i.gt.1) kinit = poipt(i-1)+1
c       kend = poipt(i)
c       cm(1) = 0.d0
c       cm(2) = 0.d0
c       cm(3) = 0.d0
c       do k=kinit,kend
c        cm(1) = cm(1) + coor(1,i) 
c        cm(2) = cm(2) + coor(2,i)
c        cm(3) = cm(3) + coor(3,i)
c       enddo
c       cm(1) = cm(1)/(dble(kend-kinit+1))
c       cm(2) = cm(2)/(dble(kend-kinit+1))
c       cm(3) = cm(3)/(dble(kend-kinit+1))
c       coor(1,i) = coor(1,i) - cm(1) + cm(1)*ratio
c       coor(2,i) = coor(2,i) - cm(2) + cm(2)*ratio
c       coor(3,i) = coor(3,i) - cm(3) + cm(3)*ratio
c      enddo

      write(uwcrd,200) 
      write(uwcrd,201)
200   format('* mixture')
201   format('*')
      write(uwcrd,202) npt
202   format(1i7)

      do 800 i=1,npt
       write(uwcrd,203) i,poimon(i),moname(poimon(i)),
     1       ptnm(i),coor(1,i),coor(2,i),coor(3,i)
 800  continue

203   format(i7,i7,1x,a4,1x,a4,3(f10.5))

      close(uwcrd)

      write(uwpol,'(a/a,i7.7/a)')
     $     "~","MOLC=(MIXT)   #mon=",totmon,"~"
      write(uwpol,'(10(a4,1x))') (moname(i), i=1,totmon)
      write(uwpol,'(a4)') "*EOD"
      close(uwpol)
c
      stop
      end
c--------------------------------------------------------






