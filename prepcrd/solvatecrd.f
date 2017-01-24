      Program SolvateCRD
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
c
      character*4 styl1,styl2
      character*4 name
      integer namel
      integer ucon,urcrd,uwcrd,urwbx,uwpol,uwcon
c
      integer geti,of,ierr
      double precision getd
      logical find,fopen
c
      character*6 getchar
      integer level,lc
      integer lpst,lpend
      character*5 ctyp
      data ucon,urcrd,uwcrd,urwbx,uwpol,uwcon/6*99/
c
      double precision xx,yy,zz,cg,cg1,cg2,cg3,ww,xb2,yb2,zb2
      double precision xx1,xx2,yy1,yy2,zz1,zz2,box,xb,yb,zb,two
c     zva : exclusion region from centered box
      double precision rvdw
      double precision xbex,ybex,zbex
      double precision totmass
      integer i,n,j,numwats,numatoms,nrec
      integer nmcon,nacon
      integer nmon1,nmon2,nmon3,na1,na2,na3
      character*4 satom
      character*2 a1,a2,a3
      integer maxlines,maxblock,nsym,nsymb,maxp
      parameter(maxlines=100000,maxblock=10,nsym=10,nsymb=5)
      parameter(rvdw=3.0)
      integer natom,nmono,numlines
      integer nblock,icount,iblock,k,icter,ipoly,iatom
      double precision totmas,xcm,ycm,zcm,ro2,rh2,ovdw
      double precision x,y,z,hvdw,rd2,crdcent(3)
      character*6 sym(nsym),ban(nsymb)
      character*1 db
      character atom*6,wat*3,elem*4,mono*4,prvmono*4,last*14,long*80
      integer indmono(maxpt)
      character*10 crdline1(maxpt),crdline2(maxpt)
      character*4 mononame
      character*80 fileicrd, fileocrd,filecon,filewat,filepol
      logical flag,ok,select
      integer npick,pointr(maxpt),ipick(maxpt)
      character*1 junk

      data sym(1),sym(2),sym(3),sym(4),sym(5),sym(6),sym(7),sym(8),
     &     sym(9),sym(10)/
     &     'ATOM  ','HETATM','NEXT  ','NEXT  ','NEXT  ','NEXT  '
     &     ,'NEXT  ','NEXT  ','NEXT  ','NEXT  '/
      data ban(1),ban(2),ban(3),ban(4),ban(5)/
     &     'CONECT','TER   ','      ','MASTER','END   '/
c
      data xb,yb,zb/3*1.0d0/
      data xbex,ybex,zbex/3*0.0d0/


c
c
c     General initialization
c
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

      select = .false.

      jnkf = 25
      open (unit=jnkf,status='scratch')

 1    continue

      call rline(name,namel,stdi)
      if (find('acti')) goto 2
      xb     = getd('xwbx',  xb  )
      yb     = getd('ywbx',  yb  )
      zb     = getd('zwbx',  zb  )
      xbex   = getd('xbex',xbex)
      ybex   = getd('ybex',ybex)
      zbex   = getd('zbex',zbex)
c
      if (find('debu')) debug = .true.
      if (find('file')) then
	 if (find('wcon')) then
            uwcon = of()
	 end if
	 if (find('wpol')) then
            uwpol = of()
	 end if
	 if (find('rwbx')) then
            urwbx = of()
	 end if
	 if (find('conn')) then
            ucon = of()
            call rconn(ucon)
	 end if
	 if (find('rcrd')) then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','CHARM',lc)
            urcrd = of()
	 end if
	 if (find('wcrd'))then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','CHARM',lc)
            uwcrd = of()
ccc   call putcrd(uwcrd,'CHARM')
	 end if
      end if
      if (find('selc')) then
         select = .true.
         call pick(ipick,npick)
         npick=0
         do 3 i=1,npt
            npick=npick+ipick(i)
 3       continue
         j=0
         do 4 i=1,npt
            if (ipick(i).gt.0) then
               j=j+1
               pointr(j)=i
            end if
 4       continue
      endif

      goto 1

 2    continue



      if ( .not. select ) then
         do 303 i=1,npt
            pointr(i) = i
 303     continue
         npick = npt
      endif




c     use getcrd for the solute
      call getcrd(urcrd,ctyp)


c
      do 80 n=1,3
         crdcent(n)=0.0
 80   continue
      totmass = 0.d0
      do 100 j=1,npick
         i=pointr(j) 
         totmass = totmass + ptms(i)
         do 90 n=1,3
            crdcent(n)= crdcent(n) + coor(n,i) * ptms(i)
 90      continue
 100  continue

      do 180 n=1,3
         crdcent(n)= crdcent(n)/totmass
 180  continue

      if (xbex.ne.0.0d0.or.ybex.ne.0.0d0.or.zbex.ne.0.0d0) then
       write(stdo,*) 'THE BOX CENTER IS x=',xbex,'y=',ybex,'z=',zbex
     
      do i=1,npt
            coor(1,i) = coor(1,i) - xbex
            coor(2,i) = coor(2,i) - ybex
            coor(3,i) = coor(3,i) - zbex
      enddo
      else
    

      do 200 i=1,npt
         do 190 n=1,3
            coor(n,i) = coor(n,i) - crdcent(n)
 190     continue
 200  continue
      endif

      open(unit=10,file='solvtem',status='unknown')
      rewind 10

      call putcrd(10,'CHARM')


c     now get the box sizes and prepare waters
      xb=abs(xb)
      yb=abs(yb)
      zb=abs(zb)
      two=2.d0
      xb2=xb/two
      yb2=yb/two
      zb2=zb/two
c     here adjust the vdrw radius
      hvdw=1.d0*rvdw
      ovdw=hvdw+0.2
      ovdw=ovdw*ovdw
      hvdw=hvdw*hvdw
c
      natom=npt
      nmono=totmon
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
      do 625 i = 1,nrec
c         read(urwbx,'(20x,3f10.5)',end=630) xx,yy,zz
c         read(urwbx,'(20x,3f10.5)') xx1,yy1,zz1
c         read(urwbx,'(20x,3f10.5)') xx2,yy2,zz2
         read(urwbx,'(24x,3f10.5)',end=630) xx,yy,zz
         read(urwbx,'(24x,3f10.5)') xx1,yy1,zz1
         read(urwbx,'(24x,3f10.5)') xx2,yy2,zz2
         if (.not.((abs(xx).lt.xb).and.(abs(yy).lt.yb)
     &        .and.(abs(zz).lt.zb))) go to 625
         xx=xx-xb2
         yy=yy-yb2
         zz=zz-zb2
         xx1=xx1-xb2
         yy1=yy1-yb2
         zz1=zz1-zb2
         xx2=xx2-xb2
         yy2=yy2-yb2
         zz2=zz2-zb2

         do 623 j=1,npt
            x=coor(1,j)
            y=coor(2,j)
            z=coor(3,j)
            ro2=(xx-x)**2+(yy-y)**2+(zz-z)**2
            if (ro2.lt.ovdw) go to 625
c     avoid hydrogens when making decision if preferable
            rh2=(xx1-x)**2+(yy1-y)**2+(zz1-z)**2
            if (rh2.lt.ovdw) go to 625
            rd2=(xx2-x)**2+(yy2-y)**2+(zz2-z)**2
            if (rd2.lt.ovdw) go to 625
 623     continue

         nmono=nmono+1
c
         natom=natom+1
         write(10,'(2i7,a10,3f10.5,a10)')
c        write(10,'(2i6,a10,3f10.5,a10)')
     >        natom,nmono," TIP3 OH2 ",xx,yy,zz
         natom=natom+1
         write(10,'(2i7,a10,3f10.5,a10)')
c        write(10,'(2i6,a10,3f10.5,a10)')
     >        natom,nmono," TIP3 H1  ",xx1,yy1,zz1
         natom=natom+1
         write(10,'(2i7,a10,3f10.5,a10)')
c        write(10,'(2i6,a10,3f10.5,a10)')
     >        natom,nmono," TIP3 H2  ",xx2,yy2,zz2

 625  continue
 630  continue
c
      close(urwbx)
c
c     write(stdo,'(/10x,a,i6/)') "Water box file was closed !!! rec=",j
      write(stdo,'(/10x,a,i6/)') "Water box file was closed !!!"
c
      rewind(10)
      rewind(uwcrd)
c
c     The file was formed by putcrd so there are only 2 coment lines
c
 700  continue
      read(10,'(a80)',end=810) long
      if(long(1:1).ne.'*') go to 750
      write(uwcrd,'(a80)') long
      goto 700
c
 750  write(uwcrd,'(i7)') natom
c750  write(uwcrd,'(i6)') natom
c
      do 800 i=1,natom
         read(10,'(a80)',end=810) long
         write(uwcrd,'(a80)') long
 800  continue
c
 810  close(10)
      close(uwcrd)
c
      numwats=(natom-npt)/3
      mononame=BULK(1)
      write(uwpol,'(a/a,a4,a,i5.5/a)')
     $     "~","MOLC=(",mononame,")   #mon=",nmono,"~"
      write(uwpol,'(10a5)') (moname(i), i=1,totmon)
      write(uwpol,'(10a5)') ("TIP3", i=1,numwats)
      write(uwpol,'(a4)') "*EOD"
      close(uwpol)
c
      write(stdo,'(/3(3x,a,f6.1)/)') ' xb=',xb,' yb=',yb,
     >     ' zb=',zb
      write(stdo,'(//3x,a,i7//3x,a,i7/)') " natom=",natom,
c     write(stdo,'(//3x,a,i6//3x,a,i6/)') " natom=",natom,
     >     " nmono=",nmono
      stop
      end
c--------------------------------------------------------






