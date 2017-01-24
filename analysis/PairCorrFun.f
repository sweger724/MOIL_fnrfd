C this program computes the g(r) between two sets of atoms picked in input
C and it also computes the integral of G(r) 4pir^2(g(r)-1)dr, useful for Kirkwood-Buff theory

       program PairCorrFun

       implicit none

           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/DEBUG.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/FREEZ.BLOCK'
           include 'COMMON/CONVERT.BLOCK'
           include 'COMMON/CCRD.BLOCK'
           include 'COMMON/SYMM.BLOCK'
           include 'COMMON/MUTA.BLOCK'
       
c --- I/O variables
c --- --- ucon: wcon file
c --- --- urcrd: crd file
c --- --- ugofr: file with gofr
c --- --- uint: file with integral of 4pir^2(g(r)-1)dr
c --- --- ugoft: file with integral as a function of time
c --- --- ugtcf: time correlation function of the integral 
c --- --- of: open file function
c --- --- namel: length of name
c --- --- name: name of file for alert
c --- --- geti: function to read an integer from input
c --- --- getd: function to read a double precision from input
c --- --- find: function to read input
c --- --- nshells: number of shells in which the sphere is divided
c --- --- rmax: maximum distance for gofr
c --- --- rbin: 1 if binary, 0 otherwise
c --- --- npick: type of particles picked
c --- --- ipick: npick if particle picked, 0 otherwise
c --- --- nstru: number of structures in the crd file
c --- --- maxstru: maximum number of structures
c --- --- ncor: maximum number of structures in correlation function
c --- --- dtstep: difference in time (in femtoseconds) between two structures
       integer ucon,urcrd,ugofr,uint,ugoft,ugtcf 
       integer of
       integer namel
       character*4 name
       integer geti
       double precision getd
       logical find
       integer nshells
       double precision rmax
       integer rbin
       integer npick,ipick(maxpt)
       integer nstru
       integer maxstruct
       parameter (maxstruct=1000000)
       integer ncor
       double precision dtstep
c --- indeces
       integer istru,i,j,k,l
c --- picked particles
c --- --- icentral: returns the particle number of the ith central atom selected
c --- --- ncentral: number of atoms selected to be central
c --- --- iperipheral: returns the particle number of the ith peripheral atom selected
c --- --- nperipheral: number of atoms selected to be peripheral
       integer icentral(maxpt),ncentral
       integer iperipheral(maxpt),nperipheral
c --- distances and gofr
c --- --- dx = xi - xj
c --- --- dy = yi - yj
c --- --- dz = zi - zj
c --- --- r = sqrt(dx**2 + dy**2 + dz**2)
c --- --- maxshells = maximum number of shells allowed
c --- --- rbound= lower boundary of a shell (smaller radius)
c --- --- rshell= distance between two shells = rmax/nshells
c --- --- gofr= g of r
c --- --- Gint=int dr 4pir**2(g(r)-1)
c --- --- Goft=V/Nper*(<Nper>-Vsphere/V Nper)-Gint at each instant of time
c --- --- GTCF=time correlation function of the G integral
c --- --- intCF= integral of TCF
       double precision dx,dy,dz,r
       integer maxshells
       parameter (maxshells = 1000)
       double precision rbound(maxshells)
       double precision rshell
       double precision gofr(maxshells),Gint(maxshells)
       double precision Goft(maxstruct),GTCF(maxstruct),intCF
c --- CELL variables
c --- --- indexes of cell, before or after PBC imposition
c --- --- MAXCELL: maximum number of cells per side
c --- --- MAXPERIPHERAL: maximum number of peripheral atoms in a cell
c --- --- NCellPerSide: number of cells per side
c --- --- CentCellID: given the number of central pt returns the x, y and z coordinate of cell
c --- --- PeriNptPerCell: returns the number of peripheral particles in one cell
c --- --- PeriCellList: has the list of peripheral particles belonging to that cell
c --- --- fCellList: flag, T if this method is used
       integer ix,iy,iz,jx,jy,jz,jxPBC,jyPBC,jzPBC,jlist
       integer MAXCELL,MAXPERIPHERAL
       parameter (MAXCELL=10)
       parameter (MAXPERIPHERAL=5000)
       integer NCellPerSide
       integer CentCellID(3,maxpt)
       integer PeriNptPerCell(MAXCELL,MAXCELL,MAXCELL)
       integer PeriCellList(MAXCELL,MAXCELL,MAXCELL,MAXPERIPHERAL)
       logical fCellList

c --- other variables
c --- --- pi: 3.1415...
c --- --- symm: true if symmetry is used
c --- --- vol: volume of the box
c --- --- volshell: volume of a spherical shell
c --- --- volsphere: volume of the sphere till rmax
c --- --- twait: time to wait before averaging
c --- --- dneffective: change in the number of peripheral atoms: 0 if they are a diffeferent species from central, -1 if the same
       integer dneffective
       !integer twait
       double precision pi
       parameter (pi = 3.1415926535d0)
       double precision vol,volshell,volsphere
       logical symm

       stdi = 5
       stdo = 6
       jnkf = 25
       open(unit=jnkf,status='scratch') 

       rbin = 1

       name = "gofr"
       namel = 4

       lpstr = 1

c default values
      nstru = 1
      twait = 0
      rmax = 0.d0
      nshells = 0
      ncentral = 0
      nperipheral = 0
      do i=1,maxpt
       icentral(i) = 0
       iperipheral(i) = 0
      enddo
      symm = .false.
      norew = .false.
      ncor = 0
      dtstep = 0.d0
      fCellList = .FALSE.


1     continue
      call rline(name,namel,stdi)
c --- I/O files
      if (find('file')) then
c --- --- wcon
       if (find('conn')) then
        muta = .false.
        if (find('muta')) muta = .true.
        ucon = of()
        call rconn(ucon)
c --- --- rcrd
       else if (find('rcrd')) then
        urcrd = of()
c --- --- gofr
       else if (find('gofr')) then
        ugofr = of()
c --- --- G = integral of gofr
       else if (find('gint')) then
        uint = of()
c --- --- G as a function of time
       else if (find('goft')) then
        ugoft = of()
c --- --- TCF of G
!       else if (find('gtcf')) then
!        ugtcf = of()
       endif
      else
       nstru = geti('#str',nstru)
       nshells = geti('#she',nshells)
       ncor = geti('#cor',ncor)
       dtstep = getd('step',dtstep)
       rmax = getd('rmax',rmax)
       twait = geti('twait',twait)
       if (find('norw')) norew = .true.
       if (find('cent')) then
        call pick (ipick,npick)
        do i=1,npt
         if (ipick(i).ne.0) then
          ncentral = ncentral + 1
          icentral(ncentral) = i
         endif
        enddo
       endif
       if (find('peri')) then
        call pick(ipick,npick)
        do i=1,npt
         if (ipick(i).ne.0) then
          nperipheral = nperipheral + 1
          iperipheral(nperipheral) = i
         endif
        enddo
       endif
       if (find('symm')) then
        symm = .true.
        a = getd('xtra',0.d0)
        b = getd('ytra',0.d0)
        c = getd('ztra',0.d0)
       endif
       if (find('cell')) fCellList = .TRUE.
       if (find('acti')) goto 2
      endif
      goto 1
2     continue

      if (nstru.gt.maxstruct) then
       call alert(name,namel,'#str>maxstr...',14,1)
      endif

      if (fCellList) then
       if ((a.ne.b) .or. (a.ne.c) .or. (b.ne.c)) then
        call alert(name,namel,'cell list only with square box',30,1)
       endif
       NCellPerSide = dint(a/rmax)
       if (NCellPerSide .lt. 3) then
        write (*,*) 'Cells work only if you make 27 or more'
        write (*,*) 'reduce the rmax or remove cell'
        call alert(name,namel,'not enough cells!',17,1)
       endif
       if (nperipheral.gt.MAXPERIPHERAL) then
        call alert(name,namel,'MAXPERIPHERAL too small',23,1)
       endif
      endif
      if (maxshells.lt.nshells) then
       call alert(name,namel,'maxshells too small',19,1)
      endif 

      do i=1,nshells
       gofr(i) = 0.d0
      enddo

      rshell = rmax/dble(nshells)
      rbound(1) = 0.d0
      do i=2,nshells
       rbound(i) = rbound(i-1) + rshell
      enddo

      vol = a*b*c
      volsphere = 4.d0/3.d0*pi*rmax**3

      inofrz=npt
      do i=1,npt
       nofreez(i)=i
      enddo

      dneffective = 0

c      DO i=1,ncentral
c       DO j=1,nperipheral
c        if (icentral(i).eq.iperipheral(j)) then
c         write (*,*) 'same chemical specie is considered for g(r)'
c         write (*,*) 'the effective number of peripheral atoms is'
c         write (*,*) 'reduced by 1'
c         dneffective = -1
c         goto 25
c        endif
c       ENDDO
c      ENDDO

25    continue
c --- loop on structures
      DO istru=1,nstru    
       if (.not.norew) rewind urcrd
       call rdyncrd(urcrd,istru,inofrz,nofreez,rbin) 
       Goft(istru) = 0.d0
       if (istru.gt.twait) then
       if (fCellList) then 
        call MakeCellList(coor,a,b,c,                         !input
     1                    icentral,ncentral,                  !input
     2                    iperipheral,nperipheral,            !input
     3                    NCellPerSide,MAXCELL,maxpt,         !input
     4                    CentCellID,PeriCellList,            !output
     5                    PeriNptPerCell)                     !output
       endif
c --- --- loop on central atoms
       DO i=1,ncentral
c --- --- loop on peripheral atoms
        if (.not.fCellList) then
c --- --- --- NO CELL LIST 
        DO 3 j=1,nperipheral
         if (icentral(i).eq.iperipheral(j)) goto 3
         dx = coor(1,icentral(i)) - coor(1,iperipheral(j))
         dy = coor(2,icentral(i)) - coor(2,iperipheral(j))
         dz = coor(3,icentral(i)) - coor(3,iperipheral(j))   
         if (symm) then
          dx = dx - a*anint(dx/a)
          dy = dy - b*anint(dy/b)
          dz = dz - c*anint(dz/c)
         endif
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r) 
         if (r.lt.0.5d0) then
          write (*,*) icentral(i),'and',iperipheral(j),'are at',r
          write (*,*) coor(1,icentral(i)),coor(2,icentral(i)),
     1                coor(3,icentral(i))
          write (*,*) coor(1,iperipheral(j)),coor(2,iperipheral(j)),
     1                coor(3,iperipheral(j))
          write (*,*) dx,dy,dz
          write (*,*) epsgm12(icentral(i))*epsgm12(iperipheral(j))/
     1 r**12 - epsgm6(icentral(i))*epsgm6(iperipheral(j))/r**6
          stop
         endif
         if (r.lt.rmax) then
          k = aint(r/rshell) + 1
          gofr(k) = gofr(k) + 1.d0
          Goft(istru) = Goft(istru) + 1.d0
         endif
3       CONTINUE
        else
c --- --- --- CELL LIST
        ix = CentCellID(1,i)
        iy = CentCellID(2,i)
        iz = CentCellID(3,i)
        DO jx=ix-1,ix+1
         jxPBC = jx
         if (jx.lt.1) jxPBC = NCellPerSide
         if (jx.gt.NCellPerSide) jxPBC = 1
         DO jy=iy-1,iy+1
          jyPBC = jy
          if (jy.lt.1) jyPBC = NCellPerSide
          if (jy.gt.NCellPerSide) jyPBC = 1
          DO jz=iz-1,iz+1
           jzPBC = jz
           if (jz.lt.1) jzPBC = NCellPerSide
           if (jz.gt.NCellPerSide) jzPBC = 1
           DO 4 jlist=1,PeriNptPerCell(jxPBC,jyPBC,jzPBC)
            j = PeriCellList(jxPBC,jyPBC,jzPBC,jlist)
            if (icentral(i).eq.iperipheral(j)) goto 4

         dx = coor(1,icentral(i)) - coor(1,iperipheral(j))
         dy = coor(2,icentral(i)) - coor(2,iperipheral(j))
         dz = coor(3,icentral(i)) - coor(3,iperipheral(j))
         if (symm) then
          dx = dx - a*anint(dx/a)
          dy = dy - b*anint(dy/b)
          dz = dz - c*anint(dz/c)
         endif
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         if (r.lt.0.5d0) then
          write (*,*) icentral(i),'and',iperipheral(j),'are at',r
          write (*,*) coor(1,icentral(i)),coor(2,icentral(i)),
     1                coor(3,icentral(i))
          write (*,*) coor(1,iperipheral(j)),coor(2,iperipheral(j)),
     1                coor(3,iperipheral(j))
          write (*,*) dx,dy,dz
          write (*,*) epsgm12(icentral(i))*epsgm12(iperipheral(j))/
     1 r**12 - epsgm6(icentral(i))*epsgm6(iperipheral(j))/r**6
          stop
         endif
         if (r.lt.rmax) then
          k = aint(r/rshell) + 1
          gofr(k) = gofr(k) + 1.d0
          Goft(istru) = Goft(istru) + 1.d0
         endif

4          CONTINUE
          ENDDO
         ENDDO
        ENDDO
        endif
       ENDDO
       Goft(istru) = Goft(istru)/dble(ncentral) - 
     1               volsphere/vol*dble(nperipheral+dneffective)
       Goft(istru) = Goft(istru)*vol/dble(nperipheral+dneffective)
       write (ugoft,'(1i10,1f12.5)') istru,Goft(istru)
       endif
      ENDDO
c --- normalize gofr
      DO i=1,nshells 
c --- --- normalize by the number of structures and the number of central atoms
       gofr(i) = gofr(i)/(dble(nstru-twait)*dble(ncentral))    
c --- --- normalize by the number of peripheral atoms that on average would be in the shell if evely distributed in the box
       if (i.lt.nshells) then
        volshell = 4.d0/3.d0*pi*(rbound(i+1)**3 - rbound(i)**3)
       else
        volshell = 4.d0/3.d0*pi*(rmax**3 - rbound(i)**3)
       endif
       gofr(i) = gofr(i)/(dble(nperipheral+dneffective)/vol*volshell)
       if (i.eq.1) then
        Gint(1) = (gofr(1)-1.d0)*
     1            4.d0*pi*((rbound(2))*0.5d0)**2*rshell
       else if (i.lt.nshells) then
        Gint(i) = Gint(i-1) + (gofr(i)-1.d0)*
     1            4.d0*pi*((rbound(i+1)+rbound(i))*0.5d0)**2*
     2            rshell
       else
        Gint(i) = Gint(i-1) + (gofr(i)-1.d0)*
     1           4.d0*pi*((rmax+rbound(i))*0.5d0)**2*
     2           rshell
       endif
       if (i.lt.nshells) then
 
       write (ugofr,'(1f10.5,1x,1f10.5)')(rbound(i)+rbound(i+1))*0.5d0,
     1                                    gofr(i)
       write (uint,'(1f10.5,1x,1f12.5)')(rbound(i)+rbound(i+1))*0.5d0,
     1                                    Gint(i) 
       else

       write (ugofr,'(1f10.5,1x,1f10.5)')(rbound(i)+rmax)*0.5d0,
     1                                    gofr(i)
       write (uint,'(1f10.5,1x,1f12.5)')(rbound(i)+rmax)*0.5d0,
     1                                    Gint(i)
       endif
      ENDDO

      !call CORFUN(ncor,nstru,dtstep,Goft,Goft,GTCF,intCF,.TRUE.)
      !DO i=1,ncor
      ! write (ugtcf,'(2f15.5)') dble(i)*dtstep,GTCF(i)
      !ENDDO
              
      stop
      end

      subroutine MakeCellList(coor,a,b,c,                         !input
     1                        icentral,ncentral,                  !input
     2                        iperipheral,nperipheral,            !input
     3                        NCellPerSide,MAXCELL,maxpt,         !input
     4                        CentCellID,PeriCellList,            !output
     5                        PeriNptPerCell)                     !output

      integer MAXCELL,maxpt
      integer icentral(maxpt),ncentral
      integer iperipheral(maxpt),nperipheral
      integer CentCellID(3,maxpt) 
      integer PeriCellList(MAXCELL,MAXCELL,MAXCELL,maxpt)
      integer PeriNptPerCell(MAXCELL,MAXCELL,MAXCELL)
      integer NCellPerSide
      double precision a,b,c
      double precision coor(3,maxpt)
   
      integer i,ipt,icellx,icelly,icellz,cnt

      !write (*,*) 'NCellPerSide',NCellPerSide

      do icellx=1,MAXCELL
       do icelly=1,MAXCELL
        do icellz=1,MAXCELL
         PeriNptPerCell(icellx,icelly,icellz) = 0
        enddo
       enddo
      enddo

      do i=1,ncentral
       ipt=icentral(i)
       icellx = dint((coor(1,ipt)+0.5d0*a)/a*dble(NCellPerSide)) + 1
       icelly = dint((coor(2,ipt)+0.5d0*b)/b*dble(NCellPerSide)) + 1
       icellz = dint((coor(3,ipt)+0.5d0*c)/c*dble(NCellPerSide)) + 1
       CentCellID(1,i) = icellx
       CentCellID(2,i) = icelly
       CentCellID(3,i) = icellz
      ! write (*,*) 'ID',CentCellID(1,i),CentCellID(2,i),CentCellID(3,i)
      enddo
 
      do j=1,nperipheral
       ipt=iperipheral(j)
       icellx = dint((coor(1,ipt)+0.5d0*a)/a*dble(NCellPerSide)) + 1
       icelly = dint((coor(2,ipt)+0.5d0*b)/b*dble(NCellPerSide)) + 1
       icellz = dint((coor(3,ipt)+0.5d0*c)/c*dble(NCellPerSide)) + 1
       PeriNptPerCell(icellx,icelly,icellz) = 
     1  PeriNptPerCell(icellx,icelly,icellz) + 1
       cnt = PeriNptPerCell(icellx,icelly,icellz)
       PeriCellList(icellx,icelly,icellz,cnt) = j
      ! write (*,*) 'LIST',icellx,icelly,icellz,cnt,j
      enddo

      return
      end
