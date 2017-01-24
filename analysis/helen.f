      program helen

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

      integer ucon,urcrd,uwlen,uwhis,namel,uwnrt
      integer imon,fmon,lmon,ipt,lpt,fpt
      integer i,j,k,of,nstru,istru,nbin,geti
      double precision s(3,maxmono),length,maxl,getd,dbin
      double precision hist(1000000),nhist,ds,CA(3,maxmono)
      integer iO(maxpt),iH(maxpt),Hcnt,Ocnt
      double precision d3,d4,nrpht,f(3),hbcut
      double precision tors,phi(maxpt),psi(maxpt)
      integer iC(maxpt),iCA(maxpt),iNa(maxpt),CAcnt,Ncnt,Ccnt
      integer iphi1(maxpt),iphi2(maxpt),iphi3(maxpt),iphi4(maxpt)
      integer ipsi1(maxpt),ipsi2(maxpt),ipsi3(maxpt),ipsi4(maxpt)
c nrpht = number of residues per helical turn
      logical find
      character*4 name

      stdi=5
      stdo=6

      namel=5
      name='helen'
c initialization
      norew=.false.
      nstru = 0
      nbin = 0
      maxl=0.d0
      hbcut=3.d0

c READ INPUT

c  open junk file for rline
c
      jnkf=25
      open(unit=jnkf,status='scratch')

1     continue
            call rline(name,namel,stdi)

            if (find('norw')) norew=.true.
C READ FILES
            if (find('file')) then

               if (find ('conn')) then
                ucon=of()
                call rconn(ucon)
               else if (find('rcrd')) then
                urcrd=of()
               else if (find('wlen')) then
                uwlen=of()
               else if (find('whis')) then
                uwhis=of()
               else if (find('wnrt')) then
                uwnrt=of()
               endif
C READ VARIABLES
             else
             nstru = geti('#str',nstru)
             maxl = getd('maxl',maxl)
             nbin = geti('#bin',nbin)
             fmon = geti('fmon',fmon)
             lmon = geti('lmon',lmon)
             hbcut= getd('hbcu',hbcut)
             if (find('acti')) goto 2
             endif
             goto 1
2     continue 
      if (fmon.lt.2) then
       call alert(name,namel,'fmon must be >2',15,1)
      endif
      if (lmon.gt.(totmon-1)) then
       call alert(name,namel,'lmon must be <totmon',20,1)
      endif 
c initi nofreez
             inofrz=npt
             do  i=1,npt
                nofreez(i)=i
             enddo

      dbin = maxl/dble(nbin)

      DO i=1,nbin
       hist(i) = 0.d0
      ENDDO

       rewind urcrd
C loop on structures
      DO 100 istru=1,nstru
       if (.not.norew) rewind urcrd
       call rdyncrd(urcrd,istru,inofrz,nofreez,1)
       DO imon=fmon,lmon
        fpt = poipt(imon)+1
        lpt = poipt(imon+1)
c        write (6,*) 'fpt=',fpt
c        write (6,*) 'lpt=',lpt
        iCA(1) = 0
        DO i=fpt,lpt
         if (ptnm(i).eq.'CA') then
          iCA(1) = i
          goto 101
         endif
c         write (6,*) imon,i,ptnm(i)
        ENDDO
        if (iCA(1).eq.0) call alert(name,namel,'no CA',5,1)
101     continue        
c        write (6,*) 'iCA(1)=',iCA(1)
        CA(1,imon-fmon+1) = coor(1,iCA(1))
        CA(2,imon-fmon+1) = coor(2,iCA(1))
        CA(3,imon-fmon+1) = coor(3,iCA(1))
c        write (777,'(1i5,3f10.5)') imon-fmon+1,
c     1  coor(1,iCA(1)),coor(2,iCA(1)),coor(3,iCA(1))
       ENDDO
       DO imon=1,lmon-fmon-2
        s(1,imon) = (CA(1,imon)+CA(1,imon+1)+CA(1,imon+2)+CA(1,imon+3))
        s(2,imon) = (CA(2,imon)+CA(2,imon+1)+CA(2,imon+2)+CA(2,imon+3))
        s(3,imon) = (CA(3,imon)+CA(3,imon+1)+CA(3,imon+2)+CA(3,imon+3))
        s(1,imon) = s(1,imon)*0.25d0
        s(2,imon) = s(2,imon)*0.25d0
        s(3,imon) = s(3,imon)*0.25d0
c        write (6,'(1i5,3f10.5)') imon,s(1,imon),s(2,imon),s(3,imon)
       ENDDO
       length = 0.d0
       DO imon=1,lmon-fmon-2
        ds = (s(1,imon+1)-s(1,imon))**2 + (s(2,imon+1)-s(2,imon))**2 +
     1       (s(3,imon+1)-s(3,imon))**2
        ds = ds**(0.5d0)
        length = length + ds 
       ENDDO
c       write (6,'(1f10.5)') length
       write (uwlen,'(1i10,1f10.5)') istru,length
       if (length.gt.maxl) then
        write (6,*) 'ERROR!'
        call alert (name,namel,'maxl<len',8,1)
       endif
       DO i=1,nbin
        if (length.ge.dble(i-1)*dbin .and. length.lt.dble(i)*dbin) then
         hist(i) = hist(i) + 1.d0
         nhist = nhist + 1.d0
         go to 105
        endif
       ENDDO
       
105   CONTINUE
      Hcnt=0
      Ocnt=0
      Ccnt=0
      CAcnt=0
      Ncnt=0
      DO i=1,totmon
       iH(i)=0
       iO(i)=0
       iNa(i)=0
       iC(i)=0
       iCA(i)=0
      ENDDO
      nrpht=0.d0
      f(1)=0.d0
      f(2)=0.d0
      f(3)=0.d0
      DO i=1,npt
       if (ptnm(i).eq.'H   ') then
        Hcnt = Hcnt + 1
        iH(Hcnt) = i
       else if (ptnm(i).eq.'O   ') then
        Ocnt = Ocnt + 1
        iO(Ocnt) = i
       else if (ptnm(i).eq.'C   ') then
        Ccnt = Ccnt + 1
        iC(Ccnt) = i
       else if (ptnm(i).eq.'CA  ') then
        CAcnt = CAcnt + 1
        iCA(CAcnt) = i
       else if (ptnm(i).eq.'N   ') then
        Ncnt = Ncnt + 1
        iNa(Ncnt) = i
       endif
      ENDDO
      DO i=1,totmon
       iphi1(i)=0
       iphi2(i)=0
       iphi3(i)=0
       iphi4(i)=0
       ipsi1(i)=0
       ipsi2(i)=0
       ipsi3(i)=0
       ipsi4(i)=0
      ENDDO
      DO i=2,totmon-1
         
       if (iC(i-1).ne.0 .and.
     1     iNa(i).ne.0 .and.
     1     iCA(i).ne.0 .and.
     1     iC(i).ne.0 ) then
        iphi1(i)=iC(i-1)
        iphi2(i)=iNa(i)
        iphi3(i)=iCA(i)
        iphi4(i)=iC(i)
        phi(i)=
     1  tors(coor,maxpt,iphi1(i),iphi2(i),iphi3(i),iphi4(i))
        write (777,*) phi(i)
       endif
       if (iNa(i).ne.0 .and.
     1     iCA(i).ne.0 .and.
     1     iC(i) .ne.0 .and.
     1     iNa(i+1).ne.0 ) then
        ipsi1(i)=iNa(i)
        ipsi2(i)=iCA(i)
        ipsi3(i)=iC(i)
        ipsi4(i)=iNa(i+1)
        psi(i)=
     1  tors(coor,maxpt,ipsi1(i),ipsi2(i),ipsi3(i),ipsi4(i))
        write (777,*) psi(i)
       endif
      ENDDO
c      write (777,*) 'nstru=',nstru
      DO 11 i=fmon,lmon-3
       if (.NOT.(phi(i)  .lt.0.d0 .and. phi(i+1).lt.0.d0 .and.
     1           phi(i+2).lt.0.d0 .and. phi(i+3).lt.0.d0 .and.
     2           psi(i)  .lt.0.d0 .and. psi(i)  .gt.-170.d0 .and.
     2           psi(i+1).lt.0.d0 .and. psi(i+1).gt.-170.d0 .and.
     2           psi(i+2).lt.0.d0 .and. psi(i+2).gt.-170.d0 .and.
     2           psi(i+3).lt.0.d0 .and. psi(i+3).gt.-170.d0)) then
        f(1) = f(1) + 1.d0
        goto 11
       endif
       d3 = (coor(1,iH(i+3))-coor(1,iO(i)))**2 +
     1      (coor(2,iH(i+3))-coor(2,iO(i)))**2 +
     2      (coor(3,iH(i+3))-coor(3,iO(i)))**2
       d3 = d3**(0.5)
       d4 = (coor(1,iH(i+4))-coor(1,iO(i)))**2 +
     1      (coor(2,iH(i+4))-coor(2,iO(i)))**2 +
     2      (coor(3,iH(i+4))-coor(3,iO(i)))**2
       d4 = d4**(0.5)
       write (999,*) d3,d4
       if (d3.lt.d4) nrpht = nrpht + 3.d0
       if (d4.lt.d3) nrpht = nrpht + 3.6d0 
       if (d3.gt.hbcut .and. d4.gt.hbcut) f(1)=f(1)+1.d0
       if (d3.lt.d4 .and. d3.lt.hbcut) f(2)=f(2)+1.d0
       if (d4.lt.d3 .and. d4.lt.hbcut) f(3)=f(3)+1.d0
c       write (777,*) iO(i),iH(i+3),d3
c       write (777,*) iO(i),iH(i+4),d4
11    CONTINUE
      nrpht = nrpht/dble(lmon-3-fmon+1)
      f(1) = f(1)/dble(lmon-3-fmon+1)
      f(2) = f(2)/dble(lmon-3-fmon+1)
      f(3) = f(3)/dble(lmon-3-fmon+1)
      write (uwnrt,'(1i10,4f10.5)') istru,nrpht,f(1),f(2),f(3)
100   CONTINUE

      DO i=1,nbin
       write (uwhis,'(2f10.5)') dble(2*i-1)*dbin*0.5d0,hist(i)/nhist
      ENDDO

      stop
      end

      double precision function tors(coor,maxpt,i1,i2,i3,i4)
      integer maxpt,i1,i2,i3,i4
      double precision coor(3,maxpt),pi180
      double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
      double precision ux,uy,uz,vx,vy,vz,uu,vv,uv
      pi180 = 180.d0/(4.d0*datan(1.d0))
c      write (777,*) coor(1,i1),coor(1,i2),coor(1,i3),coor(1,i4)
      dx1 = coor(1,i2) - coor(1,i1)
      dy1 = coor(2,i2) - coor(2,i1)
      dz1 = coor(3,i2) - coor(3,i1)

      dx2 = coor(1,i3) - coor(1,i2)
      dy2 = coor(2,i3) - coor(2,i2)
      dz2 = coor(3,i3) - coor(3,i2)

      dx3 = coor(1,i4) - coor(1,i3)
      dy3 = coor(2,i4) - coor(2,i3)
      dz3 = coor(3,i4) - coor(3,i3)

      ux  = dy1*dz2 - dz1*dy2
      uy  = dz1*dx2 - dx1*dz2
      uz  = dx1*dy2 - dy1*dx2

      vx  = dy2*dz3 - dz2*dy3
      vy  = dz2*dx3 - dx2*dz3
      vz  = dx2*dy3 - dy2*dx3

      uu  = (ux*ux+uy*uy+uz*uz)
      vv  = (vx*vx+vy*vy+vz*vz)
      uv  = (ux*vx+uy*vy+uz*vz)/dsqrt(uu*vv)

      tors = dacos(uv)*pi180

      dx1 = uy*vz - uz*vy
      dy1 = uz*vx - ux*vz
      dz1 = ux*vy - uy*vx

      if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) tors = - tors
c               phi = phi - 180.d0
      if (phi.lt.-180) tors = tors + 360
        write(777,'(4i5,1f10.5)') i1,i2,i3,i4,tors
      end function     
