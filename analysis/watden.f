       program watden
C This program computes the density of water in 8 corners of the simulation
C box and compares it with the expected water number density, which is 0.033121/A^3.
C where A=Angstrom
C atomic mass of water = 18.0
C density of water at 20C 0.998 g/cm^3 (Nelson, Biological Physics)
C density of water at 20C 998g/dm^3/18.0=55.4M
C In this program we compute the number of water atoms per A^3.
C This number times 10^(27)/[6.02*10^(23)] gives the Molarity.
C The molarity times 18.0 gives the concentration in g/cm^3.
C Variables:
C xbox,ybox.zbox = size of the box
C xwat,ywat,zwat = size of the corners where the number of water particles is computed
C                  if a particle different from water is found here the program stops.
C vol = total volume of each one of the 8 corners
C nwat = number of water particles in the 8 corners
C den = density of water particles computed in the 8 corners (nwat/vol) [den]=A^(-3)
C wden = density of water particles expected 0.033121/A^3 unless changed in input
C mH20 = mass of water = 18.0
C dH20 = density of water = 0.998 g/cm^3
C dmH20 = density of water in M 55.4M=dH20*1000/mH20
C cntm = conversion factor from #pt/A^3 to Molarity
C cntd = cntm*18/1000 = conversion factor from #pt/A^3 to g/cm^3
C na = Avogadro Number /10^(23)

       implicit none

       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/COORD.BLOCK'
       include 'COMMON/LINE.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/MUTA.BLOCK'

       integer nwat,i,n,nstr,namel,ucon,ucrd,udyn,of,lnwat(8)
       integer inofrz,geti,nofreez(maxpt),rbin
       integer stdi,stdo,uden

       double precision xbox,ybox,zbox
       double precision xwat,ywat,zwat
       double precision wden,den,vol,getd,lden(8),den2
       double precision cntm,cntd,mH20,dH20,dmH20,na

       logical find,rcrd,rdyn

       character*6 name

       parameter (mH20=18.0d0,dH20=0.998d0,dmH20=dH20*1000.d0/mH20)
       parameter (na=6.02d0)
       parameter (cntm=10000.d0/na,cntd=cntm*mH20/1000.d0)

       stdi=5
       stdo=6
       rbin = 1

       name='watden'
       namel=6

       jnkf=25
       open(unit=jnkf,status="scratch")

C INITIALIZATION

       rcrd = .FALSE.
       rdyn = .FALSE.
 
       nstr = 1

       xbox=0.d0
       ybox=0.d0
       zbox=0.d0
 
       xwat=0.d0
       ywat=0.d0
       zwat=0.d0

       wden=0.033121d0

       muta = .FALSE.

1      continue
       call rline(name,namel,stdi)
       
       if (find('file')) then

       if (find('rcon')) then
        ucon=of()
       if (find('muta'))  muta=.TRUE.
        call rconn(ucon)
       else if (find('rcrd')) then
        ucrd=of()
        call getcrd(ucrd,'CHARM')
        rcrd = .TRUE.
        if (rcrd .and. rdyn) then
         call alert(name,namel,'crd or dcd?',11,1)
        endif
       else if (find('rdyn')) then
        udyn=of()
        rdyn = .TRUE.
        if (rcrd .and. rdyn) then
         call alert(name,namel,'crd or dcd?',11,1)
        endif
       else if (find('dens')) then
        uden=of()
       endif

       else if (find('acti')) then
          go to 2

       else

       xbox=getd('xbox',xbox)
       ybox=getd('ybox',ybox)
       zbox=getd('zbox',zbox)

       xwat=getd('xwat',xwat)
       ywat=getd('ywat',ywat)
       zwat=getd('zwat',zwat)

       wden = getd('wden',wden)

       nstr = geti('#str',nstr)

       endif
       go to 1
2      continue

       write (stdo,*) xbox,ybox,zbox
       write (stdo,*) xwat,ywat,zwat

       nwat = 0
       DO i=1,8
        lnwat(i)=0
       ENDDO

       inofrz = npt
       DO i=1,npt
        nofreez(i)=i
       ENDDO

       if (rcrd .and. (nstr.gt.1)) then
        call alert(name,namel,'only 1 str in crd files!',24,1)
       endif

       DO n=1,nstr

       if (rdyn) then
        rewind udyn
        call rdyncrd(udyn,n,inofrz,nofreez,rbin)
       endif

       DO i=1,npt

        if (moname(poimon(i)).eq.'TIP3' .or.
     1      moname(poimon(i)).eq.'SPCE') then
C        write (6,*) moname(poimon(i)), coor(1,i),coor(2,i),coor(3,i)
        endif

        if ((moname(poimon(i)).eq.'TIP3' .or.
     1       moname(poimon(i)).eq.'SPCE').and.ptnm(i).eq.'OH2') then
         if ( (dabs(coor(1,i)).gt.(xbox*0.5d0-xwat)) .and.
     1        (dabs(coor(2,i)).gt.(ybox*0.5d0-ywat)) .and.
     2        (dabs(coor(3,i)).gt.(zbox*0.5d0-zwat)) ) then

C          write (77,*) coor(1,i),coor(2,i),coor(3,i)

          nwat = nwat+1

          if (coor(1,i) .gt. 0.d0) then
           if (coor(2,i) .gt. 0.d0) then
            if (coor(3,i) .gt. 0.d0) then
C + + +
             lnwat(1) = lnwat(1) + 1
            else
C + + -
             lnwat(2) = lnwat(2) + 1
            endif
           else 
            if (coor(3,i) .gt. 0.d0) then
C + - +
             lnwat(3) = lnwat(3) + 1
            else
C + - -
             lnwat(4) = lnwat(4) + 1
            endif
           endif
          else
           if (coor(2,i) .gt. 0.d0) then
            if (coor(3,i) .gt. 0.d0) then
C - + +
             lnwat(5) = lnwat(5) + 1
            else
C - + -
             lnwat(6) = lnwat(6) + 1
            endif
           else 
            if (coor(3,i) .gt. 0.d0) then
C - - +
             lnwat(7) = lnwat(7) + 1
            else
C - - -
             lnwat(8) = lnwat(8) + 1
            endif
           endif
          endif

         endif
        else if (moname(poimon(i)).ne.'TIP3' .and.
     1           moname(poimon(i)).ne.'SPCE' .and.
     2           moname(poimon(i)).ne.'NAO' .and. 
     3           moname(poimon(i)).ne.'CL') then 
         if ( (dabs(coor(1,i)).gt.(xbox*0.5d0-xwat)) .and.
     1            (dabs(coor(2,i)).gt.(ybox*0.5d0-ywat)) .and.
     2            (dabs(coor(3,i)).gt.(zbox*0.5d0-zwat)) ) then
 
          write (stdo,*) 'particle different from water found!!!'
          write (stdo,*) 'name of the particle=',moname(poimon(i)) 
          call alert(name,namel,'prt != wat found',16,0)
         endif
        endif

       ENDDO

       ENDDO

       vol = xwat*ywat*zwat

       den = dble(nwat)/vol/8.d0/(dble(nstr))
 
       den2 = 0.d0

       DO i=1,8
        lden(i) = dble(lnwat(i))/vol/(dble(nstr))
        den2 = den2 + lden(i)**2
       ENDDO
        den2 = den2/8.d0

       write (uden,1000) den,(den2-den**2)**(0.5)/8.d0**(0.5)
       write (uden,1001) dmH20/cntm

1000   format ('density of water computed den [#pt/A^3]=',
     1          1f10.5,3x,'+/-',1f10.5)
1001   format ('density of water expected wden [#pt/A^3]=',1f10.5)

       write (uden,2000) den*cntm,(den2-den**2)**(0.5)/8.d0**(0.5)*cntm
       write (uden,2001) dmH20

2000   format ('density of water computed den [molar]=',
     1          1f10.3,3x,'+/-',1f10.3)
2001   format ('density of water expected wden [molar]=',1f10.1)

       write (uden,3000) den*cntd,(den2-den**2)**(0.5)/8.d0**(0.5)*cntd
       write (uden,3001) dH20

3000   format ('density of water computed den [g/cm^3]=',
     1          1f10.5,3x,'+/-',1f10.5)
3001   format ('density of water expected wden [g/cm^3]=',1f10.3)

       write (uden,*) 'If den >> wden the system is too squeezed'
       write (uden,*) 'If den << wden the box is not filled'

       write (uden,*) 'total number of particles found:',nwat

       write (uden,4000) dble(nwat)/8.d0
4000   format ('average number of particle per corner:',1f10.1)

       write (uden,*) 'local densities'
       write (uden,5000)
       write (uden,5001) lden(1),lnwat(1)
       write (uden,5002) lden(2),lnwat(2)
       write (uden,5003) lden(3),lnwat(3)
       write (uden,5004) lden(4),lnwat(4)
       write (uden,5005) lden(5),lnwat(5)
       write (uden,5006) lden(6),lnwat(6)
       write (uden,5007) lden(7),lnwat(7)
       write (uden,5008) lden(8),lnwat(8)

5000   format ('x y z   [#pt/A^3]   #pt')
5001   format ('+ + +',1x,1f10.5,2x,1i5)
5002   format ('+ + -',1x,1f10.5,2x,1i5)
5003   format ('+ - +',1x,1f10.5,2x,1i5)
5004   format ('+ - -',1x,1f10.5,2x,1i5)
5005   format ('- + +',1x,1f10.5,2x,1i5)
5006   format ('- + -',1x,1f10.5,2x,1i5)
5007   format ('- - +',1x,1f10.5,2x,1i5)
5008   format ('- - -',1x,1f10.5,2x,1i5)


       stop
       end

