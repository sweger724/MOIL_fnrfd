       program pressure

       implicit none

c common blocks

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/LD.BLOCK'
        include 'COMMON/PT.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PRESS.BLOCK'

c declare the variables you need
c loops
        integer i,istr
c pressure
c nmol            - number of (covalently) distinct molecules in the system
c patom(imol)     - pointer to the last particle of molecule imol
c pmol(iatm)      - pointer to the molecule to which atom iatm belongs
c xmol(imol),...  - coordinates of the CM of molecule imol 
c vxmol(imol),... - velocity of the CM of the molecule imol
c fxmol(imol),... - force on the CM of the molecule imol
c massmol(imol)   - mass of the molecule imol
c press           - istantaneus pressure
c ave_press       - time average of the pressure
        integer imol,nstru
        integer ipt,ib,jb,k,fixlist(maxpt)

c I/O input
        integer ucon,urcrd,urvel,of,namel,rbin

        double precision uave,ave_temp,corr,ap2
        double precision owat12,owat6
        double precision wden,awat,bwat,cwat
        double precision input_temp
        logical find
        character*8 name

c initialization

        stdi = 5
        stdo = 6
        rbin = 1
 
        name = 'pressure'
        namel = 8

        do i=1,maxpt
         patom(i) = 0
        enddo

        pressON = .true.
        ave_press = 0.d0
        ap2 = 0.d0

        call readinput(ucon,urcrd,urvel,name,namel,nstru,awat,bwat,cwat,
     1                 input_temp)

        call init_pressure()

c create the list of molecules 
c       call dmolecule(npt,nb,ib1,ib2,ptnm,moname,poimon,nmol,patom,pmol)
c        write (*,*) 'nmol=',nmol
c        do imol=1,nmol
c         write (*,*) 'patom(',imol,')=',patom(imol)
c        enddo
c        do i=1,npt
c         write (*,*) 'pmol(',i,')=',pmol(i)
c        enddo
c create the molecular masses array
c        call compute_massmol(npt,ptms,patom,massmol)
c        do imol=1,nmol
c         write (*,*) 'massmol(',imol,')=',massmol(imol)
c        enddo

        inofrz = npt

        DO i=1,inofrz
         nofreez(i) = i
        ENDDO

        uave = 0.d0
        ave_temp = 0.d0

        !do imol=1,nmol
        ! if (moname(imol).eq.'TBUT') write (111,*) imol,nmol,'TBUT'
        !enddo

        DO 1 istr=1,nstru

        if (urvel.ne.(-1)) then
        rewind urvel
        call rdyncrd(urvel,istr,inofrz,nofreez,rbin)

        do i=1,npt
         velo(1,i) = coor(1,i) 
         velo(2,i) = coor(2,i)
         velo(3,i) = coor(3,i)
c         write (6,*) 'velo',i,'=',velo(1,i),velo(2,i),velo(3,i)
        enddo
        else
        do i=1,npt
         velo(1,i) = 0.d0
         velo(2,i) = 0.d0
         velo(3,i) = 0.d0
        enddo
        endif

        rewind urcrd
        call rdyncrd(urcrd,istr,inofrz,nofreez,rbin)

c adjust the coordinates in such a way that there are no breaks
        k = 0
        do ipt=1,npt
         fixlist(ipt)=0
        enddo
        DO imol=1,nmol
         ipt = k + 1
         fixlist(ipt)=1
         DO WHILE (ipt.lt.patom(imol))
          do ib=1,nb
           if (ib1(ib).eq.ordLISTatm(ipt)) then
            if (fixlist(ipt).eq.1) then
             call checkANDmove(ib1(ib),ib2(ib),coor,a,b,c)
             fixlist(ib2(ib)) = 1
            else if (fixlist(ib2(ib)).eq.1) then
             call checkANDmove(ib2(ib),ib1(ib),coor,a,b,c)
             fixlist(ib1(ib)) = 1
            else
             write (*,*) 'ERROR IN FIXING!'
             write (*,*) 'ipt=',ipt
             write (*,*) 'ordLISTatm(ipt)=',ordLISTatm(ipt)
             write (*,*) 'ib2=',ib2(ib)
             stop
            endif
           endif
          enddo
          ipt = ipt + 1
         ENDDO
        k = patom(imol)
       ENDDO 

c create molecular positions and velocities
c        call atom2molecule(npt,nmol,patom,pmol,coor
c     &,velo,ptms,xmol,ymol,zmol,xrel,yrel,zrel
c     &,vxmol,vymol,vzmol,massmol)

c        write (*,*) 'coor',301,'=',coor(1,301),coor(2,301),coor(3,301)
c        write (*,*) 'coor',302,'=',coor(1,302),coor(2,302),coor(3,302)
c        write (*,*) 'coor',303,'=',coor(1,303),coor(2,303),coor(3,303)

c        write (*,*) xmol(pmol(301)),ymol(pmol(301)),zmol(pmol(301))

        if (esymyes) call squeeze(a,b,c)
      call wdyncrd(126,nstru,istr,inofrz,nofreez,1)
        call nbondm()
        if (esymyes) call syminit()
        if (specl) call nbondm_spcl()


c create molecular positions and velocities
        call atom2molecule(npt,nmol,patom,ordLISTatm,pmol,coor
     &,velo,ptms,xmol,ymol,zmol,xrel,yrel,zrel
     &,vxmol,vymol,vzmol,massmol)

        virial = 0.d0

        call eforce()
        call wener(stdo)

        uave = uave + e_total

        write (*,*) 'virial =',virial

c convert dpot to a proper force (units? sign?)

c        do i=1,npt
c         dpot(1,i) = -dpot(1,i)
c         dpot(2,i) = -dpot(2,i)
c         dpot(3,i) = -dpot(3,i)
c         write (6,*) 'dpot',i,'=',dpot(1,i),dpot(2,i),dpot(3,i)
c        enddo

c       DO imol = 1,nmol
c        write (*,*) 'mol position(',imol,')=',
c     1  xmol(imol),ymol(imol),zmol(imol)
c        write (*,*) 'mol velocity(',imol,')=',
c     1  vxmol(imol),vymol(imol),vzmol(imol)
c        write (*,*) 'mol force(',imol,')=',
c     1  fxmol(imol),fymol(imol),fzmol(imol)
c       ENDDO

c compute the pressure
        press = 0.d0
        if (urvel.ne.(-1)) then
        do imol=1,nmol
         press = press + massmol(imol)*
     1                   (vxmol(imol)**2 + vymol(imol)**2 +
     2                                     vzmol(imol)**2) 
        enddo
        else
         press = input_temp*3.d0*dble(nmol)*0.001987d0
        endif
        write (*,*) 'temp=',1.d0/3.d0*press/(dble(nmol))/0.001987d0
        temp = 1.d0/3.d0*press/(dble(nmol))/0.001987d0
        ave_temp=ave_temp + 1.d0/3.d0*press/(dble(nmol))/0.001987d0

        write (6,*) 'virial= ',virial,'trace= ',virXX+virYY+virZZ

        virial = virial + V_PIdirXX + V_PIdirYY + V_PIdirZZ +
     1                    V_PIrecXX + V_PIrecYY + V_PIrecZZ -
     2                    V_PIrelXX - V_PIrelYY - V_PIrelZZ +
     3                    V_PIcorXX + V_PIcorYY + V_PIcorZZ 

        write (6,*) '            virXX               virYY  
     1        virZZ'
        write (6,99) virXX,virYY,virZZ
        write (6,*) '            V_PIXX              V_PIYY       
     1        V_PIZZ'         
        write (6,100) V_PIdirXX,V_PIdirYY,V_PIdirZZ
        write (6,101) V_PIrecXX,V_PIrecYY,V_PIrecZZ
        write (6,102) V_PIrelXX,V_PIrelYY,V_PIrelZZ
        write (6,103) V_PIcorXX,V_PIcorYY,V_PIcorZZ
99      format('vdw  ',3(f20.5,1x))
100     format('dir  ',3(f20.5,1x))
101     format('rec  ',3(f20.5,1x))
102     format('rel  ',3(f20.5,1x))
103     format('cor  ',3(f20.5,1x))

        press = press + virial
        press = press/(3.d0*a*b*c)

        write (6,*) 'virial=',virial

        write (6,*) 'pressure',istr,press

        ave_press = ave_press + press
        ap2 = ap2 + press**2
 
        write (6,*) 'ave_press',istr,ave_press/dble(istr)

        write (6,*) 'compressibility factor=',
     1  press/dble(nmol)*(a*b*c)/temp/0.001987d0

1       CONTINUE


        write (6,*) 'density=',dble(nmol)/(a*b*c)
        write (6,*) 'density*=',dble(nmol)/(a*b*c)*
     1 (epsgm12(1)/epsgm6(1))
        ave_press = ave_press/dble(nstru)
        ap2 = ap2/dble(nstru) - ave_press**2
        ap2 = (ap2/dble(nstru))**(0.5d0)
        ave_temp = ave_temp/dble(nstru)
        uave = uave/dble(nstru)

        write (6,*) 'density g/cm3 =',
     1  massmol(1)*dble(nmol)/(a*b*c)*1.66053886d0

        write (6,*) 'average pressure =',ave_press,'+/-',ap2
        write (6,*)' ... in atm=',ave_press*68569.d0,'+/-',ap2*68569.d0
c compute the correction term
        pi = 3.14159265d0
        owat12 = 0.d0
        owat6 = 0.d0
        do i=1,npt
         if ( ptnm(i).eq.'OH2' .and. 
     1       (moname(poimon(i)).eq.'TIP3' .or.
     2        moname(poimon(i)).eq.'SPCE')    ) then
             owat12 = epsgm12(i)
             owat6  = epsgm6(i)
             goto 333
         endif 
        enddo
333     continue
c        write (*,*) pi,epsgm12(1),epsgm6(1),cutvdw2
        corr = 0.d0
c estimate the bulk water density in the simulation
c this is approximate... it is based upon the fact
c that the corners of the box are bulk water and no
c non-water particles are there
        wden = 0.d0
        do i=1,npt
         if (ptnm(i).eq.'OH2'.and.
     1      (moname(poimon(i)).eq.'TIP3'.or.
     2       moname(poimon(i)).eq.'SPCE')) then
          if (dabs(coor(1,i)-0.5d0*a).lt.(awat) .and.
     1        dabs(coor(2,i)-0.5d0*b).lt.(bwat) .and.
     2        dabs(coor(3,i)-0.5d0*c).lt.(cwat) ) then
              wden = wden + 1.d0
          endif
         endif
        enddo
        wden = wden/(8.d0*awat*bwat*cwat)
        do i=1,npt
        corr = corr-24.d0*pi*
     1         (2.d0/9.d0*epsgm12(i)*owat12/((cutvdw2**(0.5))**9)-
     2          1.d0/3.d0*epsgm6(i) *owat6 /((cutvdw2**(0.5))**3))
        enddo
        corr = (wden*1.d0/(a*b*c))/6.d0*corr
        write (6,*) 'correction =',corr
c        ave_press = ave_press - corr
        write (6,*) 'pressure - correction =',(ave_press-corr)*68569.d0
        write (6,*) 'average temperature =',ave_temp
        write (6,*) 'uave =',uave

        write (6,*) 'compressibility factor=',
     1  ave_press/dble(nmol)*(a*b*c)/ave_temp/0.001987d0
        write (6,*) 'abc',a,b,c
        write (6,*) 'nmol=',nmol

      stop
      end


c         subroutine dmolecule(npt,nb,ib1,ib2,ptnm,moname,poimon,
c     1                        nmol,patom,pmol)
c
c       implicit none
c
c      Program takes the bond information and classify
c      the particles as molecules 
c      the atoms covalently connected to each other defines th emolecule
c      patom(i) - is the last atom of the i th molecule 
c      hwat - integer flag to explicitly take into account that H1 and H2
c             of water molecules are bound to the oxygen.

c         integer ib1(*),ib2(*),ir(nb)
c         integer i,npt,nb,imol,nmol
c         integer patom(*),poimon(*),pmol(*)
c         integer hwat
c         character*4 moname(*),ptnm(*)
c
c          do i=1,npt
c          ir(i)=0
c          enddo
cc    first bond is assigned to first molecule 
c          ir(1)=1
c          imol=0
c
c          do i=1, nb
c          ir(ib2(i))=1
c          enddo
c
c          do i=1,npt
cc create exeption for water:
cc if mshk is used in building the connectivity file, then 
cc the bonds within the water particles are not reported.
cc to make this code working, we explicitly tell to the 
cc code to skip the hydrogens of the water, which otherwise
cc would be read as new molecules
c         hwat = 0
c         if ((moname(poimon(i)).eq.'TIP3' .or.
c    1         moname(poimon(i)).eq.'SPCE') .and.
c     2         ptnm(i)(1:1).eq.'H') hwat = 1
c          if(ir(i).eq.0 .and. hwat.eq.0) then
c          imol=imol +1
c          patom(imol)=i-1
c          endif
c          pmol(i) = imol+1
c          enddo
cc     the last atom of the last molecule is the last atom of all 
c          imol=imol +1
c          patom(imol)=npt
c
cc          do j=1,imol
cc          write(sto,*) j,patom(j)
cc          enddo
c         nmol = imol
c
c          return
c          end


c        subroutine compute_massmol(nmol,ptms,patom,massmol)
c
c        implicit none
c
c        integer iatm,imol,k,nmol,patom(*)
c        double precision ptms(*),massmol(*)
c
c        k = 1
c        DO imol=1,nmol
c         massmol(imol) = 0.d0
c         DO iatm=k,patom(imol)
c          massmol(imol) = massmol(imol) + ptms(iatm)
c         ENDDO
c         k = patom(imol)+1
c        ENDDO
c
c        return
c        end
c
         subroutine atom2molecule(npt,nmol,patom,ordLISTatm,pmol,coor
     &,velo,mass,xmol,ymol,zmol,xrel,yrel,zrel
     &,vxmol,vymol,vzmol,massmol)

        implicit none
 
c       subroutine atom2mole(npt,imol,patom,x,y,z
c     &,mass,massmol)
c      Program atom number,position,force and velocities 
c      and convert it to molecule , x(yz)mol,fx(fyfz)mol,vx(vyvz)mol
c      atomic arrays 
       double precision coor(3,*)
       double precision velo(3,*)
       double precision mass(*)
cc       molecular arrays 
       double precision xmol(nmol),ymol(nmol),zmol(nmol)
       double precision vxmol(nmol),vymol(nmol),vzmol(nmol)
       double precision xrel(*),yrel(*),zrel(*)
       double precision massmol(nmol)
c       inverse of molecular mass 

       integer patom(*),pmol(*),npt,nmol,ordLISTatm(*)

       integer i

c       coordinates          
        call a2m(0,patom(1),ordLISTatm,mass,coor,1,massmol(1),xmol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,coor,1,
     1           massmol(i),xmol(i))
        enddo

        call a2m(0,patom(1),ordLISTatm,mass,coor,2,massmol(1),ymol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,coor,2,
     1           massmol(i),ymol(i))
        enddo

        call a2m(0,patom(1),ordLISTatm,mass,coor,3,massmol(1),zmol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,coor,3,
     1           massmol(i),zmol(i))
        enddo

c       forces          
c        call a2mf(0,patom(1),mass,dpot,1,massmol(1),fxmol(1)) ! first one needs special treatment
c        do i=2,nmol
c       call a2mf(patom(i-1),patom(i),mass,dpot,1,massmol(i),fxmol(i))
c        enddo

c        call a2mf(0,patom(1),mass,dpot,2,massmol(1),fymol(1)) ! first one needs special treatment
c        do i=2,nmol
c        call a2mf(patom(i-1),patom(i),mass,dpot,2,massmol(i),fymol(i))
c        enddo

c        call a2mf(0,patom(1),mass,dpot,3,massmol(1),fzmol(1)) ! first one needs special treatment
c        do i=2,nmol
c        call a2mf(patom(i-1),patom(i),mass,dpot,3,massmol(i),fzmol(i))
c        enddo

c       velocities
        call a2m(0,patom(1),ordLISTatm,mass,velo,1,massmol(1),vxmol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,velo,1,
     1           massmol(i),vxmol(i))
        enddo

        call a2m(0,patom(1),ordLISTatm,mass,velo,2,massmol(1),vymol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,velo,2,
     1           massmol(i),vymol(i))
        enddo

        call a2m(0,patom(1),ordLISTatm,mass,velo,3,massmol(1),vzmol(1)) ! first one needs special treatment
        do i=2,nmol
        call a2m(patom(i-1),patom(i),ordLISTatm,mass,velo,3,
     1           massmol(i),vzmol(i))
        enddo

c relative coordinates
        do i=1,npt
         xrel(i) = coor(1,i) - xmol(pmol(i))
         yrel(i) = coor(2,i) - ymol(pmol(i))
         zrel(i) = coor(3,i) - zmol(pmol(i)) 
        enddo


        return
        end

        subroutine a2m(ibeg,iend,ordLISTatm,mass,parray,xyz,
     1                 massmol,marray)
c       sums up all elements of a molecule by multiplying with its mass

        implicit none

        integer xyz
        integer ordLISTatm(*)
        double precision mass(*),parray(3,*),marray
        double precision massmol
        integer i,ibeg,iend

        marray=0.d0
        do i=ibeg+1,iend
        marray=marray+parray(xyz,ordLISTatm(i)) * mass(ordLISTatm(i))
        !write (*,*) ibeg,iend,i,xyz,parray(xyz,i),mass(i),massmol
        enddo
c       divide to molecular mass 
        marray = marray/massmol

c        write(*,*) xxmol,mm,m+1
        return
        end

        subroutine a2mf(ibeg,iend,mass,parray,xyz,massmol,marray)
c       sums up all elements of a molecule by multiplying with its mass

        implicit none

        integer xyz
        double precision mass(*),parray(3,*),marray
        double precision massmol
        integer i,ibeg,iend

        marray=0.d0
        do i=ibeg+1,iend
        marray=marray+parray(xyz,i)
        enddo
c       divide to molecular mass 
        marray = marray

c        write(*,*) xxmol,mm,m+1
        return
        end
        subroutine readinput(ucon,urcrd,urvel,name,namel,nstru,
     1             awat,bwat,cwat,input_temp)

        implicit none

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
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'

        character*8 name
        character*4 coortyp
        double precision getd
        integer namel,level,nstru
        integer of,geti
        integer i,n,nstrub
        integer rbin
        logical find

        integer ucon,urcrd,urvel
        integer ipick(maxpt)
        integer ucon1,ucon2

        double precision awat,bwat,cwat,input_temp

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

        my_pe = 0
        num_pes = 1

        coortyp = 'CHAR'
        jnkf = 25
        open (unit=jnkf,status='scratch')

        lcent = .false.
        ctrue = .true.
        shift = .false.
        ucon  = -1
        urcrd  = -1
        urvel  = -1
        debug = .false.
        hydro_scale = 1.d0
        surften=0.005d0
        input_temp = 300.d0

        call init_ef()

        awat = 5.d0
        bwat = 5.d0
        cwat = 5.d0

1       continue

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
31        continue
         end if
         if (find('rcrd')) urcrd = of()
         if (find('rvel')) urvel = of()

C   read coarse grained model parameters (associated files)
         call Read_CG_input()

        end if

C   read coarse grained model parameters (except files)
        call Read_CG_input2()

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
56       continue
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

        input_temp = getd('temp',input_temp)
        gbsu    = geti('gbsu',gbsu)
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
        awat = getd('awat',awat)
        bwat = getd('bwat',bwat)
        cwat = getd('cwat',cwat)
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

c repulsive wall
         if (find('rwal')) then
          rwall = .TRUE.
          nwalls = geti('#wal',nwalls)
          w0(1) = getd('wpo1',w0(1))
          normw(1) = geti('nwa1',normw(1))
          w0(2) = getd('wpo2',w0(2))
          normw(2) = geti('nwa2',normw(2))
          w0(3) = getd('wpo3',w0(3))
          normw(3) = geti('nwa3',normw(3))
          w0(4) = getd('wpo4',w0(4))
          normw(4) = geti('nwa4',normw(4))
          w0(5) = getd('wpo5',w0(5))
          normw(5) = geti('nwa5',normw(5))
          w0(6) = getd('wpo6',w0(6))
          normw(6) = geti('nwa6',normw(6))
          weps = getd('weps',weps)
                  write(stdo,*) 'rwal = ',rwall
          do i=1,nwalls
           if (normw(i).eq.1) then
            write (stdo,1045) i
1045        format(1x,' Repulsive wall',1i5,' perpedicular to x axis')
           else if (normw(i).eq.2) then
            write (stdo,1046) i
1046        format(1x,' Repulsive wall',1i5,' perpedicular to y axis')
           else if (normw(i).eq.3) then
            write (stdo,1047) i
1047        format(1x,' Repulsive wall',1i5,' perpedicular to z axis')
           else
            call alert(name,namel,'wall out of 3 dimensions!',25,1)
           endif
          enddo
          do i=1,nwalls
           write (stdo,1048) i,normw(i),w0(i)
1048       format(1x, 'wall #',1i5,'is perpendicular to axis',1i5,
     1    'and located at ',f9.3)
           if (normw(i).eq.0)
     1      call alert(name,namel,'no norm for a wall!',19,1)
          enddo
          write (stdo,1049) weps
1049      format(1x, 'rep wall is weps/r^6 --- weps= ',f9.3)
          endif

c electric field
          if (find('efie')) then
           efield_yes = .TRUE.
           DV = getd('DelV',DV)
           efield = getd('elec',efield)
           nefield = geti('nefi',nefield)
           DVtype = geti('DVty',DVtype)
           if (efield.ne.0.d0 .and. DV .ne. 0.d0) then
            call alert(name,namel,'both E and DV defined!',22,1)
           else if (efield.eq.0.d0 .and. DV.eq.0.d0) then
            call alert(name,namel,'neither E nor DV defined!',25,1)
           endif
          endif

        if (find('acti')) go to 2
        go to 1
2       continue


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
21              continue
        end if

c jmjm
       if (ewaldyes) call ewald_init()
c end of jmjm
        if (vp_flag) call vp_init()
c end of vp

         if (eCGyes) then
           call CGinit()
         endif

        if (gbsabool) call make_rborn
        if(nmb.gt.0) then
         if (D(nmb).le.0.d0 .or. alpha(nmb).le.0.d0) then
          level=1
          call alert(name,namel,'Dmor or alph is 0',16,level)
         end if
        end if
 
        return
        end 

        subroutine checkANDmove(fix,move,coor,a,b,c)
        implicit none
        integer fix,move
        double precision coor(3,*)
        double precision a,b,c
c check X
        if ((coor(1,fix)-coor(1,move)).gt.0.5d0*a)
     1                   coor(1,move) = coor(1,move) + a
        if ((coor(1,fix)-coor(1,move)).lt.-0.5d0*a)
     1                   coor(1,move) = coor(1,move) - a
c check Y
        if ((coor(2,fix)-coor(2,move)).gt.0.5d0*b)
     1                   coor(2,move) = coor(2,move) + b
        if ((coor(2,fix)-coor(2,move)).lt.-0.5d0*b)
     1                   coor(2,move) = coor(2,move) - b
c check Z
        if ((coor(3,fix)-coor(3,move)).gt.0.5d0*c)
     1                   coor(3,move) = coor(3,move) + c
        if ((coor(3,fix)-coor(3,move)).lt.-0.5d0*c)
     1                   coor(3,move) = coor(3,move) - c
        return
        end
