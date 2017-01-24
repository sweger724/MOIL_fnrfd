C This program returns a crd and a poly file with the mutant system
C INPUT:
C       *) native wcon and crd AND
C          (number of the monomer where the mutant should be attached OR
C          set of atoms)
C       *) (mutant wcon and crd AND
C           monomer to attach to the native OR 
C           set of atoms) OR
C          (entry for a library)
C
C OUTPUT: 
C       *) mutant PDB
C
       program PMNCrdPoly

       implicit none

       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/UNITS.BLOCK'
       include 'COMMON/DEBUG.BLOCK'
       include 'COMMON/LINE.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/MONOMERS.BLOCK'
       include 'COMMON/COORD.BLOCK'
       include 'COMMON/OVERLAP.BLOCK'
       include 'COMMON/MUTA.BLOCK'

c indexes 
       integer i,j,k,l

c units
       integer urncon,urmcon,urncrd,urmcrd,urlib,
     1         uwpol,uwcrd

c differnt types of input
       logical AminoAcid_yes,NucleicAcid_yes,
     1         OtherMoiety_yes
       logical MutConCrd_yes,MutRotLib_yes
       logical NatMuta,MutMuta

c variables and arrays for Nat and Mut
       integer iNatCA,iMutCA
       integer NatNpt,MutNpt
       integer NatTotmon,MutTotmon
       integer NatPoimon(maxpt),MutPoimon(maxpt)
       integer NatAlign(maxpt),MutAlign(maxpt)
       integer NatPoipt(maxmono),MutPoipt(maxmono)
       character*4 NatPtnm(maxpt),MutPtnm(maxpt)
       character*4 NatMoname(maxmono),MutMoname(maxmono)
       character*4 MutName,MutMonomerName
       double precision NatCoor(3,maxpt),MutCoor(3,maxpt)

c selction and rotation
       integer NatPick(maxpt),MutPick(maxpt)
       integer NNatPick,NMutPick,iMutConn,iNatConn
       double precision NatMass(maxpt),MutMass(maxpt)
       double precision NatCoorAlign(3,maxpt),MutCoorAlign(3,maxpt)
       double precision rms

c variables for reading input
       integer ipick(maxpt),npick
       integer of
       logical find
       character*4 getchar

c variable for alert
       integer namel
       character*10 name

c standard input and output
       stdi=5
       stdo=6

c open junk file for rline
       jnkf=25
       open(unit=jnkf,status='scratch')

c name for alert
       name='PMNCrdPoly'
       namel=10

c initialize the variables
       do i=1,maxpt
        ipick(i) = 0
       enddo
       npick = 0
c --- options for generation of the mutant crd
       MutConCrd_yes = .false.
       MutRotLib_yes = .false.
c --- options for typology of substituent
       AminoAcid_yes = .false.
       NucleicAcid_yes = .false.
       OtherMoiety_yes = .false.
c --- options to read the Nat and Mut conn files
       NatMuta = .false.
       MutMuta = .false.
c --- default names
       MutName = "XXXX"
       MutMonomerName = "XXXX"

c start reading
1      continue
       call rline(name,namel,stdi)
c --- read the files
       if (find('file')) then   
        if (find('ncon')) then
         if (find('muta')) NatMuta = .true.
         urncon = of() 
         muta = NatMuta
         call rconn(urncon)
c --- --- store Native info
	         NatNpt = npt
	         NatTotmon = totmon
	         do i=1,NatNpt
	          NatPtnm(i) = ptnm(i)
	          NatPoimon(i) = poimon(i)
	         enddo
	         do i=1,NatTotmon
	          NatMoname(i) = moname(i)
	          NatPoipt(i) = poipt(i)
	         enddo
        muta = .false.
        else if (find('mcon')) then
         if (find('muta')) MutMuta = .true.
         muta = MutMuta
         MutConCrd_yes = .true.
         urmcon = of()
         call rconn(urmcon)
c --- --- store Mutant info
	         MutNpt = npt
	         MutTotmon = totmon
	         do i=1,MutNpt
	          MutPtnm(i) = ptnm(i)
	          MutPoimon(i) = poimon(i)
	         enddo
	         do i=1,MutTotmon
	          MutMoname(i) = moname(i)
	          MutPoipt(i) = poipt(i)
	         enddo
        muta = .false.
        else if (find('ncrd')) then
         urncrd = of()
c --- --- copy back Native info
	 npt = NatNpt
	 totmon = NatTotmon
	 do i=1,NatNpt
	  ptnm(i) = NatPtnm(i)
	  poimon(i) = NatPoimon(i)
	 enddo
	 do i=1,totmon
	  moname(i) = NatMoname(i)
	  poipt(i) = NatPoipt(i)
	 enddo
         call getcrd(urncrd,'CHARM')
c --- --- store Native crd
	         do i=1,npt
	          NatCoor(1,i) = coor(1,i)
	          NatCoor(2,i) = coor(2,i)
	          NatCoor(3,i) = coor(3,i)
	         enddo
        else if (find('mcrd')) then
         urmcrd = of()
c --- --- copy back Mutant info
	 npt = MutNpt
	 totmon = MutTotmon
	 do i=1,MutNpt
	  ptnm(i) = MutPtnm(i)
	  poimon(i) = MutPoimon(i)
	 enddo
	 do i=1,totmon
	  moname(i) = MutMoname(i)
	  poipt(i) = MutPoipt(i)
	 enddo
         call getcrd(urmcrd,'CHARM')
c --- --- store Mutant crd
	         do i=1,npt
	          MutCoor(1,i) = coor(1,i)
	          MutCoor(2,i) = coor(2,i)
	          MutCoor(3,i) = coor(3,i)
	         enddo
        else if (find('mlib')) then
         urlib = of()
         if (MutConCrd_yes) then 
          call alert(name,namel,'cannot have crd and library',27,1)
         endif
         MutRotLib_yes = .true.
        else if (find('wcrd')) then
         uwcrd = of()
        else if (find('wpol')) then
         uwpol = of()
        endif
c --- read the variables  
c --- --- pick the native atoms
       else if (find('Nati')) then
c --- --- --- copy back the Native info
		             npt = NatNpt
		             totmon = NatTotmon
		             do i=1,npt
		              ptnm(i) = NatPtnm(i)
		              poimon(i) = NatPoimon(i)
		             enddo
		             do i=1,totmon
		              moname(i) = NatMoname(i)
		              poipt(i) = NatPoipt(i)
		             enddo
        call pick(ipick,npick)
        do i=1,npt
         NatPick(i) = ipick(i)
        enddo
       else if (find('Muta')) then
c --- --- --- copy back the Mutant info
		             npt = MutNpt
		             totmon = MutTotmon
		             do i=1,npt
		              ptnm(i) = MutPtnm(i)
		              poimon(i) = MutPoimon(i)
		             enddo
		             do i=1,totmon
		              moname(i) = MutMoname(i)
		              poipt(i) = MutPoipt(i)
		             enddo
        call pick(ipick,npick)
        do i=1,npt
         MutPick(i) = ipick(i)
        enddo
        else if (find('AmAc')) then
         AminoAcid_yes = .true.
        else if (find('NuAc')) then
         NucleicAcid_yes = .true.
        else if (find('OtMo')) then
         OtherMoiety_yes = .true. 
        else
        MutMonomerName=getchar('MMoN',MutMonomerName,4)
        MutName=getchar('MuNa',MutName,4)
c --- --- stop reading
        if (find('acti')) goto 2
        endif
        goto 1
2      continue

c check that proper names were introduced
       if (MutName.eq.'XXXX') then
        call alert(name,namel,'no mut name found in input',26,1) 
       endif
       if (MutMonomerName.eq.'XXXX') then
        call alert(name,namel,'no mut mol name found in input',30,1)
       endif

c check that a valid option is selected

       if ((.not.AminoAcid_yes).and.(.not.NucleicAcid_yes).and.
     1     (.not.OtherMoiety_yes)) then
        call alert(name,namel,'no moiety selected!',19,1)
       endif
       if (AminoAcid_yes.and.(NucleicAcid_yes.or.OtherMoiety_yes)) then
        call alert(name,namel,'too many moieties!',18,1)
       endif
       if (NucleicAcid_yes.and.OtherMoiety_yes) then
        call alert(name,namel,'too many moieties!',18,1)
       endif

       if (OtherMoiety_yes) then
        call alert(name,namel,'Other Moiety not funct',22,1)
       endif

c       if (NucleicAcid_yes) then
c        call alert(name,namel,'Nuclein Acid not funct',22,1)
c       endif

c select atoms for alignment
       if (AminoAcid_yes) then
c --- if AminoAcid select N, CA, CB, C 
        call findNCACBC(NatPick,NNatPick,NatPtnm,NatNpt,NatAlign,
     1                  NatCoor,NatCoorAlign,NatMass,iNatConn)
        call findNCACBC(MutPick,NMutPick,MutPtnm,MutNpt,MutAlign,
     1                  MutCoor,MutCoorAlign,MutMass,iMutConn)
       else if (NucleicAcid_yes) then
c --- if NucleicAcid select O4, C1, NX, C2
        call findO4C1NXC2(NatPick,NNatPick,NatPtnm,NatNpt,NatAlign,
     1                  NatCoor,NatCoorAlign,NatMass,iNatConn)
        call findO4C1NXC2(MutPick,NMutPick,MutPtnm,MutNpt,MutAlign,
     1                  MutCoor,MutCoorAlign,MutMass,iMutConn)
       endif

       if (NMutPick.gt.NNatPick) then
        do i=NNatPick+1,NMutPick
         NatCoorAlign(1,i) = 0.d0
         NatCoorAlign(2,i) = 0.d0
         NatCoorAlign(3,i) = 0.d0
         MutMass(i) = 0.d0
        enddo
       endif

c perform alignment
       call rmsd_weight(NMutPick,NatCoorAlign,MutCoorAlign,
     1                  rms,.TRUE.,MutMass)

c prepare the output --- switch back to the native info
       npt = NatNpt
       totmon = NatTotmon
       do i=1,npt
          ptnm(i) = NatPtnm(i)
          poimon(i) = NatPoimon(i)
          coor(1,i) = NatCoor(1,i)
          coor(2,i) = NatCoor(2,i)
          coor(3,i) = NatCoor(3,i)
       enddo
       do i=1,totmon
          moname(i) = NatMoname(i)
          poipt(i) = NatPoipt(i)
       enddo

c copy the coordinates 
c the translation is performed in such a way that the iNatConn and iMutConn particles are the same

       totmon = totmon + 1
       do i=4,NMutPick
        npt = npt + 1
        coor(1,npt)=MutCoorAlign(1,i)+coor(1,iNatConn)-
     1              MutCoorAlign(1,2)
        coor(2,npt)=MutCoorAlign(2,i)+coor(2,iNatConn)-
     1              MutCoorAlign(2,2)
        coor(3,npt)=MutCoorAlign(3,i)+coor(3,iNatConn)-
     1              MutCoorAlign(3,2)
        ptnm(npt)=MutPtnm(MutAlign(i))
        poimon(npt)=totmon
        moname(totmon)=MutMonomerName
       enddo

c now create the crd file with the mutant
       call putcrd(uwcrd,'CHARM')

c now create the poly file
      write(uwpol,'(3a,i6.6)') 'MOLC=(',MutName,') #mon=',totmon
      if (totmon.gt.10) then
      do  i=1,totmon,10
         if ((totmon-i).ge.10) then
         write(uwpol,101) (moname(j),j=i,i+9)
         else
         write(uwpol,101) (moname(j),j=i,totmon)
         endif
      enddo
      else
         write(uwpol,101) (moname(j),j=1,totmon)
      endif
      write(uwpol,102)
101   format(10(a4,1x))
102   format('*EOD')

       stop
       end

       subroutine findNCACBC(ipick,npick,ptnm,npt,align,
     1                       Coor,CoorAlign,mass,iCa)
       implicit none

       integer i 

       integer npt,iCa,npick
       integer ipick(*),align(*)
       double precision mass(*),Coor(3,*),CoorAlign(3,*)
       character*4 ptnm(*)

       do i=1,npt
        mass(i) = 0.d0
        align(i) = 0
       enddo

       do i=1,npt
        if (ipick(i).eq.1) then
         if (ptnm(i).eq.'N') then
          CoorAlign(1,1) = Coor(1,i)
          CoorAlign(2,1) = Coor(2,i)
          CoorAlign(3,1) = Coor(3,i)
          mass(1) = 1.d0
          ipick(i) = 0
          align(1) = i
         else if (ptnm(i).eq.'CA') then
          CoorAlign(1,2) = Coor(1,i)
          CoorAlign(2,2) = Coor(2,i)
          CoorAlign(3,2) = Coor(3,i)
          mass(2) = 1.d0
          ipick(i) = 0
          iCa = i
          align(2) = i
         else if (ptnm(i).eq.'C') then
          CoorAlign(1,3) = Coor(1,i)
          CoorAlign(2,3) = Coor(2,i)
          CoorAlign(3,3) = Coor(3,i)
          mass(3) = 1.d0
          ipick(i) = 0
          align(3) = i
         else if (ptnm(i).eq.'CB') then
          CoorAlign(1,4) = Coor(1,i)
          CoorAlign(2,4) = Coor(2,i)
          CoorAlign(3,4) = Coor(3,i)
          mass(4) = 1.d0
          ipick(i) = 0
          align(4) = i
         endif
        endif
       enddo

       npick = 4

       do i=1,npt
        if (ipick(i).eq.1) then
c --- here we want only the sidechains
         if (ptnm(i).ne.'H' .and. ptnm(i).ne.'HA' .and.
     1       ptnm(i).ne.'O') then
          npick = npick + 1
          CoorAlign(1,npick) = Coor(1,i)
          CoorAlign(2,npick) = Coor(2,i)
          CoorAlign(3,npick) = Coor(3,i)
          mass(npick) = 0.d0
          align(npick) = i
         endif
        endif
       enddo

       if      (mass(1).lt.0.9d0) then
        write (*,*) 'N not found!'
        stop
       else if (mass(2).lt.0.9d0) then
        write (*,*) 'CA not found!'
        stop
       else if (mass(3).lt.0.9d0) then
        write (*,*) 'C not found!'
        stop
       else if (mass(4).lt.0.9d0) then
        write (*,*) 'CB not found!'
        stop
       endif

       return
       end

       subroutine findO4C1NXC2(ipick,npick,ptnm,npt,align,
     1                       Coor,CoorAlign,mass,iC1)
       implicit none

       integer i

       integer npt,iC1,npick
       integer ipick(*),align(*),NHeavy
       double precision mass(*),Coor(3,*),CoorAlign(3,*)
       character*4 ptnm(*),NX
       logical check

       do i=1,npt
        mass(i) = 0.d0
        align(i) = 0
       enddo

       NHeavy = 0
       do i=1,npt
        if (ipick(i).eq.1) then
         check = .false.
         if (ptnm(i)(1:1).eq.'O' .and.
     1       (ptnm(i)(3:3).eq.'A'.or.
     2        ptnm(i)(3:3).eq.'B'.or.
     3        ptnm(i)(3:3).eq.'G')) then
           check = .true.
         endif
         if (ptnm(i)(1:1).eq.'C' .or. 
     1       ptnm(i)(1:1).eq.'N' .or.
     2       ptnm(i)(1:1).eq.'O') then
          if (ptnm(i)(1:1).ne.'P' .and.
     1       ptnm(i)(2:2).ne.'P' .and.
     2       ptnm(i)(3:3).ne."'") then
             if (.not.check) then
             write (*,*) ptnm(i),ptnm(i)(1:1),ptnm(i)(2:2),ptnm(i)(3:3)
             NHeavy = NHeavy + 1
             endif
          endif
         endif
        endif
       enddo

       if (NHeavy.le.9) then 
        NX="N1  "
        write (*,*) 'PYRIMIDINE',NHeavy
       else
        NX="N9  "
        write (*,*) 'PURINE',NHeavy
       endif

       do i=1,npt
        if (ipick(i).eq.1) then
         if (ptnm(i).eq."O4'") then
          CoorAlign(1,1) = Coor(1,i)
          CoorAlign(2,1) = Coor(2,i)
          CoorAlign(3,1) = Coor(3,i)
          mass(1) = 1.d0
          ipick(i) = 0
          align(1) = i
         else if (ptnm(i).eq."C1'") then
          CoorAlign(1,2) = Coor(1,i)
          CoorAlign(2,2) = Coor(2,i)
          CoorAlign(3,2) = Coor(3,i)
          mass(2) = 1.d0
          ipick(i) = 0
          iC1 = i
          align(2) = i
         else if (ptnm(i).eq."C2'") then
          CoorAlign(1,3) = Coor(1,i)
          CoorAlign(2,3) = Coor(2,i)
          CoorAlign(3,3) = Coor(3,i)
          mass(3) = 1.d0
          ipick(i) = 0
          align(3) = i
         else if (ptnm(i).eq.NX) then
          CoorAlign(1,4) = Coor(1,i)
          CoorAlign(2,4) = Coor(2,i)
          CoorAlign(3,4) = Coor(3,i)
          mass(4) = 1.d0
          ipick(i) = 0
          align(4) = i
         endif
        endif
       enddo

       npick = 4

       do i=1,npt
        if (ipick(i).eq.1) then
c --- here we want only the base
         check = .false.
         if (ptnm(i)(1:1).eq.'O' .and.
     1       (ptnm(i)(3:3).eq.'A'.or.
     2        ptnm(i)(3:3).eq.'B'.or.
     3        ptnm(i)(3:3).eq.'G')) then
           check = .true.
         endif
         if (ptnm(i)(1:1).ne.'P' .and. 
     1       ptnm(i)(2:2).ne.'P' .and.
     2       ptnm(i)(3:3).ne."'") then
          if (.not.check) then
          npick = npick + 1
          CoorAlign(1,npick) = Coor(1,i)
          CoorAlign(2,npick) = Coor(2,i)
          CoorAlign(3,npick) = Coor(3,i)
          mass(npick) = 0.d0
          align(npick) = i
          endif
         endif
        endif
       enddo

       if      (mass(1).lt.0.9d0) then
        write (*,*) 'O4 not found!'
        stop
       else if (mass(2).lt.0.9d0) then
        write (*,*) 'C1 not found!'
        stop
       else if (mass(3).lt.0.9d0) then
        write (*,*) 'C2 not found!'
        stop
       else if (mass(4).lt.0.9d0) then
        write (*,*) 'NX not found!'
        stop
       endif

       return
       end
