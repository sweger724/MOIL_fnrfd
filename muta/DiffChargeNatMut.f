C this code reads the native wcon and the mutant wcon
C and returns a set of atoms whose charge has changed
C 
C this information is then used in muta to perform PME
C on uncharged systems

        program DiffChargeNatMut

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/LINE.BLOCK'

        integer namel
        character*16 name

        logical find
        logical Nmuta,Mmuta

c loops
        integer i,j,istart,jstart,iend,jend
c I/O
        integer ipick(maxpt),npick
        logical pickpt
        integer of,uncon,umcon,uwdif

c special variables for comparison
        integer Nnpt,Mnpt,Nnmono,Mnmono
        integer Npoimon(maxpt),Mpoimon(maxpt)
        integer Nmutaid(maxpt),Mmutaid(maxpt)
        integer im,nm,MutantMonomer(maxmono)
        double precision Nptchg(maxpt),Mptchg(maxpt)
        character*4 Nptnm(maxpt),Mptnm(maxpt)
        character*4 Nmoname(maxmono),Mmoname(maxmono)
c extra set of atoms
        integer Nextra(maxpt),Mextra(maxpt)
        stdi = 5
        stdo = 6
        stderr = 0  
 
        pickpt = .false.

        uncon = -9999
        umcon = -9999
        uwdif = -9999

        name  = 'DiffChargeNatMut'
        namel = 16

        jnkf = 25
        open (unit=jnkf,status='scratch')

        Nmuta = .false.
        Mmuta = .false.
        muta  = .false.

1       continue
        call rline(name,namel,stdi)
        if (find('file')) then
c native conn file
         if (find('ncon')) then
          uncon  = of()
          if (find('muta')) Nmuta = .true.
          if (Nmuta) muta = .true.
          call rconn(uncon) 
          muta = .false.
          Nnpt = npt
          Nnmono = totmon
          do i=1,Nnpt
           Nptnm(i) = ptnm(i)
           Npoimon(i) = poimon(i)
           Nptchg(i) = ptchg(i)
           Nmutaid(i) = mutaid(i)
           Nextra(i) = 0
          enddo
          do i=1,Nnmono
           Nmoname(i) = moname(i)
          enddo
         endif
c mutant conn file
         if (find('mcon')) then
          umcon = of()
          if (find('muta')) Mmuta = .true.
          if (Mmuta) muta = .true.
          call rconn(umcon)
          muta = .false.
          Mnpt = npt
          Mnmono = totmon
          do i=1,Mnpt
           Mptnm(i) = ptnm(i)
           Mpoimon(i) = poimon(i)
           Mptchg(i) = ptchg(i)
           Mmutaid(i) = mutaid(i)
           Mextra(i) = 0
          enddo
          do i=1,Mnmono
           Mmoname(i) = moname(i)
          enddo
         endif
c output file
         if (find('wdif')) then
          uwdif = of()
         endif
        endif

C select extra atoms to use to keep PN and PM uncharged
        if (find('Nxtr')) then
         call pick(ipick,npick)
         do i=1,Nnpt
          if (ipick(i).gt.0) write (*,*) 'N picked',i
          Nextra(i) = ipick(i)
         enddo
        endif
        if (find('Mxtr')) then
         call pick(ipick,npick)
         do i=1,Mnpt
          if (ipick(i).gt.0) write (*,*) 'M picked',i
          Mextra(i) = ipick(i)
         enddo
        endif

        if (find('acti')) goto 2
        goto 1

2       continue

        if (uncon .lt. 0) then
         call alert(name,namel,'nat con file not open!',22,1)
        else if (umcon .lt. 0) then
         call alert(name,namel,'mut con file not open!',22,1)
        else if (uwdif .lt. 0) then
         call alert(name,namel,'diff file not found!!!',22,1)
        endif

        if (((.not.Nmuta).or.(.not.Mmuta))) then
         call alert(name,namel,'conn files are not of mutants',29,1)
        endif

        if (Nnmono .ne. Mnmono) then
         call alert(name,namel,'cannot handle differnt nmmono',29,1)
        endif

        nm = 0
C find the mutant monomer
        do i=1,Nnmono
         if (Nmoname(i).ne.Mmoname(i)) then
          nm = nm + 1
          MutantMonomer(nm) = i
          write (*,*) 'Mutant Monomer #',nm,'=',i
         endif
        enddo

        istart=1
        jstart=1
        iend=1
        jend=1

        do im=1,nm
         do i=iend,Nnpt
          if (Npoimon(i).eq.MutantMonomer(im)) goto 3
         enddo
3        continue
         istart = i
         do i=istart,Nnpt
          if (Npoimon(i).gt.MutantMonomer(im)) goto 30
         enddo
30       continue
         iend=i
         do j=jend,Mnpt
          if (Mpoimon(j).eq.MutantMonomer(im)) goto 4
         enddo
4        continue
         jstart = j
         do j=jstart,Mnpt
          if (Mpoimon(j).gt.MutantMonomer(im)) goto 40
         enddo
40       continue
         jend=j

        do 5 i=istart,iend-1
C select all the P particles that belong to the same #mono and
C have the same name ...          
         do j=jstart,jend-1 
          if ( (Npoimon(i).eq.Mpoimon(j)) .and.
     1         (Nptnm(i)(1:4).eq.Mptnm(j)(1:4))     .and.
     2         (Nmutaid(i).eq.0)           .and.
     3         (Mmutaid(j).eq.0)          ) goto 6
         enddo
         goto 5
6        continue
C ... but different charge OR are explicitly selected
         if ((Nptchg(i).ne.Mptchg(j)).or.(Nextra(i).eq.1).or.
     1                                   (Mextra(j).eq.1)) then
          write (uwdif,'(2i7,1x,1a4,1x,1a4,2f10.5)')
     1    i,Npoimon(i),Nmoname(Npoimon(i)),Nptnm(i),Nptchg(i),Mptchg(j)
         endif
5       continue

        enddo

        write (uwdif,*) 'END'

        stop
        end
