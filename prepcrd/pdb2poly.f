      Program PDB2puth
c
c     make some editorial changes to pdb file to use it for moil
c     prepare poly.xxx file
c
c     this is crude version tested only for myoglobine - adjust it to your
c     needs: any specific monomer that is supposed to be kept in the simula-
c     tions has to be added to sym; terminus types have to be ad-
c     justed; rules for choosing alternative side chains can be also
c     added;
c     to run the program one has to remove output.pdb and poly.xxx files first,
c     as well as make a copy of the original pdb file called input.pdb
c     JM 01.1997
c
      implicit none

      include 'COMMON/LENGTH.BLOCK'
c      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
c      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c
      character*4 styl1,styl2,mtyp
      character*8 name
      integer namel
      integer urcrd,uwcrd,uwpol
c
      integer of
      logical find
c
      character*6 getchar
      integer level,lc
      integer lpst,lpend
      character*5 ctyp
      data urcrd,uwcrd,uwpol/3*99/

      integer maxblock,nsym,nsymb,maxp,maxlines
      parameter(maxblock=10,nsym=10,nsymb=5,maxp=50000)
      double precision xx,yy,zz,cg,cg1,cg2
      double precision xx1,xx2,yy1,yy2,zz1,zz2
      integer i,n,j,pointline(2*maxblock),numwats,numlines
      integer nmon1,nmon2,na1,nblock,icount,iblock,k
      integer icter,ipoly
      character*1 db
      character*6 sym(nsym),ban(nsymb)
      character*4 poly(maxp),prvelem
      character atom*6,wat*3,elem*4,mono*4,prvmono*4,last*14,long*80
      logical flag,okey,nolig,colig,ok

      data sym(1),sym(2),sym(3),sym(4),sym(5),sym(6),sym(7),sym(8),
     &     sym(9),sym(10)/
     &     'ATOM  ','HETATM','NEXT  ','NEXT  ','NEXT  ','NEXT  '
     &     ,'NEXT  ','NEXT  ','NEXT  ','NEXT  '/
      data ban(1),ban(2),ban(3),ban(4),ban(5)/
     &     'CONECT','TER   ','      ','MASTER','END   '/
      data db /' '/


      open(unit=10,file='inter',status='unknown')
c
c
c     General initialization
c
      lc     = 5
      stdi   = 5
      stdo   = 6
      stderr = 0
      debug  = .false.
      name   = 'pdb2puth'
      namel  = 8

      jnkf = 25
      open (unit=jnkf,status='scratch')

 1    continue

      call rline(name,namel,stdi)
      if (find('acti')) goto 2
c
      if (find('debu')) debug = .true.
      mtyp = getchar('MOLC','XXXX',lc)
      if (find('file')) then
	 if (find('wpol')) then
            uwpol = of()
	 end if
	 if (find('rcrd')) then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','pdb',lc)
            urcrd = of()
	 end if
	 if (find('wcrd'))then
            if (find('bina')) then
               level = 1
               call alert(name,namel,
     $              'Binary crd not supported',24,level)
            end if
            ctyp = getchar('ctyp','pdb',lc)
            uwcrd = of()
	 end if
      end if
      goto 1

 2    continue

      do 15 i=1,maxp
         poly(i)=ban(3)
 15   continue

c     remove blank lines (and other undesired) from pdb file

      maxlines=0
 25   continue
      ok=.false.
      read(urcrd,411,end=26)long
      atom=long(1:6)
      if ((atom.eq.sym(1)).or.(atom.eq.sym(2)).or.
     $     (atom.eq.ban(2))) then
         ok=.true.
      endif
      if (ok) then
         write(10,411) long
         maxlines=maxlines+1
      endif
      goto 25
 26   continue
      close(urcrd)
      rewind 10

c     read pdb file to construct pointer pointline
c     block of lines starting with allowed symbols are pointed by pointline
      do 125 i=1,2*maxblock
         pointline(i)=0
 125  continue
      numlines=0
      nblock=0
      iblock=0
      icter=0
      icount=-1
      flag=.false.
      prvelem='    '
      seq_idx = 0
      prev_res = 'XXX'
      prev_idx = 0
      do 225 i=1,maxlines
       read(10,500)head,atom_id,atom_name,res_name,chain_nm,res_id,x,y,z
500    format(a6,i5,a4,1x,a3,1x,a1,i4,1x,3f8.3)
       do i=1,20
       if (res_name.eq.aminoacid(i)) then
        if (prev_res.eq.'XXX') then
         seq_idx = seq_idx + 1
	 sequence(seq_idx) = NTER
         prev_res = aminoacid(i)
         prev_idx = prev_idx + 1
        end if
        seq_idx = seq_idx + 1
	sequence(seq_idx) = res_name
	if (atom_name.eq.'CA  ') then
          if 
       end if
         read(10,412,end=226)atom,prvmono,elem,mono
         if (icter.eq.0) then
            if (atom.eq.ban(2)) icter=i
            if ((prvelem(1:2).eq.' O').and.(mono(1:3).eq.'HEM')) icter=i
         end if
         prvelem=elem
         if (flag) icount=i-1
         flag=.false.
         do 221 j=1,nsym
            if (atom.eq.sym(j)) then
               flag=.true.
               go to 222
            end if
 221     continue
 222     continue
         if (flag) then
            if (icount+1.ne.i) then
               nblock=nblock+1
               if (nblock.gt.maxblock) then
                  write(stdo,*)' MAXBLOCK too small '
                  stop
               end if
               pointline(nblock)=i
               if (nblock.gt.1) pointline(maxblock+nblock-1)=icount
            end if
         end if
         numlines=numlines+1
 225  continue
c
 226  rewind(10)
c
      pointline(maxblock+nblock)=icount
c     write the first section containing various comments
      iblock=pointline(1)
      if (iblock.gt.1) then
         do 325 i=1,iblock-1
            read(10,411)long
            write(uwcrd,411)long
 325     continue
      end if

c     write N-terminus - change its type if needed
      atom='ATOM '
      elem=' HX2'
      mono='NTER'
      na1=1
      nmon1=0
      xx=9999.999
      yy=xx
      zz=xx
      cg1=1.00
      cg2=00.00
c      write(uwcrd,422)atom,na1,elem,db,mono,nmon1,xx,yy,zz,cg1,cg2
      ipoly=1
      poly(ipoly)=mono
      prvmono=mono
      nmon2=nmon1

c     write some information onto screen
      write(stdo,*) 'no of blocks    no of lines  ',nblock,numlines
      write(stdo,*) 'good lines blocks pointers - start  end '
      do 335 i=1,nblock
         write(stdo,*) pointline(i),pointline(maxblock+i)
 335  continue
c      write(stdo,*) 'N-terminus has been added - check its type'
      write(stdo,422)atom,na1,elem,db,mono,nmon1,xx,yy,zz,cg1,cg2
c      if (icter.eq.0) then
c         write(stdo,*)
c     $        'There is problem with C-terminus - add it manually'
c      else
c         write(stdo,*) 'C-terminus will be added - check its type'
c      end if
c     write all the rest excluding lines starting with "banned" symbols
c     skip all hydrogens
      colig=.false.
      nolig=.false.
      k=1
      do 425 i = iblock,maxlines
         j=pointline(maxblock+k)
         if ((k.eq.nblock).and.(j.lt.numlines)) j=numlines
         if (i.le.j) then
            read(10,419,end=430)atom,na1,elem,db,mono,nmon1,xx,yy,zz,
     &           cg1,cg2,last
            okey=.true.
            if ((elem(1:2).eq.' H').or.(elem(1:2).eq.' D')) go to 425
            if ((elem(1:2).eq.'1H').or.(elem(1:2).eq.'2H')) go to 425
            if (elem(1:2).eq.'3H') go to 425
c
c     remove double side chains
            if ((elem(2:3).eq.prvelem(2:3)).and.(nmon2.eq.nmon1)) then
               if ((db.eq.'B').or.(db.eq.'2')) go to 425
               if ((db.eq.'C').or.(db.eq.'3')) go to 425
               if ((db.eq.'D').or.(db.eq.'4')) go to 425
            end if
            db=' '
c
         else
            okey=.false.
            read(10,411,end=430)long
            if (i+1.eq.pointline(k+1)) k=k+1
         end if
         if (okey) then
c            if (i+1.eq.icter) then
c              if (mono.eq.'GLY ') then
c                  elem=' OX2'
c                  mono='CTRG'
c               else
c                  elem=' OX2'
c                  mono='CTER'
c               end if
c               nmon1=0
c            end if
            if ((mono.eq.'GLN ').or.(mono.eq.'ASN ')) then
               if (elem.eq.' AD1') elem=' OD1'
               if (elem.eq.' AD2') elem=' ND2'
               if (elem.eq.' AE1') elem=' OE1'
               if (elem.eq.' AE2') elem=' NE2'
               if ((elem(2:3).eq.prvelem(2:3))
     $              .and.(nmon2.eq.nmon1)) then
                  if ((db.eq.'B').or.(db.eq.'2')) go to 425
               end if
            end if
            if ((mono.eq.'HOH ').or.(mono.eq.'WAT ').or.
     &           (mono.eq.'WTR ')) then
               if((elem(1:3).eq.' O ').or.
     $              (elem(1:3).eq.' O1')) elem=' OH2'
               mono='TIP3'
               if(elem(1:3).eq.' O2') goto 425
               if(elem(1:3).eq.' O3') goto 425
               if(elem(1:3).eq.' O4') goto 425
            end if
            if (mono(1:1).eq.' ') mono=mono(2:4)//' '
            if ((atom.eq.'HETATM').or.(mono(1:3).eq.'HEM')) then
c     here ligand is not bonded
               if (mono.eq.'SO4 ') mono='SUL '
               if (mono.eq.'HEM ') mono='HEM1'
               if (elem.eq.'FE  ') elem=' FE '
               if (elem.eq.' N A ') elem=' NA '
               if (elem.eq.' N B ') elem=' NB '
               if (elem.eq.' N C ') elem=' NC '
               if (elem.eq.' N D ') elem=' ND '
c     if N or C atom found next, O is assumed to belong to ligand
               if (elem.eq.' N  ') then
                  nolig=.true.
                  mono='NO  '
                  write(stdo,*) 'Ligand found: ',mono,' - check it'
               end if
               if (elem.eq.' C  ') then
                  colig=.true.
                  mono='CO  '
                  write(stdo,*) 'Ligand found: ',mono,' - check it'
               end if
               if ((nolig).and.(mono.eq.'HEM1').and.(elem.eq.' O  '))
     &              mono='NO  '
               if ((colig).and.(mono.eq.'HEM1').and.(elem.eq.' O  '))
     &              mono='CO  '
            end if
c     prepare list of monomers for poly.xxx file
            if (.not.(mono.eq.prvmono)) then
               ipoly=ipoly+1
               poly(ipoly)=mono
            else
               if ((mono.eq.'TIP3').and.(elem.eq.' OH2')) then
                  ipoly=ipoly+1
                  poly(ipoly)=mono
               else
                  if (.not.(nmon1.eq.nmon2)) then
                     ipoly=ipoly+1
                     poly(ipoly)=mono
                  end if
               end if
            end if
            if (ipoly.gt.maxp-9) then
               write(stdo,*) 'MAXP is too small '
               stop
            end if
            prvmono=mono
            prvelem=elem

            write(uwcrd,419)atom,na1,elem,db,mono,nmon1,xx,yy,zz,
     &           cg1,cg2,last
            nmon2=nmon1
         endif


 425  continue
 430  continue

      write(uwcrd,'(a)') "END"
      close(uwcrd)

      write(uwpol,'(3a,i6.6)') 'MOLC=(',mtyp,') #mon=',ipoly
      do 435 i=1,ipoly,10
         write(uwpol,424) (poly(j),j=i,i+9)
 435  continue
      write(uwpol,426)

 410  format(a6)
 411  format(a80)
 412  format(a6,1x,a4,1x,a4,1x,a4)

c     420  format(a6,i5,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,a14)
c     puth	   format(a4,2x,i5,2x,a4,a4,1x,i4,4x,3(f8.4),6x,f6.2)
 419  format(a6,i5,1x,a4,a1,a4,1x,i4,4x,3f8.3,2f6.2,a14)
 420  format(a6,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,a14)
 421  format(a6,i5,1x,a4,1x,a4,2x,i4,3x,3f8.3,2f6.2,a14)
 422  format(a6,i5,1x,a4,a1,a4,1x,i4,4x,3f8.3,2f6.2)
 424  format(10(a4,1x))
 426  format('*EOD')

      write(stdo,'(//10x,a//)') " PDB2PUTH finish editing  !!!!"

      stop
      end
c--------------------------------------------------------
