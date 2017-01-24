        subroutine getcrd(ucrd,type)
c
c subroutine to read coordinates
c currently very simple (almost no checks) version
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer ucrd
        character*5 type
c local
        character*4 char1,char2,char3,char4
        character*1 res_dup,prv_res_dup
        character*6 name
        character*80 zevel
        character*4 aminoacid(20),nucleicacid(8)
        character*8 nucleic_start(8),nucleic_end(8)
        character*1 chain_id,conformer
        character*2 hatom(9)
        integer leng,namel,i,j,imoncmp,imon
        integer istart,i1,j1,k,l,level
        integer nucleic_id
        double precision xtmp,ytmp,ztmp,moretmp
        logical nucleic_true
        data name/'getcrd'/
        data namel/6/
      data aminoacid/'ALA ','ARG ','ASN ','ASP ','CYS ','GLU ','GLN ',
     &     'GLY ','HIS ','ILE ','LEU ','LYS ','MET ','PHE ','PRO ',
     &     'SER ','THR ','TRP ','TYR ','VAL '/
      data nucleicacid/'A   ','G   ','C   ','U   ',
     &      'DA  ','DG  ','DC  ','DT  '/
      data nucleic_start/'A5  ','G5  ','C5  ','U5  ',
     &      'DA5 ','DG5 ','DC5 ','DT5 '/
      data nucleic_end/'A3  ','G3  ','C3  ','U3  ',
     &      'DA3 ','DG3 ','DC3 ','DT3 '/
      data hatom/'1H','2H','3H','4H','5H','6H','7H','8H','9H'/

        rewind ucrd
c initialize all coordinates to 9999.0
        do 1 i=1,npt
         chain(i) = '*'
         do 1 j=1,3
         coor(j,i) = 9999.d0
1       continue
        if (type.eq.'CHARM') then
        write(stdo,*)' Reading CHARMm coordinate file '
c
c               read the title and find out its length
c
        leng = 0
2       continue
        read(ucrd,105)zevel
105     format(80a)
        if (zevel(1:1).eq.'*') then
         do 21 k=80,1,-1
           if (zevel(k:k).ne.' ') then
            call echo(name,namel,zevel,k)
            go to 22
           end if
21       continue
22       continue
         leng = leng + 1
         go to 2
        end if
c
c               check that title exist
c
        if (leng.eq.0) then
         level = 1
         call alert(name,namel,'Crd file with no title',22,level)
         return
        end if
c
c               do not forget to backspace the file !!
c
        backspace ucrd
        read(ucrd,*) i
        if (i.ne.0 .and. i.ne.npt) then
         level = 1
         call alert(name,namel,
     1   'Number of particles do not match',32,level)
         return
        end if
        istart = 1
        if (debug) then
         write(stdo,*)' totmon npt poipt '
     1   ,totmon,npt,(poipt(i),i=1,totmon)
        end if
        do 5 i=1,totmon
         do 4  j=istart,poipt(i)
c read one line of coordinate
          read(ucrd,100,err=7,end=55)i1,j1,char1,char2,xtmp
     1          ,ytmp,ztmp,char3,char4,moretmp
100       format(i7,i7,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)   
c check that this is the right monomers name
          if (char1.ne.moname(i)) then
c         if (char1.ne.moname(i).or.j1.ne.i) then
c
c may be we just miss some atoms (this is ok sometime..
c for example when hydrogens are not included
c
           if (j1.ne.i) then
            level = 0
            write(stdo,101)i,moname(i)
101         format(1x,' Missing pts for monomer ',i7,' type = ',
     1a4)
            call alert(name,namel,'Missing particles',17,level)
            backspace ucrd
            go to 4
           else if (j1.eq.totmon) then
            level = 0
            write(stdo,102)i,moname(i)
102         format(1x,' Missing pts for last monomer ',i7,' type = ',
     1a4)
            call alert(name,namel,'Missing particles',17,level)
            backspace ucrd
            go to 5
           else
            level = 1
            write(stdo,103)i,char1,moname(i)
103         format(1x,'mono ',i7,' Read = ',a4,' Expected = ',a4)
            call alert(name,namel,'Monomers do not match',21,level)
           end if
          else
            do 3 k=istart,poipt(i)
             if (char2 .eq.  ptnm(k))then 
              if (cplbl(k).eq.0) then
               coor(1,k) = xtmp
               coor(2,k) = ytmp
               coor(3,k) = ztmp
               more(k) = moretmp
               go to 35
              else 
               l = 0
27             continue
               if (coor(1,k+l).gt.9998.d0 .and. char2.eq.ptnm(k)) then
                coor(1,k+l) = xtmp
                coor(2,k+l) = ytmp
                coor(3,k+l) = ztmp
                more(k+l)  = moretmp
                go to 35
               else if (cplbl(k+l+1).gt.cplbl(k+l)) then
                l = l +1
                go to 27
               end if
              end if
             end if
3           continue
35          continue
          end if
4         continue
          istart=poipt(i) + 1
5        continue
c check that there are not any unidentified atoms
55       continue
         do 6 i=1,npt
          if (coor(1,i).gt.9998.0) then
           if(ptnm(i)(1:1).ne.'H') then
           write(stdo,104)i,ptnm(i),poimon(i),moname(poimon(i))
104        format(1x,' Coordinates of particle ',i7,2x,a4,
     1' not defined in monomer ',i7,2x,a4)
            end if
c@           level = 0
c@           call alert(name,namel,'Undefined coordinates',21,level)
          end if
6        continue
         return
        else if (type(1:3).eq.'pdb') then
        write(stdo,*)' Reading pdb coordinate file '
c
c check that the current line is not ATOM line. If not simply echo the
c line and keep reading. Stop this loop when ATOM is detected backspace
c the file and start read thee coordinates
c
         imoncmp = 0
         imon    = 0
         istart = 1
         prv_res_dup = ' '
         do 11 i=1,totmon

c
c first check. If the current monomer we attempt to read is N-terminal skip
c the process and jump to next monomers (we build the hydrogens).
c
           if (moname(i).eq.'NTER' .or. moname(i).eq.'NTRG'
     1      .or. moname(i).eq.'NTRP')  then
               imon = imon + 1
                do j=poipt(i-1)+1,poipt(i)
                k = poipt(j)
                more(k)   = 0.
                chain(k)  = '*'
                go to 125
               end do
           end if
check that the present line is ATOM or HETA
          do 12 j=istart,poipt(i)
8          continue
           read(ucrd,105,end=75,err=7)zevel
           if (zevel(1:4).ne.'ATOM' .and. zevel(1:4).ne.'HETA') goto 8
           backspace ucrd
c read one line of coordinates
           read(ucrd,106,end=75,err=7)
     1     zevel(1:4),j1,char2,conformer,char1,chain_id,i1,res_dup,
     2     xtmp,ytmp,ztmp,moretmp
106        format(a4,2x,i5,1x,a4,a1,a4,a1,i4,a1,3x,3(f8.4),6x,f6.2)
c
c sometimes the sequence starts with negative amino acid (to include a leading peptide)
c
           if (j1.lt.0) then
            imon = j1 -1
            imoncmp = j1 -1
           end if
c
c The residue name is in char1. There are 4 charcters for residue name.
c make sure that there are no leading spaces
c

c
c shift the residue & atom names to the left (avoid spaces on the left)
c
       do k=1,4
       if (char1(1:1).eq.' ') then
        do l=1,3
         char1(l:l)=char1(l+1:l+1)
        end do
        char1(4:4)=' '
       end if
       if(char2(1:1).eq.' ') then
        do l=1,3
         char2(l:l) = char2(l+1:l+1)
        end do
        char2(4:4) = ' '
       end if
       end do
c
c check if current residue is a nucleic acid if yes get nucleic_id
c and set nucleic_true to ture
c
       nucleic_true = .false.
       do k=1,8
        if (char1.eq.nucleicacid(k)) then
         nucleic_true = .true.
         nucleic_id = k
         go to 14  
        end if
       end do
14     continue
c
c check for side chain rotamers. If more than one exists keep only
c  rotamer with the index "1"
c
	   if (.not.(conformer.eq.' ' .or. conformer.eq.'1'
     1         .or. conformer.eq.'A')) go to 8

c
c check for water molcules and repalce all PDB naming by TIP3
c
           if (char1.eq.'HOH ' .or. char1.eq.'WAT ' .or.
     1       char1.eq.'WTR ') then
             char1 = 'TIP3'
             if (char2(1:1).eq.'O') then
              char2 = 'OH2 '
             end if
           end if
c
c Hydrogens in the PDB are more pain than gain. Here we just skip them and we will build
c them back by puth. This code needs to changed if the original hydorgen positions
c are needed.
c
	 if (char2(1:1).eq.'H') then
c
c Hg is silver. A typical heavy atoms used the crystallography, that we keep.
c
          if (char1.ne.'HG') go to 8
         else
          do k=1,9
           if (char2(1:2).eq.hatom(k)) go to 8
          end do
         end if

c
c check for HIS / HIP definition in the wcon
c           write(*,*)'char2 char1 mononame = ',char2,' ',char1,
c     1  ' ',moname(i)
           if (moname(i).eq.'HIP ' .and. char1.eq.'HIS ') char1='HIP '
c           write(*,*)' char1 moname(i) AFTER = ',char1,' ',moname(i)

c
c check for C-terminal
c
           if (moname(i).eq.'CTER' .or. moname(i).eq.'CTRG'
     1          .or. moname(i).eq.'CTRP') then
             if (char2.eq.'OXT') then
              k = poipt(i)
              coor(1,k) = xtmp
              coor(2,k) = ytmp
              coor(3,k) = ztmp
              more(k)   = moretmp
              chain(k)  = chain_id
              imon = imon + 1
             else
              backspace ucrd
             end if
             go to 125
           end if
           if (char2.eq.'OXT') then
            if (moname(i).eq.'CTER' .or. moname(i).eq.'CTRG'
     1       .or. moname(i).eq.'CTRP') then    
                k = poipt(i)
                coor(1,k) = xtmp
                coor(2,k) = ytmp
                coor(3,k) = ztmp
                more(k)   = moretmp
                chain(k)  = chain_id
                imon = imon + 1
                go to 125
            else if (moname(i+1).eq.'CTER'
     1       .or. moname(i+1).eq.'CTRG'
     2       .or. moname(i+1).eq.'CTRP') then    
                backspace ucrd
                go to 125
            end if
           end if

c
c check if new residue res_dup is for sequences with numbers
c that looks like 166 166A 166B 167 etc.
c
c
           if (imoncmp.ne.i1 .or. res_dup.ne.prv_res_dup) then
            imon    = imon + 1
            imoncmp = i1
            prv_res_dup = res_dup
            if (j.ne.istart) then
             backspace ucrd
             go to 125
            end if
           end if
           if (char1.ne.moname(i)) then

           if (nucleic_true) then
              if (moname(i).eq.nucleicacid(nucleic_id) .or.
     1           moname(i).eq.nucleic_start(nucleic_id) .or.
     2           moname(i).eq.nucleic_end(nucleic_id)) then
                do k=istart,poipt(i)
                 if (char2.eq.ptnm(k)) then
                  coor(1,k) = xtmp
                  coor(2,k) = ytmp
                  coor(3,k) = ztmp
                  more(k) = moretmp
                  chain(k)  = chain_id
                 end if
                end do
               else
                 level = 1
                 write(stdo,103)i,char1,moname(i)
                 call alert(name,namel,'Monomers do not match',21,level)
               end if
c
c check for change of chains
c
            if (char2.eq.'OXT') then
                k = poipt(i)
                coor(1,k) = xtmp
                coor(2,k) = ytmp
                coor(3,k) = ztmp
                more(k)   = moretmp
                chain(k)  = chain_id
                go to 125
            end if
            else if (imon.ne.i .and. (.not.nucleic_true)) then 
c --		 char1.ne.moname
c
             level = 0
             write(stdo,101)i,moname(i)
             call alert(name,namel,'Missing particles',17,level)
             backspace ucrd
             go to 12
            else if (imon.eq.totmon) then
c --             char1.ne.moname
c
             level = 0
             write(stdo,102)i,moname(i)
             call alert(name,namel,'Missing particles',17,level)
             backspace ucrd
             go to 12
            end if
           else if (char1.eq.'HIS' .and. moname(i).eq.'HIP')  then
            do k=istart,poipt(i)
             if (char2.eq.ptnm(k)) then
              coor(1,k) = xtmp
              coor(2,k) = ytmp
              coor(3,k) = ztmp
              chain(k)  = chain_id
             end if
            end do
           else
c --            char1.ne.moname ( so here it is equal)
c
            do 10 k=istart,poipt(i)
             if (char2.eq.ptnm(k)) then
              coor(1,k) = xtmp
              coor(2,k) = ytmp
              coor(3,k) = ztmp
              more(k) = moretmp
              chain(k)  = chain_id
             end if
10          continue
           end if
12        continue
125       continue
          istart = poipt(i) + 1
11       continue
c check that ther are not any unidentified particles
c Hydrogen are trivial missing particles and are not reported
c
75      continue
         do 13 i=1,npt
          if (coor(1,i).gt.9998.0d0) then
           if (.not.(ptnm(i)(1:1).eq.'H' .or. ptnm(i)(1:2).eq.'1H'
     1       .or. ptnm(i)(1:2).eq.'2H')) then
            write(stdo,104)i,ptnm(i),poimon(i),moname(poimon(i))
           end if
          end if
13       continue
         return
        end if
7       continue
        level = 0
        call alert(name,namel,'Error during read',17,level)
        return
        end
