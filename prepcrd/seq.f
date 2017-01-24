      program pdbseq
c
c extract sequence from a PDB file and generate a poly file
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
      data urcrd,uwpol/2*99/

      integer maxlines,no_aa,no_na
      parameter(maxlines=100000)
      parameter(no_aa=22)
      parameter(no_na=8)
      character*4 aminoacid(no_aa)
      character*4 nucleicacid(no_na)
      character*4 nucleic_start(no_na)
      character*4 nucleic_end(no_na)
      character*6 head
      character*4 atom_name,res_name,prev_res
      character*1 chain_nm,chain_nm_prev,res_dup
      character*1 prev_res_dup
      character*1 chain_id(maxlines+1)
      character*4 sequence(0:maxlines)

      integer atom_id,res_id,prev_idx
      integer i,j,k,seq_idx
      real x,y,z,prev_x,prev_y,prev_z,dist

      logical aminoacid_ok,amino_prev_ok
      logical nucleic_ok,nucleic_prev_ok
      logical test,lTER

      data aminoacid/'ALA ','ARG ','ASN ','ASP ','CYS ','GLU ','GLN ',
     &     'GLY ','HIS ','ILE ','LEU ','LYS ','MET ','PHE ','PRO ',
     &     'SER ','THR ','TRP ','TYR ','VAL ','MSE ','CMT'/
      data nucleicacid/'A   ','G   ','C   ','U   ',
     &      'DA  ','DG  ','DC  ','DT  '/
      data nucleic_start/'A5  ','G5  ','C5  ','U5  ',
     &      'DA5 ','DG5 ','DC5 ','DT5 '/
      data nucleic_end/'A3  ','G3  ','C3  ','U3  ',
     &      'DA3 ','DG3 ','DC3 ','DT3 '/
c
c
c     General initialization
c
      lc     = 4
      stdi   = 5
      stdo   = 6
      stderr = 0
      debug  = .false.
      name   = 'pdbseq'
      namel  = 6
	write(*,*)" This program generates a sequence (poly) file "
	write(*,*)" the conn program. The information is extracted "
	write(*,*)" from the ATOM record (not SEQRES) to avoid missing "
        write(*,*)" coordinates. If the pdb entry is from NMR "
	write(*,*)" this program will terminate after the first model "
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
            ctyp = getchar('ctyp','pdb',lc)
            urcrd = of()
	 end if
      end if
      goto 1

 2    continue
 25   continue

      chain_id(maxlines + 1)= 'X'
      sequence(0)(1:4)= 'XXXX'
c
c make a "quick" read to figure out where to start counting amino acids
c
3     continue
      read(urcrd,499)head
      if (head.eq.'ATOM  ' .or. head.eq.'HETATM') then
       backspace(urcrd)
       read(urcrd,500)head,atom_id,atom_name,res_name,
     &  chain_nm,res_id,res_dup,x,y,z
       backspace(urcrd)
      else
       go to 3
      end if
      prev_res = 'XXXX'
      seq_idx = 0
      prev_idx = res_id - 1
      prev_res_dup = ' '
      chain_nm_prev = 'X'
      prev_x = 1.e7
      prev_y = 1.e7
      prev_z = 1.e7
      prev_res_dup = ' '
      amino_prev_ok = .false.
      aminoacid_ok = .false.
      do 225 i=1,maxlines
       read(urcrd,499,end=226)head
499    format(a6)
       if (head.eq.'ENDMDL') go to 226
       if (head.eq.'ATOM  ' .or. head.eq.'HETATM') then
       backspace(urcrd)
       read(urcrd,500)head,atom_id,atom_name,res_name,
     &  chain_nm,res_id,res_dup,x,y,z
500    format(a6,i5,1x,a4,1x,a4,a1,i4,a1,3x,3f8.3)
      
       do k=1,4
       if (res_name(1:1).eq.' ') then
	do j=1,3
         res_name(j:j)=res_name(j+1:j+1)
        end do
        res_name(4:4)=' '
       end if
       end do

       if (head.eq.'TER '.and. amino_prev_ok) then
c
c add C-terminal
c
            seq_idx = seq_idx + 1
            if (prev_res.eq.'GLY') then
             sequence(seq_idx) = 'CTRG'
            else if (prev_res.eq.'PRO') then
             sequence(seq_idx) = 'CTRP'
            else if (prev_res.eq.'CMT') then
             write(*,*)' CMT includes terminal'
             seq_idx = seq_idx - 1
            else
             sequence(seq_idx) = 'CTER'
            end if
            seq_idx = seq_idx + 1
            amino_prev_ok = .false.
            lTER = .true.
            go to 224 
        else
            lTER = .false.
        end if
       if (res_id.ne.prev_idx .or.
     1     res_dup.ne.prev_res_dup) then
       if (res_id-prev_idx.gt.1) then
        write(*,*)' ** WARNING A LIKELY BREAK AT RESIDUE : ',res_id
       end if
c
c if the residue index is different from previous read this is a new residue
c
        if(res_dup.ne.prev_res_dup) prev_res_dup=res_dup
        amino_prev_ok = aminoacid_ok
        aminoacid_ok = .false.
       do j=1,no_aa
        if (res_name.eq.aminoacid(j)) then
		aminoacid_ok = .true.
		go to 21
	end if
       end do
21     continue
c
c if the current residue is an amino acid, make a few checks.
c check if it is the firt residue (requires adding n terminal)
c check if the chain is broken (distance between C-alpha too large
c check if this is a last resiude in a chain and add C terminal
c
	if (aminoacid_ok) then
c
c XXXX is the first resiude, so if prev_ res is XXXX this is the beginning of the chain
c
         if (prev_res.eq.'XXXX' .or. lTER) then
          seq_idx = seq_idx + 1
          if (res_name.eq.'GLY ') then
           sequence(seq_idx) = 'NTRG'
          else if (res_name.eq.'PRO ') then
	   sequence(seq_idx) = 'NTRP'
          else
           sequence(seq_idx) = 'NTER'
          end if
          seq_idx = seq_idx + 1
          sequence(seq_idx) = res_name
          chain_id(seq_idx) = chain_nm
          prev_res = sequence(seq_idx)
          prev_idx = res_id
          chain_nm_prev = chain_nm
          go to 225
c
c below are indicator for a new chain (a change in the chain indicator)
c
         else if ((chain_nm.ne.chain_nm_prev
     &        .and.seq_idx.gt.2) .and. .not. lTER) then
c
c One chain ended and a new chain started (CTER detected and a new amino acid
c opening a new chain was detected). Add appropriate NTER and current residue
c
           if (amino_prev_ok) then
            seq_idx = seq_idx + 1

            if (prev_res.eq.'GLY') then
             sequence(seq_idx) = 'CTRG'
            else if (prev_res.eq.'PRO') then
             sequence(seq_idx) = 'CTRP'
            else if (prev_res.eq.'CMT') then
c       CMT includes terminal
             seq_idx = seq_idx - 1
            else
             sequence(seq_idx) = 'CTER'
            end if
            seq_idx = seq_idx + 1
            if (res_name.eq.'GLY ') then
             sequence(seq_idx) = 'NTRG'
            else if (res_name.eq.'PRO ') then
	     sequence(seq_idx) = 'NTRP'
            else
             sequence(seq_idx) = 'NTER'
            end if
           else 				! amino_prev_ok is NOT OK
            seq_idx = seq_idx + 1
            if (res_name.eq.'GLY ') then
             sequence(seq_idx) = 'NTRG'
            else if (res_name.eq.'PRO ') then
	     sequence(seq_idx) = 'NTRP'
            else
             sequence(seq_idx) = 'NTER'
            end if
           end if
           end if
           seq_idx = seq_idx + 1
           sequence(seq_idx) = res_name
           chain_id(seq_idx) = chain_nm
           prev_res = sequence(seq_idx)
           prev_idx = res_id
           chain_nm_prev = chain_nm
          
         if (atom_name.eq.'CA  ') then
          if (prev_x .gt.1.e6) then
             prev_x = x
             prev_y = y
             prev_z = z
          else
             dist = 0
	     dist = dist + (prev_x - x)**2
	     dist = dist + (prev_y - y)**2
	     dist = dist + (prev_z - z)**2
             dist = sqrt(dist)
             if (dist.gt.4.0) then
		write(*,*)' ** AMINO ACIDS i,i+1 too separated '
		write(*,*)' ** i i+1 r ',prev_idx,seq_idx,dist
             end if
             prev_x = x
             prev_y = y
             prev_z = z
           end if
          seq_idx = seq_idx + 1
          sequence(seq_idx) = res_name
          prev_res = sequence(seq_idx)
          chain_id(seq_idx) = chain_nm
          prev_idx = res_id
          chain_nm_prev = chain_nm
          prev_res_dup = res_dup
         end if
	else                            ! amino acid (not) OK
         if (amino_prev_ok) then
c
c the last group was an amino acid. The current read group is not OXT
c  A C terminal is required.
c
           
           if (sequence(seq_idx).eq.'GLY') then
             seq_idx = seq_idx + 1
             sequence(seq_idx) = 'CTRG'
           else if (sequence(seq_idx).eq.'PRO') then
             seq_idx = seq_idx + 1
             sequence(seq_idx) = 'CTRP'
           else
             seq_idx = seq_idx + 1
             sequence(seq_idx) = 'CTER'
           end if
           amino_prev_ok = .false.
          end if
c end of if amino_prev_ok
c
		seq_idx = seq_idx + 1
                prev_res = sequence(seq_idx)
                prev_idx = res_id
                chain_nm_prev = chain_nm
	 if (res_name.eq.'HOH '.or.res_name.eq.'WAT '.or.
     &		res_name.eq.'WTR ') then
		sequence(seq_idx) = 'TIP3'
         else if(res_name.eq.'HEM ') then
                sequence(seq_idx) = 'HEME'
         else
                sequence(seq_idx) = res_name
                chain_id(seq_idx) = chain_nm
         end if

         end if
        end if
       end if
224    continue
225    continue
226    continue


c end of read file was reached. check if cter was inserted properly
c
       if (sequence(seq_idx).eq.'CTER' .or. 
     1     sequence(seq_idx).eq.'CTRG' .or.
     2     sequence(seq_idx).eq.'CTRP' .or.
     3     sequence(seq_idx).eq.'NH2 ') go to 1001
        aminoacid_ok = .false.
        do i=1,no_aa
         if (sequence(seq_idx).eq.aminoacid(j)) then
	  aminoacid_ok = .true.
          go to 1002
         end if
        end do
1002	continue
        if (aminoacid_ok) then
         seq_idx = seq_idx + 1
         if (prev_res.eq.'GLY') then
            sequence(seq_idx) = 'CTRG'
         else if (prev_res.eq.'PRO') then
            sequence(seq_idx) = 'CTRP'
         else
            sequence(seq_idx) = 'CTER'
         end if
        end if
1001  continue
c
c DO THE SUBSTITUTIONS OF NUCLEIC ACIDS 5' and 3' HERE.
c
      do i=1,seq_idx
       nucleic_ok = .false.
       do j=1,no_na
        if (sequence(i).eq.nucleicacid(j)) then
                nucleic_ok = .true.
                go to 210
        end if
       end do
210     continue
       if (nucleic_ok) then
        if (i.eq.1) then
         do j=1,4
          if(sequence(i)(j:j).eq.' ') then
             sequence(i)(j:j)='5'
            go to 22
          end if
         end do
        else if (i.eq.seq_idx) then
         do j=1,4
          if(sequence(i)(j:j).eq.' ') then
             sequence(i)(j:j)='3'
            go to 22
          end if
         end do
        end if
        nucleic_prev_ok=.false.
        do j=1,no_na
          if (sequence(i-1).eq.nucleicacid(j)
     1     .or. sequence(i-1).eq.nucleic_start(j)) then
           nucleic_prev_ok = .true.
          end if
        end do
        if (.not.nucleic_prev_ok) then
         do j=1,4
          if (sequence(i)(j:j).eq.' ') then
            sequence(i)(j:j)='5'
            go to 22
          end if
         end do
        end if
        nucleic_prev_ok=.false.
        do j=1,no_na
          if (sequence(i+1).eq.nucleicacid(j)
     1     .or. sequence(i+1).eq.nucleic_start(j)) then
           nucleic_prev_ok=.true.
          end if
         end do
         if (.not.nucleic_prev_ok) then
          do j=1,4
           if (sequence(i)(j:j).eq.' ') then
              sequence(i)(j:j)='3'
              go to 22
           end if
          end do
         end if
         if (i.eq.seq_idx) go to 22
         if (chain_id(i).ne.chain_id(i+1)) then
          nucleic_prev_ok = .false.
          do j=1,no_na
           if (sequence(i+1).eq.nucleicacid(j)) then
            nucleic_prev_ok = .true.
           end if
          end do
          if (nucleic_prev_ok) then
             do j=1,4
              if (sequence(i)(j:j).eq.' ') then
               sequence(i)(j:j)='3'
               sequence(i+1)(j:j)='5'
               go to 22
              end if
             end do
           end if
          end if
         end if
22	continue
         end do

         
      write(uwpol,'(3a,i6.6)') 'MOLC=(',mtyp,') #mon=',seq_idx
      write(uwpol,424) (sequence(i),i=1,seq_idx)
424   format(10(a4,1x))
 435  continue
      write(uwpol,426)
 426  format('*EOD')
      stop
      end
c--------------------------------------------------------
