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

      integer maxlines
      parameter(maxlines=100000)
      character*3 aminoacid(20)
      character*6 head
      character*4 atom_name,res_name,prev_res
      character*1 chain_nm,chain_nm_prev
      character*4 sequence(maxlines)

      integer atom_id,res_id,prev_idx
      integer i,j,seq_idx
      real x,y,z,prev_x,prev_y,prev_z,dist

      logical aminoacid_ok

      data aminoacid/'ALA','ARG','ASN','ASP','CYS','GLU','GLN',
     &     'GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO',
     &     'SER','THR','TRP','TYR','VAL'/
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
	write(*,*)" This program generate a sequence (poly) file for"
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

      prev_res = 'XXXX'
      prev_idx = 1
      res_id = 0
      chain_nm_prev = 'X'
      prev_x = 1.e7
      prev_y = 1.e7
      prev_z = 1.e7
      do 225 i=1,maxlines
       read(urcrd,499,end=226)head
499    format(a6)
       if (head.eq.'ENDMDL') go to 226
       if (head.eq.'ATOM  ' .or. head.eq.'HETATM') then
       backspace(urcrd)
       read(urcrd,500)head,atom_id,atom_name,res_name,
     &  chain_nm,res_id,x,y,z
500    format(a6,i5,a4,1x,a4,a1,i4,1x,3f8.3)
      
       if (head.eq.'TER ') go to 224 
       if (res_id.ne.prev_idx) then
c
c if the residue index is different from previous read this is a new residue
c
       aminoacid_ok = .false.
       do j=1,20
        if (res_name.eq.aminoacid(j)) then
		aminoacid_ok = .true.
		go to 21
	end if
       end do
21	continue
c
c if the current residue is an amino acid, make a few checks.
c check if it is the firt residue (requires adding n terminal)
c check if the chain is broken (distance between C-alpha too large
c check if this is a last resiude in a chain and add C terminal
c
	if (aminoacid_ok) then
         if (prev_res.eq.'XXXX' .or. prev_res.eq.'CTER'
     &  .or. prev_res.eq.'CTRG' .or. prev_res.eq.'CTRP') then
          seq_idx = seq_idx + 1
          if (res_name.eq.'GLY ') then
           sequence(seq_idx) = 'NTRG'
          else if (res_name.eq.'NTRP') then
	   sequence(seq_idx) = 'NTRP'
          else
           sequence(seq_idx) = 'NTER'
          end if
         end if
         if (res_id.lt.prev_idx .or. (chain_nm.ne.chain_nm_prev
     &        .and.seq_idx.gt.1)) then
           seq_idx = seq_idx + 1
           if (prev_res.eq.'GLY') then
             sequence(seq_idx) = 'CTRG'
           else if (prev_res.eq.'PRO') then
             sequence(seq_idx) = 'CTRP'
           else
             sequence(seq_idx) = 'CTER'
           end if
          seq_idx = seq_idx + 1
          sequence(seq_idx) = res_name
          prev_res = sequence(seq_idx)
          prev_idx = res_id
          chain_nm_prev = chain_nm
         end if
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
         end if
	else
	 if (res_name.eq.'HOH '.or.res_name.eq.'WAT '.or.
     &		res_name.eq.'WTR ') then
		seq_idx = seq_idx + 1
		sequence(seq_idx) = 'TIP3'
          end if
         end if
        end if
       end if
224    continue
225    continue
226    continue

      write(uwpol,'(3a,i6.6)') 'MOLC=(',mtyp,') #mon=',seq_idx
      do 435 i=1,seq_idx,10
         write(uwpol,424) (sequence(j),j=i,i+9)
424   format(10(a4,1x))
 435  continue
      write(uwpol,426)
 426  format('*EOD')

      stop
      end
c--------------------------------------------------------
