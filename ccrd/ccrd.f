       program ccrd
c 
c converting coordinate formats
c 
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        character*4 styl1,styl2
        character*5 name
        integer nofreez(maxpt),ipick(maxpt),jpick(maxpt)
        integer nofrez1(maxpt),ipic1(maxpt),jpic1(maxpt)
        integer namel,inofrz,npick,npic1,k
        double precision rms, tmp_mass(maxpt)
        integer ucon,urcrd,uwcrd,ucmbn,urpth,uwpth,tobin
        integer geti,of
        logical find,fopen,lcmbn,wsub,ovlp,bConv
        integer i,ich,level
        integer lpend,ncmbn,ntotal

        integer l,j
c rbin = 1 if reading binary dcd, 0 otherwise
c wbin = 1 if writeing binary dcd, 0 otherwise
        integer rbin,wbin
        
        data ucon,urcrd,uwcrd/3*99/
c
c General initialization
c
        stdi   = 5
        stdo   = 6
        stderr = 0
        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        nbulk  = 0
        inofrz = 0
        debug  = .false.
        norew  = .false. 
        name   = 'ccrd'
        namel  = 4
        wsub = .false.
        ovlp = .false.
        bConv = .false.
        npic1 = 0
        urpth = 0
        uwpth = 0
        tobin = 0
c apr
c set by default r/wbin=1 which means UNFORMATTED file (binary=TRUE)
        rbin = 1
        wbin = 1
c muta is false by default
        muta = .false.
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c 
c ccrd default
c
        styl1 = 'unkw'
        styl2 = 'unkw'
        lpstr = 1
        lpend = 1
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
         if (find('conn')) then
                 if (find('muta')) muta = .true.
                 ucon  = of()
                 call rconn(ucon)
                 inofrz = npt
                 do i=1,npt
                  nofreez(i) = i
                 end do
         end if
         if (find('rcrd')) urcrd = of()
         if (find('rcr1')) ich   = of()
         if (find('wcrd')) uwcrd = of()
         if (find('rpth')) urpth = of() 
         if (find('wpth')) uwpth = of() 
         if (find('cmbn')) then 
                ucmbn = of()
                lcmbn = .true.
                ncmbn = geti('comb',0)
                styl1 = 'CMBN'
                styl2 = 'CMBN'
         end if
        end if
        if (find('wpck')) then
                call pick(ipick,i)
                npick = 0
                do 10 i=1,npt
                 if (ipick(i).ne.0) then
                  npick = npick + 1
                  jpick(npick) = i
                 end if
10              continue
        end if  
        if (find('opck')) then
                call pick(ipic1,i)
                npic1 = 0
                do 109 i=1,npt
                 if (ipic1(i).ne.0) then
                  npic1 = npic1 + 1
                  jpic1(npic1) = i
                 end if
109             continue
        end if  

cccccccccccccccccccccccccccccccccccccccc
cccccccccccc LIST OF STYLES cccccccccccc
c FROM  TO                             c
c fpth tpth --- PATH unformatted style c
c fdyn tdyn --- DYNA unformatted style c
c fchr tchr --- CHARMM formatted style c
c fgpu tgpu --- DYNA formatted style   c
c ffpt tfpt --- PATH formatted style   c
c fxyz      --- like PATH formatted st c
cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc

        if (find('fpth')) then
                styl1 = 'PATH'
        else if (find('fdyn')) then
                styl1 = 'DYNA'
        else if (find('fchr')) then
                styl1 = 'CHAR'
        else if (find('fxyz')) then
                styl1 = 'PXYZ'
c conversion gpu dcd to bina dcd
        else if (find('fgpu')) then
                styl1 = 'DGPU' 
                rbin = 0   
c conversion formatted path to unformatted path
        else if (find('ffpt')) then
                styl1 = 'FPTH'
        end if
        if (find('tpth')) then
                styl2 = 'PATH'
        else if (find('tdyn')) then
                styl2 = 'DYNA'
        else if (find('tchr')) then
                styl2 = 'CHAR'
        else if (find('tfpt')) then
                styl2 = 'FPTH'
        else if (find('tgpu')) then
                styl2 = 'DGPU'
        end if
        lpstr  = geti('lpst',lpstr)
        lpend = geti('lpen',lpend)
        if (find('wsub')) wsub = .true.
        if (find('ovlp')) ovlp = .true.

c the use of ffpt and tfpt overrides this input 
c- apr convert same file type from bina to formatted or vice versa
c        if (find('cb2f')) then
c          bConv = .true.
c          tobin = 0
c        else if (find('cf2b')) then
c          bConv = .true.
c          tobin = 1
c        endif

        if (find('acti')) go to 2
        go to 1
2       continue
        if (lpstr.lt.1) then
         level = 0
         call alert(name,namel,'lpstr must be > 0',17,level)
         lpstr = 1
        end if
c
c check that required files were opened
c
        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        endif 
c        else if (.not.fopen(urcrd) .and. (.not.lcmbn)) then
c         level = 1
c         call alert(name,namel,'urcrd not opened',16,level)
c        else if (.not.fopen(uwcrd) .and. (.not.lcmbn)) then
c         level = 1
c         call alert(name,namel,'uwcrd not opened',16,level)
c        else if (.not.fopen(urpth) .and. (.not.lcmbn)) then
c         level = 1
c         call alert(name,namel,'uwcrd not opened',16,level)
c        end if

c check that the files styles were given
c
        if (styl1.eq.'unkw' .or. styl2.eq.'unkw') then
         level = 1
         call alert(name,namel,'Missing coordinate style',24,level)
        end if


c conversion from bina path to formatted path
        !if (styl1.eq.'PATH'.and.styl2.eq.'FPTH') then
        ! call convertpath(urpth,uwpth,lpend,0)
        if (styl1.eq.'FPTH'.and.styl2.eq.'PATH') then
         call convertpath(urpth,uwpth,lpend,1)
        endif
c the use of FPTH and PATH overrides this
c        if (bConv) then
c          if (styl1.eq.'PATH') then
c            if (styl2.eq.'PATH') then
c              call convertpath(urpth,uwpth,lpend,tobin)
c            endif
c          endif
c          stop
c        endif  



c conversion from gpu dcd to regurlar dcd
        if (styl1.eq.'DGPU') then
         if (styl2.eq.'CHAR') then
           i=lpstr
           call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
           call putcrd(uwcrd,'CHARM')
           stop
         else if (styl2.eq.'DYNA') then
          do 25 i=lpstr,lpend
           norew = .true.
           call rdyncrd(urcrd,i,inofrz,nofreez,rbin)

           do k = 1, npt
             tmp_mass(k) = ptms(k)
             ptms(k) = 0.e0
           end do
           do k = 1,npic1
             ptms(jpic1(k)) = tmp_mass(jpic1(k))
           end do

           if (ovlp) call rmsd_weight(npt,coor2,coor,rms,.false.,ptms)

           do k = 1, npt
             ptms(k) = tmp_mass(k)
           end do
           call wdyncrd(uwcrd,lpend-lpstr+1,i,inofrz,nofreez,wbin)
25        continue
         else
          level=1
          call alert(name,namel,
     1    'only form to unform dcd or crd allowed',38,level)
         endif
        endif

c       ------
        if (styl1.eq.'PATH') then
         if (styl2.eq.'DYNA') then
          rewind(urcrd)
          do 3 i=lpstr,lpend
           call rpath_seq(urcrd,i)
           call wdyncrd(uwcrd,lpend-lpstr+1,i-lpstr+1,npt,
     1                  0,wbin)
3         continue
         else if (styl2.eq.'CHAR') then
          call rpath(urcrd,lpstr)
          call putcrd(uwcrd,'CHARM')
         else if (styl2.eq.'PATH') then
          do 305 i=lpstr,lpend
           call rpath(urcrd,i)
           call wpath(uwcrd)
305       continue
         else if (styl2.eq.'FPTH') then
          do 306 i=lpstr,lpend
           call rpath(urpth,i)
           call wpathform(uwcrd)
306       continue
         end if

c       ------
        else if (styl1.eq.'PXYZ') then
         if (styl2.eq.'DYNA') then
          do 31 i=lpstr,lpend
           call rpxyz(urcrd,i)
           call wdyncrd(uwcrd,lpend-lpstr+1,i-lpstr+1,inofrz,
     1                  nofreez,wbin)
31        continue
         else if (styl2.eq.'CHAR') then
          call rpxyz(urcrd,lpstr)
          call putcrd(uwcrd,'CHARM')
         end if

c       ------
c apr: add from dyna unformatted to (gpu) formatted 
        else if (styl1.eq.'DYNA') then
         if (styl2.eq.'DYNA') then
          rewind urcrd
          do 35 i=lpstr,lpend
c           norew = .true.
           call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
           do k = 1, npt
             tmp_mass(k) = ptms(k)
             ptms(k) = 0.e0
           end do
           do k = 1,npic1
             ptms(jpic1(k)) = tmp_mass(jpic1(k))
           end do
           
           if (ovlp) call rmsd_weight(npt,coor2,coor,rms,.false.,ptms)

           do k = 1, npt
             ptms(k) = tmp_mass(k) 
           end do
           
           if (wsub) then
                call wsubset(uwcrd,lpend,i,npick,jpick)
           else
c wdyncrd with five arguments is obsolete: see dynamics/wdyncrd
      call wdyncrd(uwcrd,lpend-lpstr+1,i,inofrz,nofreez,wbin)
           end if
35        continue
         else if (styl2.eq.'PATH') then
          do 4 i=lpstr,lpend
            norew = .true.
           call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
           if (npick.ne.npt) then
            call wpath_select(npick,jpick,uwpth)
           else
            call wpath(uwcrd)
           end if
4         continue
         else if (styl2.eq.'CHAR') then
          i = lpstr
          call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
          call putcrd(uwcrd,'CHARM')
         else if (styl2.eq.'DGPU') then
           norew = .true.
           rbin = 1
           wbin = 0 
           do i=lpstr,lpend
             call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
             call wdyncrd(uwcrd,lpend-lpstr+1,i,inofrz,nofreez,wbin)
           enddo    
         end if

c       ------
        else if (styl1.eq.'CHAR') then
c for styl1=CHARMM (i.e. writing from CHARMM to PATH/DYNA files
c lpstr MUST be equal 1 (any other value does not make sense)
c
         if (lpstr.ne.1) then
          level = 1
          call alert(name,namel,'lpst .ne. 1 is illegal ',23,level)
         end if
         if (styl2.eq.'PATH') then
          do 5 i=lpstr,lpend
           call rline(name,namel,stdi)
           ich = of()
           call getcrd(ich,'CHARM')
           call freeunit(ich)
           call wpath(uwcrd)
5         continue
         else if (styl2.eq.'DYNA') then
          do 6 i=lpstr,lpend
           call rline(name,namel,stdi)
           ich = of()
           call getcrd(ich,'CHARM')
           if (i .eq. lpstr) then 
c            save reference structure             
             do j=1,npt
               do l= 1, 3
                 coor2(l,j) = coor(l,j)
               end do
                 nofreez(j) = j
             end do             
           end if
           call freeunit(ich)
c          overlap with respect to the reference structure
           if (ovlp) call rmsd_weight(npt,coor2,coor,rms,.FALSE.,ptms)

cccccccccccccccccccccc
           do k = 1, npt
             tmp_mass(k) = ptms(k)
             ptms(k) = 0.e0
           end do
           do k = 1,npic1
             ptms(jpic1(k)) = tmp_mass(jpic1(k))
           end do

           if (ovlp) call rmsd_weight(npt,coor2,coor,rms,.false.,ptms)

           do k = 1, npt
             ptms(k) = tmp_mass(k)
           end do
cccccccccccccccccccccc

c wdyncrd with five arguments is obsolete: see dynamics/wdyncrd
c           call wdyncrd(uwcrd,lpend,i,npt,wbin)
           call wdyncrd(uwcrd,lpend-lpstr+1,i,npt,nofreez,wbin)
6         continue
         end if
        else if (lcmbn) then
c
c append couple of dynamics file to a single file
c no test are done. It is assumed that the all the files
c do make sense and represent structures of the same 
c connectivity file
c
              ntotal = 0
7             continue
              call rline(name,namel,stdi)
               if (find('*EOD')) then
                if (ncmbn.ne.ntotal) then
                 level = 0
                 call alert(name,namel,'ncmbn != ntotal',15,level)
                 stop
                end if
                stop
               end if
               ich = geti('lpst',0)
               if (ich.eq.0) then
                level = 1
                call alert(name,namel,'Append 0 struc?!',16,level)
               end if
               urcrd = of()
               norew = .false.
               if (ntotal.eq.0) then
                   rewind urcrd
                 do 8 i=1,ich
                   call rdyncrd(urcrd,i,inofrz,nofreez,rbin)
                   norew = .true.
                   call wdyncrd(ucmbn,ncmbn,i,inofrz,nofreez,wbin)
8                continue
                  ntotal = ich
               else
                   rewind urcrd
                 do 9 i=ntotal+1,ntotal+ich
                     call rdyncrd(urcrd,i-ntotal,inofrz,nofreez,rbin)
                   norew = .true.
                     call wdyncrd(ucmbn,ncmbn,i,inofrz,nofreez,wbin)
9                continue
                 ntotal = ntotal + ich
               end if
               go to 7
        end if
        stop
        end

