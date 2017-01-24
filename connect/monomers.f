        subroutine monomers(ucon,mono,nprt,mono_pt,nbnd,basenm,
     1  basechg,ib1,ib2,updown,iden,maxumon,totunq,numon,maxdiv,
     2  monodiv,mdivyes,m_charge, m_said)
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/LINE.BLOCK'
c
c subroutine to read monomer properties
c mono   - character*4 vector for name of monomers
c basenm - particle names according to property index
c basechg- particle charges according to property index
c nprt   - intg pointer to particle indices in the monomer.
c               nprt(i) is the last particle of the i-th monomer.
c nbnd   - intg pointer. nbnd(i) is the position of the last bond of
c               the i-th monomer in the ib? vectors.
c iden   - a pointer of a particle in a monomer to property index
c          iden(nprt(k-1)+i) is the property index of the i-th atom in
c          the k-th monomer.
c mono_pt   - a vector (char*4) of particle names. mono_pt(nprt(k-1)+i)
c          is the name of the i-th particle in the k-th monomer
c ib1    - pointer to bonds of individual monomers - ib1(nbnd(k-1)+i) is
c          the i-th bond in monomer k
c ib2    - same as above, second atom
c updown - intg vector for particle family relation
c          updown=-1 particle of the residue before
c          updown=0 particle defined at the present residue,
c          updown=1 particle defined at the residue after.
c          updown=999 remove particle from next monomer
c          updown=-999 remove particle from previous monomer.
c totunq - total number of property vector elements
c numon  - actual number of unique monomers
c m_said   - atom type id for surface area calculation
c
        integer totunq,numon,maxumon
        integer ucon
        character*4 mono(*),mono_pt(*),basenm(*)
        integer nprt(*),nbnd(*)
        integer iden(*),ib1(*),ib2(*)
        integer updown(*)
        double precision basechg(*),m_charge(*)
        integer m_said(*)
        integer maxdiv,monodiv(0:maxdiv+1,maxumon)
        logical mdivyes


c local
        integer l1,l2,k1
        integer ibnd,iprt,mainlp
        integer i,j,previous,iunq,level
        integer namel,lc
        integer curridx,lastidx,fdivnum,lastprt

        logical empty,special,test
        double precision tchrg,tchrg1
        character*8 name,chbnd
        character*4 prtc,char_tmp
        
        integer geti
        double precision getd
        logical find
        character*5 getchar

        name = 'monomers'
        namel = 8

c
c read head for MONOMERS
c
        silent = .true.
        call rline(name,namel,ucon)
        if (.not. (find('MONO').and.find('LIST')) ) then
         level = 1
         call alert(name,namel,'Missing monomers header',23,level)
        end if
c
c syntax for defining a monomer
c  name     number of particles  charge 
c MONO=(name)   #prt=intg      chrg=float
c

c previous - number of previous particles
c ibnd     - total number of current bonds

        previous = 0
        ibnd = 0

        do 7 mainlp=1,maxumon
c
c For test purposes
c iprt  - number of particles found
c tchrg - total charge read from input
c tchrg1 - total charge calculated from properties vector of
c       particles for test purposes (tchrg1.eq.tchrg)
c
        iprt = 0
        tchrg = 0.d0
        tchrg1 = 0.d0
        lastidx = 0

          call rline(name,namel,ucon)
          if (find('*EOD')) then
c
c *EOD - sign for End Of Data
c
           numon = mainlp - 1
           previous = previous + iprt
           go to 8
          end if
          lc = 4
          mono(mainlp)   =getchar('MONO','????',lc)
          iprt           =geti('#prt',-1)
          tchrg          =getd('chrg',999.0d0)
          fdivnum        =geti('#div',0)

          nprt(mainlp) = previous + iprt
          monodiv(0,mainlp) = fdivnum
          monodiv(1,mainlp) = 0
          lastprt = 0
          
          if (mono(mainlp).eq.'????') then
                level = 1
                call alert(name,namel,'No name for monomer',19,level)
          else if (iprt.lt.0) then
                level = 1
                call alert(name,namel,'Missing #prt',12,level)
          else if (tchrg.gt.998) then
                level = 1
                call alert(name,namel,'Missing chrg',12,level)
          else if (fdivnum.gt.maxdiv) then
                level = 1
                call alert(name,namel,'#div too high',13,level)
          end if
c
c syntax for defining a particle
c Unique name (per monomers)    name according to particle properties
c       UNIQ=(name)              PRTC=(name)  PCHG=X.XX UPDOWN_WORD
c
          do 2 j=previous+1,nprt(mainlp)

           call rline(name,namel,ucon)

           lc = 4
           char_tmp   = ' '

           if (mdivyes) curridx = geti('divi',0)

           char_tmp   = getchar('UNIQ','????',lc)
           mono_pt(j)(1:4) = char_tmp(1:4)
           lc = 4
           char_tmp   = ' '
           char_tmp   = getchar('PRTC','????',lc)
           prtc(1:4)  = char_tmp
           char_tmp   = ' '
           m_charge(j) = getd('PCHG',999.d0)
           m_said(j) = geti('SAID', 1)
           call get4c(char_tmp,empty)

           if (mono_pt(j).eq.'????') then
            level = 1
            call alert(name,namel,'Missing particle name',21,level)
           else if (prtc.eq.'????') then
            level =1
            call alert(name,namel,'Missing unique prtc name',24,level)
           end if

           if (empty .or. char_tmp .eq. 'HERE') then
                updown(j) = 0
           else if (char_tmp .eq. 'NEXT') then
                updown(j) = 1
           else if (char_tmp .eq. 'PREV') then
                updown(j) = -1
           else if (char_tmp .eq. 'DNXT') then
                updown(j) = 999
           else if (char_tmp .eq. 'DPRV') then
                updown(j) = -999
           end if
        
C --- Yael get the division number -----
         if (updown(j).eq.0) then  
           if (mdivyes) then
             if ((lastidx.eq.-1).and.(curridx.ne.0)) then
               level = 1
               call alert(name,namel,'wrong divi number',20,level)
             else if (curridx.gt.fdivnum) then
               level = 1
               call alert(name,namel,'divi larger than limit',25,level)
             endif
             if ((curridx.eq.lastidx+1).and.lastidx.ne.-1) then
               monodiv(curridx,mainlp) = lastprt
               lastidx = curridx
             else if (curridx.eq.0) then
               if (lastidx.ne.0) then
                monodiv(lastidx+1,mainlp) = lastprt
                lastidx = -1
               endif
             else if (curridx.ne.lastidx) then
             level = 1
             call alert(name,namel,'divi numbers not in order',30,level)
             endif
          endif
          lastprt = lastprt + 1
        endif
C ---------------------------------------

c
c get the particle property index
c
           do 1 iunq=1,totunq
                if (basenm(iunq).eq.prtc) then
                        iden(j)=iunq
                        go to 2
                end if
1          continue
           level = 1
           write(*,*)' PRTC = ',prtc
           call alert(name,namel,'PRTC not found',14,level)

2         continue
c yael
         if (mdivyes) then
           if (lastidx.gt.0) then
             if (lastidx.ne.fdivnum) then
               level = 1
               call alert(name,namel,'No. of mon parts < than defined',
     1                  21,level)
             endif

           endif
           if (fdivnum.eq.0) then
              lastidx = 1
              monodiv(0,mainlp) = 1
           endif
             if (lastidx.ne.-1) then
               monodiv(lastidx+1,mainlp) = lastprt
           else
              monodiv(0,mainlp) = -1
           endif
         endif

c
c check total charge
c
          tchrg1 = 0.d0
          do 3 i=previous+1,nprt(mainlp)
                if (m_charge(i).lt.999.d0) then
                 tchrg1 = tchrg1 + m_charge(i)
                else
                 tchrg1 = tchrg1 + basechg(iden(i))
                end if
3         continue
          if (dabs(tchrg-tchrg1).gt.1.d-6) then
                level = 1
                write(stdo,100)tchrg,tchrg1
100             format(1x,' Given monomer charge: ',f10.5,
     1                ' Calculated: ',f10.5)
                call alert(name,namel,'Charges do not match',20,level)
          end if

          call rline(name,namel,ucon)
          if (.not.find('DONE')) then
                level = 0
                call alert(name,namel,'Missing DONE',12,level)
          end if

          call rline(name,namel,ucon)
          if (.not.find('BOND')) then
                level = 0
                call alert(name,namel,'Missing BOND',12,level)
                if (mainlp.eq.1) then
                 nbnd(mainlp) = 0
                else
                 nbnd(mainlp) =nbnd(mainlp-1)
                end if
                if (find('*EOD')) then
                        numon = mainlp 
                        previous = previous + iprt
                        go to 8
                else if (find('DONE')) then
                        previous = previous + iprt
                        go to 7
                else
                 level = 1
                 call alert(name,namel,'Illegal syntax',14,level)
                end if
          end if

c
c syntax bonds is as follows: two characters separated by a dash "-"
c e.g. CA-CB CP1H-CP2H OXXX-HXXX
c maximum length of an atom name - 4 characters
c If no bond was found character length will be negative.
c     DO NOT USE "-" AS A PART OF AN ATOM NAME.
c
        call rline(name,namel,ucon)
4       continue

        call brkpair(chbnd(1:4),l1,chbnd(5:8),l2,special)
c
c Is the current line empty from bonds? If yes read new line
c
        if (l1 .lt. 0 .or. l2 .lt. 0 ) go to 5
        ibnd = ibnd + 1
c
c set default bond id's to -100 (particle not found)
c check if both particles are in the particles id list
c
        ib1(ibnd) = -100
        ib2(ibnd) = -100
        do 6 k1=1,iprt
         test = (updown(previous+k1).eq.0 .and. (.not.special))
     1          .or. (updown(previous+k1).ne.0 .and.special)
         if (debug .and .special) write(stdo,*)' special = true '
         if (mono_pt(previous+k1)(1:4).eq.chbnd(1:4).and.
     1          updown(previous+k1).eq.0) then
                ib1(ibnd) = k1
         else if (mono_pt(previous+k1)(1:4).eq.chbnd(5:8)
     1          .and.test) then
                if (debug) then
                 write(stdo,*)' bond was found'
                 write(stdo,*)' particle ',mono_pt(previous+k1)
                 write(stdo,*)' index = ',iprt
                end if
                ib2(ibnd) = k1
         end if
6       continue
        if (ib1(ibnd).lt.0 .or. ib2(ibnd).lt.0) then
                level = 1
                write(stdo,101)chbnd(1:8)
101             format(1x,' Funny input to bond ',a8)
                call alert(name,namel,'Unrecognized atom',17,level)
        end if
        if (ib1(ibnd).eq.ib2(ibnd)) then
         level = 1
         write(stdo,102)chbnd(1:8)
102      format(1x,' Atom bonded to itself ',a8)
         call alert(name,namel,'Illegal bonding',15,level)
        end if
c
c Are all the bonds extracted from line ? Try once more
        go to 4
c
c test for word DONE
c
5       continue
        call rline(name,namel,ucon)
        if (.not.find('DONE')) go to 4
        nbnd(mainlp) = ibnd
        previous = previous + iprt
7       continue
        numon = maxumon
        level = 0
        call alert(name,namel,'Mono vectors are full!!!',24,level)
        call rline(name,namel,ucon)
        if (.not.find('*EOD')) then
         level = 0
      call alert(name,namel,'Missing *EOD statement for mono',31,level)
        end if
8       continue
C
C*** END OF MONOMERS DATA
C

        return

        end
