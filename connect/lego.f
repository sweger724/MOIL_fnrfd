        subroutine lego(upoly,mono,nprt,mono_pt,nbnd,
     1          monob1,monob2,updown,iden,numon,
     2          basenm,basechg,basesgm,baseeps,basems,totunq,
     3          bondp,basekb,basereq,m_charge, m_said)
c
c subroutine lego is the place where the pieces are joined and 
c the BULK is built. Two basic data sets are required
c parameters: particle properties (elec,vdw) and covalent
c               properties (bond, angle, torsion, improper torsions)
c               (see propert.f for more details)
c monomers  : defining the pieces of the macrmolecules in terms of
c               particle properties and covalent binding
c               (see monomers.f for more details)
c
c The routine is reading a list of monomers names with rules of
c how/if to join them. Output: molecular vectors 
c moname  - moname(i) is the name of the i-th monomer in the BULK
c poimon  - poimon(i) is the number of the monomer to which the
c               i-th particle belongs
c poipt   - poipt(i) is the last particle of the i-th monomer
c particle properties vector: ptnm,ptms,ptchg,epsgm6,epsgm12,ptwei
c bond properties: ib1,ib2,kbond,req
c (FOR MORE DETAILS SEE: CONNECT.BLOCK)
c
c This routine may be called a number of times to make a number
c of molecular set
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer upoly
        character*4 mono(*),mono_pt(*),basenm(*)
        integer nprt(*),nbnd(*),monob1(*),monob2(*),updown(*)
        integer iden(*),bondp(maxubnd,2)
        integer numon,totunq
        double precision basechg(*),basesgm(*),baseeps(*),basems(*)
        double precision basekb(*),basereq(*),m_charge(*)
        integer m_said(*)

        character*80 getchar
        integer geti
        logical find,first
c local
        character*4 name,mono1,mono2
        logical empty,les1,les2
        integer level,imult,index,ityp1,i1,idb1,ipart1
        integer kb1,ityp2,ipart2,kb2
        integer i,namel
        integer copies1,copies2
        integer index_tmp,lestyp1,lestyp2

c Data base to manipulate connectivity and particle
c definition of the monomers
c _m1 refer to monomer 1 ; _m2 refer to monomer 2
c ptnm(i) is a unique name of particle i within a monomer
c id(i)   gives the i-th particle property index
c ib1(i) & ib2(i) paticle indices of the i-th bond
c updown family relation vector. See monomers.f for detail
c       explanation
c
        character*4 ptnm_m1(maxptmo),ptnm_m2(maxptmo)
        integer id_m1(maxptmo),id_m2(maxptmo)
        integer ib1_m1(maxptmo),ib2_m1(maxptmo),ib1_m2(maxptmo)
        integer ib2_m2(maxptmo),updwn_m1(maxptmo),updwn_m2(maxptmo)
        double precision weight1(maxcopy),weight2(maxcopy)
        double precision m1_charge(maxptmo),m2_charge(maxptmo)
        data weight1,weight2/maxcopy*0.d0,maxcopy*0.d0/
c YS 
        integer m1_said(maxptmo), m2_said(maxptmo) 



        name = 'LEGO'
        namel=4

c--------------------------------------------------------------------
c GET NAME OF MOLECULE TYPE OF GLUE AND COVER RESIDUES
c label 1 is where a new molecule should start
c

1       continue
        call rline(name,namel,upoly)
        if (find('*EOD')) then
c        if (debug) write(stdo,*)(bulk(i),i=1,nbulk)
         return
        end if
c
c Get a name for your BULK
c
        NBULK = NBULK + 1
        IF (NBULK.GT.MAXBULK) THEN
         level = 1
         call alert(name,namel,' Too many molecules ',20,level)
         return
        end if
        BULK(NBULK)(1:4) = getchar('MOLC','BULK',4)
c--------------------------------------------------------------------

c--------------------------------------------------------------------
c set number of copies to zero
c
        imult = 0 
c--------------------------------------------------------------------

c--------------------------------------------------------------------
c set index (current monomer number) to totmon. If this is the first
c molecular set totmon.eq.0, then
c get the number of monomers in the current molecular set
c Note that currently index is set to the old totmon (for a new
c system it is zero). Index is updated only at getmoname.
c
        index = totmon 
        totmon = totmon + geti('#mon',0)
        if (totmon.gt.maxmono) then
         level = 1
         call alert(name,namel,'Max. # of monomers exceed',26,level)
        end if
c       if (debug) then
c        write(stdo,*) 'NBULK totmon = ',NBULK,totmon
c       end if
        if (totmon-index.eq.0) then
         level = 1
         call alert(name,namel,'Missing number of monomers',26,level)
        end if
c--------------------------------------------------------------------


c--------------------------------------------------------------------
c start loop on monomers by reading first line and getting first
c monomer name
c
         call rline(name,namel,upoly)
c        if (debug) then
c         write(stdo,*)' numon mono(1) ',numon,mono(1)
c         write(stdo,*)' mono ',(mono(i),i=1,numon)
c        end if
         call getmoname(mono1,imult,les1,empty,copies1,mono
     1    ,moname,index,ityp1,numon,lestyp1,weight1)
c set a general flag on existnce of LES monomers
         if (les1) lesflag = .true.
         index_tmp = index
         if (empty) then
          level = 1
          call alert(name,namel,'Missing first mono name',24,level)
         end if
c--------------------------------------------------------------------


c--------------------------------------------------------------------
c if totmon .eq. index (the total number of monomer in the current
c molecular set) then our procedure is rather short:
c It includes only one monomer and therefore no special paerticles.
c (some of the work usually done by "dimer" need to be done here)
c
        if (totmon.eq.index) then
          if (ityp1.eq.1) then
           i1  = 1
           idb1 = 1
           ipart1 = nprt(1)
           kb1    = nbnd(1)
          else
           i1 = nprt(ityp1-1) + 1
           idb1 = nbnd(ityp1-1) + 1
           ipart1 = nprt(ityp1) - nprt(ityp1-1)
           kb1    = nbnd(ityp1) - nbnd(ityp1-1)
          end if
          call vicopy(iden(i1),id_m1,ipart1)
          call vicopy(updown(i1),updwn_m1,ipart1)
          call vccopy(mono_pt(i1),ptnm_m1,ipart1,4)
          call vicopy(monob1(idb1),ib1_m1,kb1)
          call vicopy(monob2(idb1),ib2_m1,kb1)
          call vdcopy(m_charge(i1),m1_charge,ipart1)
          call vicopy(m_said(i1), m1_said, ipart1)
          call fillpt(ipart1,ityp1,ptnm_m1,id_m1,ib1_m1,ib2_m1,
     1     kb1,les1,copies1,index,lestyp1,weight1,m1_charge,m1_said)
          pbulk(NBULK) = npt
          write(stdo,100)BULK(NBULK)
          go to 1
        end if
c--------------------------------------------------------------------
         
c------------------------------------------------------------------
c Finally THE loop on monomers (get second monomer)
c label 2 is where reading a sequence of monomers of the
c same molecule is done
c

c first is true on first call to dimer. In this case dimer initiate
c the vectors of monomer 1 too.
        first = .true.
2       continue

c       if (debug) then
c               write(stdo,*) ' before getting a second monomer '
c               write(stdo,*) ' index totmon ',index,totmon
c       end if
        call getmoname(mono2,imult,les2,empty,copies2,mono
     1   ,moname,index,ityp2,numon,lestyp2,weight2)
        if (debug) write(*,*) ' index = ',index
c
c Set a general flag on the existence of LES monomers
        if (les2) lesflag = .true.
c       if (debug) then
c        write(stdo,*) 'empty = ',empty
c       end if
c------------------------------------------------------------------

c------------------------------------------------------------------
c If the current line is empty the following reasons/actions should
c be considered:
c (a) we are done and ready to fill the last monomer properties 
c (b) we need to read a new line
c (c) one more BUG was found
c
        if (empty) then
         if (totmon.eq.index) then
c         if (debug) then
c          write(stdo,*)' in lego baseeps basesgm '
c     1                    ,baseeps(1),basesgm(1)
c         end if
c this is the last monomer (we use 1 rather than two since two was
c copied to one) and the last (unsuccessful) call to getmoname
c is likely to change some of the variables (i.e. les2)
c
          call fillpt(ipart1,ityp1,ptnm_m1,id_m1,ib1_m1,ib2_m1,kb1,
     1     les1,copies1,index_tmp,lestyp1,weight1,m1_charge,m1_said)
c         if (debug) then
c          write(stdo,*) ' in lego epsgm6 epsgm12 '
c     1         ,(epsgm6(ij),ij=1,npt),(epsgm12(ij),ij=1,npt)
c         end if
          pbulk(NBULK) = npt
          write(stdo,100)BULK(NBULK)
100       format(/,1x,'*** MOLECULE ',a4, ' : BONDS COMPLETED ***'
     1     ,//)
          go to 1
         else if (index.gt.totmon) then
          level = 1
          write(stdo,101) totmon,index
101       format(1x,' Expected (read) number of mono',i7,' found '
     1     ,i7)
          call alert(name,namel,'Number of monomers exceed',27,level)
         else
          call rline(name,namel,upoly)
          go to 2
         end if
        end if
c------------------------------------------------------------------

c------------------------------------------------------------------
c wedding. attach two monomers to form a dimer
c
c *** WARNING ***
c *** In dimer it is assumed that interface particles are for two monomers
c *** only. An interface particle for three monomers is illegal
c
c       if (debug) then
c        write(stdo,*) ' Before dimer '
c        write(stdo,*) ' mono_pt(1) ',mono_pt(1)
c        write(stdo,*) ' ib1(1) ib2(1) ',ib1(1),ib2(1)
c        write(stdo,*) ' ityp1 ityp2 ',ityp1,ityp2
c        write(stdo,*) ' nbnd = ',(nbnd(i),i=1,2)
c        write(stdo,*) ' id_m1 ',(id_m1(iii),iii=1,ipart1)
c        write(stdo,*) ' updwn_m2 ',(updwn_m2(iii),iii=1,ipart2)
c       end if
        call dimer(ityp1,ityp2,mono_pt,updown,nprt,nbnd,iden
     1          ,monob1,monob2,ptnm_m1,ptnm_m2
     1   ,id_m1,id_m2,ib1_m1,ib2_m1,ib1_m2,ib2_m2
     2   ,updwn_m1,updwn_m2,ipart1,ipart2,kb1,kb2,first,
     3   m1_charge,m2_charge,m_charge,
     4   m1_said,m2_said,m_said)
        first = .false.
c       if (debug) then
c               write(stdo,*) ' After dimer '
c               write(stdo,*)'id_m1 ',(id_m1(iii),iii=1,ipart1)
c               write(stdo,*)'id_m2 ',(id_m2(iii),iii=1,ipart2)
c        write(stdo,*) ' updwn_m2 ',(updwn_m2(iii),iii=1,ipart2)
c               write(stdo,*) ' kb1 kb2 ',kb1,kb2
c       end if
c------------------------------------------------------------------

c------------------------------------------------------------------
c fill particle properties for monomer number one 
c since this was at least the second call to getmoname
c and we are filling the previous one, the index is set to
c index_tmp
c
        call fillpt(ipart1,ityp1,ptnm_m1,id_m1,ib1_m1,ib2_m1,kb1,
     1     les1,copies1,index_tmp,lestyp1,weight1,m1_charge,m1_said)
c------------------------------------------------------------------

c------------------------------------------------------------------
c Ready to continue to read another monomer but first
c fill monmomer 1 properties with those of monomer 2
c and then go to 2 to read a new monomer 2
c
        les1    = les2
        copies1 = copies2
        ityp1   = ityp2
        ipart1  = ipart2
        kb1     = kb2
        index_tmp = index
        lestyp1   = lestyp2
        if (les1 .and. weight2(1).gt.0) then
         do 31 i=1,maxcopy
          weight1(i) = weight2(i)
31       continue
        end if
c       if (debug) then
c        write(stdo,*)' copying monomer 2 to monomer 1'
c       end if
        do 3 i=1,ipart2
         ptnm_m1(i) = ptnm_m2(i)
         id_m1(i)   = id_m2(i)
         updwn_m1(i)= updwn_m2(i)
         m1_charge(i) = m2_charge(i)
         m1_said(i) = m2_said(i)
3       continue
c       if (debug) then
c        write (stdo,*)' after copying monomer 2 to monomer 1'
c        write(stdo,*)' id_m1 = ',(id_m1(iii),iii=1,ipart1)
c       end if
        do 4 i=1,kb2
         ib1_m1(i)  = ib1_m2(i)
         ib2_m1(i)  = ib2_m2(i)
4       continue
c------------------------------------------------------------------
                
c------------------------------------------------------------------
c Go to read another monomer
c
        go to 2
        end















