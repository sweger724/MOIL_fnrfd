        subroutine dimer(ityp1,ityp2,mono_pt,updown,nprt,nbnd,iden,
     1   ib1,ib2,
     2   ptnm_m1,ptnm_m2,id_m1,id_m2,ib1_m1,ib2_m1,ib1_m2,ib2_m2,
     3   updwn_m1,updwn_m2,ipart1,ipart2,nb1,nb2,first,m1_charge,
     4   m2_charge,m_charge,
     5   m1_said,m2_said,m_said)
c
c subroutine to generate a list of connections between the monomers
c INPUT
c ityp1   - index of a monomer in the data-base (first monomer)
c ityp2   - index of a second monomer in the data base
c mono_pt - char*4 vector with unique particle names in the
c               different monomers
c nprt    - a pointer to the last particle of a monomer in the
c               monomer connectivity list. I.e. nprt(j) is the
c               last particle of the j-th monomer
c nbnd    - a pointer to the last bond in a monomer, nbnd(j)
c               is the number of the last atom of j-th mono
c iden    - a list of particle types according to nprt pointer
c ib1 ib2 - list of bonds according to nbnd pointer
c updown  - direction for connetivity to interface between
c           monomers: Particles in the j-th monomer are between
c           nprt(j-1)+1 to nprt(j) and updown(nprt(j-1))+1
c           to updown(nprt(j)) can have the following values
c             0 : leave as is (normal atom)
c             1 : virtual atom. actually exists in the next monomer
c               according to unique atom name the virtual atom is
c               identified in the next monomer 
c            -1 : virtual atom. actually exists in the previous monomer
c           999 : atom to be deleted in the next monomer
c          -999 : atom to be deleted in the previous monomer
c OUTPUT
c the following vectors get filled
c id_m1 & id_m2 - identification index for particles in the
c                 two monomers. "Virtual" particles are set to
c                 id = -id. Only Bonds of virtual particles are kept.
c                 Particles to be deleted are set to id = 0.
c ib1_m1 & ib1_m1 - the bond list is according to the atom number in
c ib2_m1 & ib2_m2 - the two monomers. To increase the confusion
c                   negative address is also used to reflect
c                   bonding to previous monomers. Negative address
c                   is added only to the second monomer but the
c                   first monomer may have negative address too
c                   from previous marriage
c ptnm_m1 & ptnm_m2 - compressed lists of particles names in the
c                     the two monomers.
c ipart1 ipart1     - number of particles in monomers 1 & 2
c nb1    nb2        - number of bonds in monomers 1 and 2
c
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        logical first
        integer ityp1,ityp2,ipart1,ipart2,kb1,kb2,nb1,nb2
        integer nprt(*),nbnd(*),updown(*),iden(*),ib1(*),ib2(*)
        integer id_m1(*),id_m2(*),ib1_m1(*)
        integer ib2_m1(*),ib1_m2(*),ib2_m2(*)
        integer updwn_m1(*),updwn_m2(*)
        character*4 mono_pt(*),ptnm_m1(*),ptnm_m2(*)
        double precision m1_charge(*),m2_charge(*),m_charge(*)
        integer m1_said(*), m2_said(*), m_said(*)
c
c local
c
        integer i1,i2,k1,k2,ikb1,ikb2,namel,level
        character*5 name
        data name/'dimer'/
        data namel/5/
        
        
c get pointers to monomer arrays
c

c The vectors of monomer 1 need to be filled only at the first call
c to dimer
        if (first) then
         if (ityp1.eq.1) then
          k1  = 1
          kb1 = 1
         else
          k1  = nprt(ityp1-1)+1
          kb1 = nbnd(ityp1-1)+1
         end if
         ipart1 = nprt(ityp1) - k1 + 1
         nb1    = nbnd(ityp1) - kb1 + 1
         call vicopy(iden(k1),id_m1,ipart1)
         call vicopy(updown(k1),updwn_m1,ipart1)
         call vccopy(mono_pt(k1),ptnm_m1,ipart1,4)
         call vdcopy(m_charge(k1),m1_charge,ipart1)
         call vicopy(m_said(k1),m1_said,ipart1)
         if (nb1.gt.0) then
          call vicopy(ib1(kb1),ib1_m1,nb1)
          call vicopy(ib2(kb1),ib2_m1,nb1)
         end if
        end if

        if (ityp2.eq.1) then
         k2  = 1
         kb2 = 1
        else
         k2  = nprt(ityp2-1)+1
         kb2 = nbnd(ityp2-1)+1
        end if
        ipart2 = nprt(ityp2) - k2 + 1
        nb2    = nbnd(ityp2) - kb2 + 1
        if (debug) then
         write(stdo,*)' in dimer kb2 ',kb2
        end if
        if (ipart1.le.0 .or. ipart2.le.0) then
         if (debug) then
          write (stdo,*)' ityp1 ityp2 ',ityp1,ityp2
         end if
         write(stdo,100)ipart1,ipart2,nb1,nb2
100      format(1x,'*** ILLEGAL VALUES ',/,1x,
     1          ' length of first and second ',i7,i7,
     2          ' bond number of first and second ',i7,i7)
         level = 1
         call alert(name,namel,' lists are not defined ',23,level)
        end if

c fill the pair id vectors and bond vectors
c vicopy(x,y,l) copy the content of an integer vector x
c with length l to y
c vccopy(x,y,l,size) copy the content of the character vector
c  x length l to y. Each element size of "size"
c
        if (debug) then
         write(stdo,*)' ipart1 ipart2 k1 k2 ',ipart1,ipart2,k1,k2
         write(stdo,*)' mono_pt(k2) '
         write(stdo,*)mono_pt(k2)
        end if
        call vicopy(iden(k2),id_m2,ipart2)
        call vicopy(updown(k2),updwn_m2,ipart2)
        call vccopy(mono_pt(k2),ptnm_m2,ipart2,4)
        call vdcopy(m_charge(k2),m2_charge,ipart2)
        call vicopy(m_said(k2),m2_said,ipart2)

        if (debug) then
         write(stdo,*)' After second vccopy '
         write(stdo,*)' ib1(kb2) ib2(kb2) '
         write(stdo,*) ib1(kb2),ib2(kb2)
        end if

        if (nb2.gt.0) then
        call vicopy(ib1(kb2),ib1_m2,nb2)
        call vicopy(ib2(kb2),ib2_m2,nb2)
        end if

c start loop on particles of monomer1 check on updown of 1
c and verify that no instruction of monomer 2 refer to present
c particle
c
        do 9 i1=1,ipart1

c
c if particle id = 0  ->  this particle was erased.
c
         if (id_m1(i1).eq.0) go to 9
c
c Is the first monomer list suggests that the i-th particle is ok?
         if (updwn_m1(i1).eq.0) then


          do 3 i2=1,ipart2
c
c Does the second monomer list have abnormal particles
c (updown(j)!=0) ?
c and does the abnormal particle concides with the first monomer
c particle?

           if ((updwn_m2(i2) .ne. 0) .and.
     1        (ptnm_m1(i1) .eq. ptnm_m2(i2))) then
c
c If yes, do we need to erase this particle (updown=-999)?
             if (updwn_m2(i2).eq.-999) then  
                if (debug) then
                 write(stdo,*) ' delete previous particle was called '
                 write(stdo,*) ' particle name ',ptnm_m1(i1)
                end if
c
c If yes take the following actions:
c (a) set id_m1(i1) & id_m2(i2) to zero
c (b) erase bonds which include particles i1/i2 and reduce bond list
                id_m1(i1) = 0
                id_m2(i2) = 0
c
c cmprs2(x,y,erase,l) :
c if (x(i) or y(i) = erase ) then the vectors x and y are comprssed
c to remove the "erase". The length l is changing on output
                if (nb1.gt.0) then
                 call cmprs2(ib1_m1,ib2_m1,i1,nb1)
                end if
                if (nb2.gt.0) then
                 call cmprs2(ib1_m2,ib2_m2,i2,nb2)
                end if
c
c The PREVious command (when updwn=-1) set the id of the particle
c in monomer *** ONE ***
c to that of monomer two....
c A tricky point, becareful
             else if (updwn_m2(i2).eq.-1) then
c
c set the particle type in monomer 1 to that of monomer 2. Then
c set the id of 2 to virtual (-1). On fillpt particle 2 will be
c ignored and only bonds will be taken into account.
                if (debug) then
                 write(stdo,*)' Previous particle found '
                 write(stdo,*)' id_m1(i1) id_m2(i2) '
                 write(stdo,*) id_m1(i1),id_m2(i2)
                end if
                id_m1(i1) = id_m2(i2)
            if (m2_charge(i2) .ne. 999.d0) m1_charge(i1) = m2_charge(i2)
                id_m2(i2) = -1
c
c     all the i2 bonds will have -i1 instead of i2. Note that i1
c     is a relative position in monomer 1 and care should be taken
c     when implementing to the molecular list.
c

                do 2 ikb2=1,nb2
                 if (ib1_m2(ikb2).eq.i2) ib1_m2(ikb2) = -i1
                 if (ib2_m2(ikb2).eq.i2) ib2_m2(ikb2) = -i1
2               continue
             end if
           end if
3         continue

c
c check if i1 requires special treatment
c if updwn=1 then i1 is a virtual particle all its bonds should
c be transfered to the next monomer
c the id of the NEXT particles is set to that of the virtual
c
         else if (updwn_m1(i1).eq.1) then
c 
          if (debug) then
           write(stdo,*) ' NeXT particle was found '
           write(stdo,*) 'Current bond number ',kb1
           write(stdo,*) ' ptnm(',i1,') = ',ptnm_m1(i1)
          end if
          do 5 i2=1,ipart2
           if (debug) then
            write(stdo,*) 'ptnm_m1 ptnm_m2 updwn_m2 '
            write(stdo,*)ptnm_m1(i1),ptnm_m2(i2),updwn_m2(i2)
           end if
c check that the two particles have identical names,
c but i2 may be special...c in this case ignore the match
           if (ptnm_m1(i1) .eq. ptnm_m2(i2) .and. 
     1     updwn_m2(i2).eq.0) then
            do 4 ikb1=1,nb1
             if (ib1_m1(ikb1).eq.i1) then
              nb2 = nb2 + 1
              ib1_m2(nb2) = i2
              ib2_m2(nb2) = -ib2_m1(ikb1)
             else if (ib2_m1(ikb1).eq.i1) then
              nb2 = nb2 + 1
              ib1_m2(nb2) = i2
              ib2_m2(nb2) = -ib1_m1(ikb1)
             end if
4           continue
c
c set id of i2 to that of i1
            id_m2(i2) = id_m1(i1)
         if (m1_charge(i1) .ne. 999.d0) m2_charge(i2) = m1_charge(i1)
c 
c delete particle i1
            id_m1(i1) = 0
c remove bonds of particle i1
            call cmprs2(ib1_m1,ib2_m1,i1,nb1)
            if (debug ) then
             write(stdo,*)' After compressing NeXT particle'
             write(stdo,*) 'Current bond number ',kb1
            end if
            go to 9
           end if
5         continue
          write(stdo,101)ptnm_m1(i1)
101       format(1x,'Virtual particle not found ',1x,a4)
c delete this particle and erase all its bonds
c
          if (debug) then
                write(stdo,*)' i1 id_m1(i1) ',i1,id_m1(i1)
          end if
          id_m1(i1) = 0
          call cmprs2(ib1_m1,ib2_m1,i1,nb1)
          level = 0
          call alert(name,namel,
     1     'Matching with second monomer unsuccesful',40,level)

c
c check if i1 particle is a virtual particle to be deleted at 2
         else if (updwn_m1(i1).eq.999) then
          do 7 i2=1,ipart2
           if (ptnm_m1(i1).eq.ptnm_m2(i2)) then
                if (debug) then
                 write(stdo,*)' DNXT was called for particles:'
                 write(stdo,*)' i1 i2 = ',i1,i2
                 write(stdo,*)ptnm_m1(i1),' of first monomer '
                 write(stdo,*)ptnm_m2(i2),' of second monomer '
                end if
c delete i1 & i2
            id_m1(i1) = 0
            id_m2(i2) = 0
c compress bond list
            if (nb1.gt.0) then
            call cmprs2(ib1_m1,ib2_m1,i1,nb1)
            end if
            if (nb2.gt.0) then
            call cmprs2(ib1_m2,ib2_m2,i2,nb2)
            end if
            go to 8
           end if
7         continue
          write(stdo,102)ptnm_m1(i1)
102       format(1x,' Particle to be deleted has no match ',a4)
          level = 1
          call alert(name,namel,
     1     'Matching with second monomer unsuccesful',40,level)
8         continue
         end if

9        continue
         return
         end
