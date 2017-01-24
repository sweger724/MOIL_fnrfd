        subroutine fillpt(ipart1,ityp1,ptnm_m1,id_m1,ib1_m1,ib2_m1
     1   ,kb1,les1,copies,index,lestype,weight1,m1_charge,m1_said)

c Given the family relation defined between the monomers -
c fill the arrays for one monomer in the molecular array.

c Input: 
c ipart1 - number of particles in monomer 1 (including deletions)
c ityp1  - type of monomer
c ptnm_m1 - a vector with unique names for particles in monomer 1
c id_m1   - a vector with atom types for monomer 1
c ib1_m1 & ib2_m1 - vectors of bonds ib1_m1(i) & ib2_m2(i) are the
c                       two bonded atoms of the i-th bond
c kb1 - number of bonds for the present monomer
c les1 - a logical flag indicating if LES is activated if activated
c        the present monomer is multiplied and added as "copies" copies
c copies - the number of monomer copies used.

c INPUT & OUTPUT (in CONNECT.BLOCK)
c npt - current number of particles in the molecule
c nb  - total number of bonds for the present molecule
c index  - current number of monomers updated before at fillpt 
c lestyp - current number of different LES type monomers

c Output:
c (Most of the output is in the CONNECT.BLOCK)
c poimon - integer pointer (poimon(ipt) the monomer id of ipt particle)
c ptnm   - (char*4) name of the particle as a function of the
c          particle number
c ptms   - (double) mass of the particle
c ptchg  - (double) charge of the particle
c epsgm6 - (double) twice the square root of van der Waals energy
c               times sigma^6
c epsgm12- (double) twice the square root of van der Waals energy
c               times sigma^12
c ptwei  - (double) partical weight, used in LES calculations 
c lesid  - integer id of LES particle. Normal particle = 0
c               LES particles = les_indx. For the same les_indx
c               particles see each other in full.
c moname - character vector with the name of monomer types in the molecule
c poipt  - a pointer. poipt(i) is the last atom of monomer i
c ib1 ib2- pointers to bonds ib?(k) is a ? particle number that participate in
c               the k-th bond
c kbond  - (double) force constant for bonds.
c req    - (double) equilibrium position for harmonic energy of bonds
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PROPERT.BLOCK'
        include 'COMMON/MONOMERS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        logical les1
        integer search
        integer ipart1,ityp1,kb1,copies,lestype,index
        integer id_m1(*),ib1_m1(*),ib2_m1(*)
        double precision weight1(*),m1_charge(*)
        character*4 ptnm_m1(*)
        integer m1_said(*)
c
c local
c
        integer poitmp(maxptmo),poitold(maxptmo)
        character*6 name
        integer namel,i,j,jcopy,indx1,indx2,indx3,itmp
        integer itmp1,itmp2,precop,premono,prept
        integer level
        double precision tmp,weight
        data name/'fillpt'/
        data namel/6/

        save poitold

        if (debug) then
         write(stdo,*) ' ipart1 ityp1 kb1 les1 copies '
         write(stdo,*)ipart1,ityp1,kb1,les1,copies
         write(stdo,*) ' id_m1 = ',(id_m1(i),i=1,ipart1)
         write(stdo,*) ' ib1_m1 = ',(ib1_m1(i),i=1,ipart1)
         write(stdo,*) ' ib2_m1 = ',(ib2_m1(i),i=1,ipart1)
        end if

        do 1 i=1,maxptmo
         poitmp(i) = 0
1       continue

c yael - update the monomer division.
      if (mdivyes) then
          mdivlist(0,index) = monodiv(0,ityp1)
        if (mdivlist(0,index).eq.-1) then
           mdivlist(1,index) = monodiv(1,ityp1) + npt
           mdivlist(2,index) = monodiv(2,ityp1) + npt
        endif
        do 9 i=1,monodiv(0,ityp1)+1
            mdivlist(i,index) = monodiv(i,ityp1) + npt
9       continue
      endif

c----------------------------------------------


        if (.not.les1) then
         do 2 i=1,ipart1
          if (debug) then
           write(stdo,*)' id_m1(',i,')    npt ',id_m1(i),npt
          end if
          if (id_m1(i).gt.0) then
           if (debug) then
           write(stdo,*)'basesgm(id_m1(',i,') = ',basesgm(id_m1(i))
           write(stdo,*)'basechg(id_m1(',i,') = ',basechg(id_m1(i))
           write(stdo,*)'baseeps(id_m1(',i,') = ',baseeps(id_m1(i))
           end if
           npt          = npt + 1
           if (npt.gt.maxpt) then
            level = 1
            call alert(name,namel,'Max. number of pt. exceed',25,level)
           end if
           poitmp(i)    = npt
           ptid(npt)    = id_m1(i)
           ptnm(npt)    = ptnm_m1(i)
           ptms(npt)    = basems(id_m1(i))
           ptsaid(npt)  = m1_said(i)
c       ileana -- should be m1_charge and not m_charge
           if (m1_charge(i).lt.998.d0) then
                ptchg(npt) = m1_charge(i)
           else if (basechg(id_m1(i)).lt.998.d0) then
                 ptchg(npt)   = basechg(id_m1(i))
           else
                level = 1
                write(stdo,*)' No charge for ptnm ',ptnm(npt)
                call alert(name,namel,'No charge data',14,level)
           end if
c       ileana

           if(.not.arith) then

           tmp          = basesgm(id_m1(i))*basesgm(id_m1(i))


           else
           tmp          = basesgm(id_m1(i))

           endif

           if (debug) then
            write(stdo,*)' tmp = ',tmp
           end if

c       ileana

           if (.not.arith) then

           epsgm6(npt)  = tmp*tmp*tmp
           tmp          = 4.d0*epsgm6(npt)*baseeps(id_m1(i))

           else

c           epsgm6(npt) = 0.5*tmp
c       because of the combination rule, sigma12 = sigma1 + sigma2, not half this sum

            epsgm6(npt) = tmp  
           tmp = 4.d0*baseeps(id_m1(i))

           endif
              
           if (debug) then
            write(stdo,*)' tmp = ',tmp
           end if

           if(.not.arith) then

           epsgm12(npt) = dsqrt(epsgm6(npt)*tmp)
           epsgm6(npt)  = dsqrt(tmp)
               if (debug) then
               write(stdo,*)' epsgm12 epsgm6 ',epsgm12(npt),epsgm6(npt)
               end if
           else
           epsgm12(npt) = dsqrt(tmp)
              if (debug) then
         write(stdo,*)' sqrt(4*eps)  sigma/2 ',epsgm12(npt),epsgm6(npt)
              end if
           endif


           
           ptwei(npt)   = 1.d0
           lesid(npt)   = 0
           poimon(npt)  = index
          end if
2        continue
         poipt(index ) = npt
         
         if (debug) then
          write(stdo,*)' kb1 = ',kb1
          write(stdo,*)' ib1_m1 ',(ib1_m1(i),i=1,kb1)
          write(stdo,*)' ib2_m1 ',(ib2_m1(i),i=1,kb1)
         end if
         do 3 i=1,kb1
          if (ib1_m1(i).gt.0) then

           if (ib2_m1(i).gt.0) then
c get id of the bonded guys
            indx1   = id_m1(ib1_m1(i))
            indx2   = id_m1(ib2_m1(i))
            if (indx1.eq.0 .or. indx2.eq.0) goto 3
            nb      = nb + 1
            if (nb.gt.maxbond) then
             level = 1
             call alert(name,namel,'Max. # of bonds exceeded',24,level)
            end if
            itmp1   = poitmp(ib1_m1(i))
            itmp2   = poitmp(ib2_m1(i))
            if (itmp1.lt.itmp2) then
             ib1(nb) = itmp1
             ib2(nb) = itmp2
            else
             ib1(nb) = itmp2
             ib2(nb) = itmp1
            end if
c sort indx1 & indx2: indx1 < indx2
            if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
            end if
c Find the index of the present bond in the bond parameter list
c
           if (debug) then
             write(stdo,*)' indx1 indx2 maxubnd bondp(1,1) bondp(1,2) '
             write(stdo,*)indx1,indx2,maxubnd,bondp(1,1),bondp(1,2)
           end if
            indx3   = search(indx1,indx2,0,0,bondp(1,1),
     1              bondp(1,2),0,0,totubnd,2)
            kbond(nb) = basekb(indx3)
            req(nb)   = basereq(indx3)
c if the index is negative it refers to a bond with a particle
c from a previous residue. Use poitold to find out the particle
c number.
c
           else if (ib2_m1(i).lt.0) then
            itmp1     = poitmp(ib1_m1(i))
            itmp2     = poitold(iabs(ib2_m1(i)))
            indx1     = ptid(itmp1)
            indx2     = ptid(itmp2)
c sort indx1 & indx2: indx1 < indx2
            if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
            end if
            indx3     = search(indx1,indx2,0,0,bondp(1,1),
     1              bondp(1,2),0,0,totubnd,2)
c check that the previous monomer was not les type
            if (lesid(itmp2).gt.0) then
c previous monomer was les...
c get the number of copies of previous monomer
                precop = int( 1.d0/ptwei(itmp2) + 0.5d0 )
c get the number of particles in the previous LES monomers
                premono = poimon(itmp2)
                prept   = poipt(premono)-poipt(premono-1)
             else
                precop = 1
                prept = 0
             end if
             tmp = basekb(indx3)*ptwei(itmp2)*ptwei(itmp1)
             do 25 j=1,precop
                 nb = nb + 1
                 ib1(nb)   = itmp2 - (precop - j)*prept
                 ib2(nb)   = itmp1
                 kbond(nb) = tmp
                 req(nb)   = basereq(indx3)
25            continue
           else if (ib2_m1(i).eq.0) then
            level = 1
            call alert(name,namel,'Zero bond pointer',17,level)
           end if

          else if (ib1_m1(i).lt.0) then
           if (ib2_m1(i).lt.0) then
            level = 1
            call alert(name,namel,'Two negative bond pntrs',23,level)
           end if
           itmp1      = poitold(iabs(ib1_m1(i)))
           itmp2      = poitmp(ib2_m1(i))
           indx1      = ptid(itmp1)
           indx2      = ptid(itmp2)
c sort indx1 & indx2: indx1 < indx2
           if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
           end if
           indx3     = search(indx1,indx2,0,0,
     1             bondp(1,1),bondp(1,2),0,0,totubnd,2)
c check that the previous monomer was not LES
           if (lesid(itmp1).gt.0) then
c previous monomer was LES
c Get the number of copies of previous monomer (precop)
                precop = int(1.d0/ptwei(itmp1) + 0.5d0)
c Get the number of particles of previous monomer (prept)
                premono = poimon(itmp1)
                prept   =  poipt(premono) - poipt(premono-1)
           else
c previous monomer was not LES set number of previous copies
c (precop) to 1
                precop = 1
           end if
           tmp = basekb(indx3)*ptwei(itmp1)*ptwei(itmp2)
           do 26 j=1,precop
            nb  = nb + 1
            ib1(nb)    = itmp1 - (precop - j)*prept
            ib2(nb)    = itmp2
            kbond(nb) = tmp
            req(nb)   = basereq(indx3)
26         continue
          else 
           level = 1
           call alert(name,namel,'Zero bond pointer',17,level)
          end if

3        continue

c
c If LES monomer
        else
         if (lestype.gt.0)  then
          lestyp = lestype
         else if (lestype.lt.0) then
          lestyp = lestyp + 1
         end if
c the current implementation is NOT efficient since the whole
c work is simply repated copies time. Will be fixed, sometime...
         do 6 jcopy=1,copies
          if (weight1(1).lt.0) then
           weight = 1.d0/dble(copies)
          else
           weight = weight1(jcopy)
          end if
          do 4 i=1,ipart1
           if (id_m1(i).gt.0) then
            npt          = npt + 1
            if (npt.gt.maxpt) then
             level = 1
             call alert(name,namel,'Max. # of pt. exceed',20,level)
            end if
            poitmp(i)    = npt
            ptid(npt)    = id_m1(i)
            ptnm(npt)    = ptnm_m1(i)
c The mass of the LES particles is multiplied by weight.
c Note that also the non-bonded paramaters are scaled.
c More interactions of LES monomer are likely to be found
c compared to particles in the same residue, therefore all
c variables are scaled down (multiplied by weight)
c Care should be taken when generating non-bonded interactions
c There the interactions between LES particles which belong
c to the same monomer should be DIVIDED by weight
c
            ptms(npt)    = basems(id_m1(i))*weight
            ptchg(npt)   = basechg(id_m1(i))*weight


c       ileana

            if (.not.arith) then
            tmp          = basesgm(id_m1(i))*basesgm(id_m1(i))
            epsgm6(npt)  = tmp*tmp*tmp
            tmp          = 4.d0*epsgm6(npt)*baseeps(id_m1(i))
            epsgm12(npt) = dsqrt(epsgm6(npt)*tmp)*weight
            epsgm6(npt)  = dsqrt(tmp)*weight

            else

            tmp = basesgm(id_m1(i))
            epsgm6(npt) = 0.5d0*tmp
            epsgm12(npt) = dsqrt(baseeps(id_m1(i)))*weight

            endif


            ptwei(npt)   = weight
            lesid(npt)   = lestyp
            poimon(npt)  = index - copies + jcopy
           end if
4         continue

          poipt(index -copies + jcopy ) = npt
         
          do 5 i=1,kb1
           if (ib1_m1(i).gt.0) then

            if (ib2_m1(i).gt.0) then
c get id of the bonded guys
            indx1   = id_m1(ib1_m1(i))
            indx2   = id_m1(ib2_m1(i))
            if (indx1.eq.0 .or. indx2.eq.0) goto 5
             nb      = nb + 1
             if (nb.gt.maxbond) then
              level = 1
              call alert(name,namel,'Max. # of bonds exceeded',24,level)
             end if
             itmp1   = poitmp(ib1_m1(i))
             itmp2   = poitmp(ib2_m1(i))
             if (itmp1.lt.itmp2) then
              ib1(nb) = itmp1
              ib2(nb) = itmp2
             else
              ib1(nb) = itmp2
              ib2(nb) = itmp1
             end if
c sort indx1 & indx2: indx1 < indx2
            if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
            end if
c Find the index of the present bond in the bond parameter list
c
             indx3   = search(indx1,indx2,0,0,
     1              bondp(1,1),bondp(1,2),0,0,totubnd,2)
             kbond(nb) = basekb(indx3)
             if (lesid(ib1(nb)).ne.0 .or. lesid(ib2(nb)).ne.0) then
              if (poimon(ib1(nb)).eq.poimon(ib2(nb))) then
               kbond(nb) = kbond(nb)*ptwei(ib1(nb))
              else
               kbond(nb) = kbond(nb)*ptwei(ib1(nb))*ptwei(ib2(nb))
              end if
             end if
               
             req(nb)   = basereq(indx3)
c if the index is negative it refers to a bond with a particle
c from a previous residue. Use poitold to find out the particle
c number.
c
            else if (ib2_m1(i).lt.0) then
             itmp1   = poitmp(ib1_m1(i))
             itmp2   = poitold(iabs(ib2_m1(i)))
             indx1     = ptid(itmp1)
             indx2     = ptid(itmp2)
c sort indx1 & indx2: indx1 < indx2
             if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
             end if
             indx3     = search(indx1,indx2,0,0,
     1                    bondp(1,1),bondp(1,2),0,0,totubnd,2)
c check if the previous monomer was of LES
             if (lesid(itmp2).gt.0) then
c previous monomer was LES
c Get the number of copies of the previous LES monomer (precop)
                precop = int(1.d0/ptwei(itmp2) + 0.5d0)
C get the number of particles in previous monomer (prept)
                premono = poimon(itmp2)
                prept   = poipt(premono)-poipt(premono-1)
             else
c previous monomer was NOT les type
                precop = 1
             end if
             tmp = basekb(indx3)*ptwei(itmp1)*ptwei(itmp2)
             do 41 j=1,precop
              nb = nb + 1
              ib1(nb) = itmp2 - (precop - j)*prept
              ib2(nb) = itmp1
              kbond(nb) = tmp
              req(nb)   = basereq(indx3)
41           continue
            else if (ib2_m1(i).eq.0) then
             level = 1
             call alert(name,namel,'Zero bond pointer',17,level)
            end if

           else if (ib1_m1(i).lt.0) then
            if (ib2_m1(i).lt.0) then
             level = 1
             call alert(name,namel,'Two negative bond pntrs',23,level)
            end if
             itmp1   = poitold(iabs(ib2_m1(i)))
             itmp2   = poitmp(ib1_m1(i))
             indx1      = ptid(itmp1)
             indx2      = ptid(itmp2)
c sort indx1 & indx2: indx1 < indx2
             if (indx1.gt.indx2) then
                itmp  = indx2
                indx2 = indx1
                indx1 = itmp
             end if
             indx3     = search(indx1,indx2,0,0,
     1                    bondp(1,1),bondp(1,2),0,0,totubnd,2)
c check if previous monomer was of les type
                if (lesid(itmp1).gt.0) then
c previous monomer = LES. Get the number of copies (precop)
                        precop = int(1.d0/ptwei(itmp1)+0.5d0)
c Get the number of particles in previous monomer (prept)
                        premono = poimon(itmp1)
                        prept   = poipt(premono) - poipt(premono-1)
                else
                        prept   = 1
                end if
                tmp = basekb(indx3)*ptwei(itmp1)*ptwei(itmp2)
                do 42 j=1,precop
                 nb = nb + 1
                 ib1(nb) = itmp1 - (precop - j)*prept
                 ib2(nb) = itmp2
                 kbond(nb) = tmp
                 req(nb)   = basereq(indx3)
42              continue

           else 
            level = 1
            call alert(name,namel,'Zero bond pointer',17,level)
           end if

5         continue
6        continue
        end if

c copy poitmp to poitold for next time
c
         do 7 i=1,ipart1
          poitold(i) = poitmp(i)
7        continue

         if (debug) then
          write(stdo,*)' epsgm6(npt) epsgm12(npt)
     1    ',epsgm6(npt),epsgm12(npt)
         end if
         return
         end







