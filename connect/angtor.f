      subroutine angtor(shakm)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/BUFFER.BLOCK'
        include 'COMMON/PROPERT.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
      character*6 name
      integer namel,level
      integer search
      logical shakm
C This subroutine generates lists of all angles, torsions, an
C exclusion list (comprised of all bonds and all 1-3 interactions).
C and a list of special interactions (1-4)
C
C  Written by R.F. Goldstein, UIC, 6-25-91
C  Modified by RE
C
C
C loop variables
      integer ind1,ind2,igrp,l,ix
      integer i, j, k, i1, j1, k1, l1, kkk
      integer itors
      integer temp
C flag
      logical flag
      logical ok1,ok2,ok3
C indices
      integer idtmp(5)
C arrays of all atoms bonded to a given atom
c maxptbn - maximum of particle's bonds set in LENGTH.BLOCK
      integer allbond1(maxptbn)
      integer curbond1, curbond2
      integer tmp_ntors, count_tors
      integer tmp_itor1(maxtors),tmp_itor2(maxtors)
      integer tmp_itor3(maxtors),tmp_itor4(maxtors)

        name = 'angtor'
        namel = 6
        if (nb.eq.0) return
        if (debug) then
                write(stdo,*)' totuagl totutrs totuimp '
                write(stdo,*)totuagl,totutrs,totuimp
        end if
      if (maxbond .lt. nb) then
       level = 1
       call alert(name,namel,'Number of bonds exceeded',24,level)
       return
      end if

c Sort bonds array and set up reverse sorted bonds array in "exclude"
c as a work array. 

      do 1 i=1, nb
c an atom should not be bonded to itself
          if (ib1(i) .eq. ib2(i) ) return
          if (ib1(i) .gt. ib2(i)) then
              temp = ib1(i)
              ib1(i) = ib2(i)
              ib2(i) = temp
          endif
c It is useful to have another copy of bonds sorted according
c to the second atom list, we put them in buffer1 and buffer2
c as temporary arrays.
          buffer1(i) = ib2(i)
          buffer2(i) = ib1(i)
1     continue
      call sort4(nb, ib1, ib2, kbond, req)
      if (debug) then
       write(stdo,*) ' ib1 ib2 kbond req '
       do 18 i=1,nb
        write(stdo,*)ib1(i),ib2(i),kbond(i),req(i)
18     continue
      end if
      call sort2(nb, buffer1, buffer2)
      if (debug) then
       write(stdo,*) 'buffer1 buffer2 '
       do 19 i=1,nb
        write(stdo,*) buffer1(i),buffer2(i)
19     continue
      end if

C Now generate angles.  Consider each atom in turn as a possible central
C atom.  First find all atoms bonded to the central atom and 
C place in "allbond1".
        ind1 = 1
        ind2 = 1
        do 10 i1=ib1(1),buffer1(nb)

c curbond1/2 on input is the number of bonds for which a search
c is required, on output is the number of bonds actually found
c
         curbond1 = nb - ind1 + 1
c the search for bonds is done only from temp up since tha bonds
c were sorted. Also since ib2(i) > ib1(i) bonds on 2 start also
c from temp
         call pickup2(curbond1,i1,ib1(ind1),ib2(ind1),allbond1)
         ind1 = ind1 + curbond1
c         write(stdo,*)'i1   ib1  ib2 '
c          write(stdo,*)i1,ib1(ind1),ib2(ind1)
         
         curbond2 = nb - ind2 + 1
         call pickup2(curbond2,i1,buffer1(ind2),buffer2(ind2)
     1      ,allbond1(curbond1+1))
         curbond1 = curbond1 + curbond2
         if (curbond1.eq.0) go to 10
         ind2 = ind2 + curbond2
C Now form all combinations of atoms in "allbond1" to get angles.
          call piksrt(curbond1, allbond1)
          do 4 j=1, curbond1-1
           do 3 k=j+1, curbond1
             if (nangl .gt. maxangl) then
                level = 1
                call alert(name,namel,'overflow of angles',18,level)
                return
             end if
             j1 = allbond1(j)
             k1 = allbond1(k)
c nangl is initialized to zero in connect, stored in CONNECT.BLOCK
                nangl = nangl + 1
                iangl1(nangl) =  j1
                iangl2(nangl) =  i1
                iangl3(nangl) =  k1
                if (debug) then
                 write(stdo,*)' ANGLE FOUND FOR ',j1,i1,k1
                end if
c pick id's of particles to the present angle
                idtmp(5) = search(ptid(j1),ptid(i1),ptid(k1),0,
     1          anglp(1,1),anglp(1,2),anglp(1,3),0,totuagl,3)
                kangl(nangl) = basekt(idtmp(5))
                angleq(nangl) = baseaeq(idtmp(5))
                if (debug) then
                 write(stdo,*)' nangl = ',nangl
                 write(stdo,*)' basekt baseaeq'
     1           ,basekt(idtmp(5)),baseaeq(idtmp(5))
                 write(stdo,*)' idtmp(5) kangl',idtmp(5),kangl(nangl)
                 write(stdo,*)' ptwei(i,j,k)=',
     1           ptwei(i1),ptwei(j1),ptwei(k1)
                end if
3          continue
4         continue
c
c here we can also do the improper torsions. If i1 has only
c 3 bonds then an improper torsion should be set
c If les particles present then then improper torsion may exist
c of the type i'- i'- r- r -i -i -i -i where all the i's are of the
c same lesid. Waht is done below is to check how many group of
c different lesid are. If there are three, all possible combinations of
c improper torsions BETWEEN THE GROUPS must be included.
c

c if curbond1.lt.3 improper torsion not possible
c
        if (debug) write(stdo,*)' ptnm = '
     1  ,ptnm(i1),ptnm(i1)(1:1)
        if (curbond1.ne.3) go to 10
         j1 = allbond1(1)
         k1 = allbond1(2)
         l1 = allbond1(3)
         idtmp(1) = ptid(i1)
         idtmp(2) = ptid(j1)
         idtmp(3) = ptid(k1)
         idtmp(4) = ptid(l1)
         call srch_imp(idtmp(1),idtmp(2),idtmp(3),idtmp(4),
     1     i1,j1,k1,l1,improp(1,1),improp(1,2),improp(1,3),
     2          improp(1,4),totuimp,4,idtmp(5))
         nimp = nimp + 1
         if (nimp.gt.maximp) then
          level = 1
          write(*,*)' maximp = ',maximp,'   nimp = ',nimp
          call alert(name,namel,'maximp too small',15,level)
          return
         end if
         iimp1(nimp) = i1
         iimp2(nimp) = j1
         iimp3(nimp) = k1
         iimp4(nimp) = l1
         impeq(nimp) = baseieq(idtmp(5))
         kimp(nimp)  = baseimp(idtmp(5))
10      continue

C Now generate torsions. Torsions are extracted from the angle list.
C Given a list of angles 1-2-3, compare all angles say 1-2-3
C and 1'-2'-3' for 1-2 = 2'-3' or for 2-3 = 1'-2'. For such
C equalities generate torsions 1'-1-2-3 or 1-1'-2'-3'. Note that
C LES restricitions should be already statisfied in the angle level
C and they do not require further checks here. We assume that there
c are no dubplicate of angles.

c Note also that after the torsions are generated we have to take
c care of another problem: for multiple rotations around a bond
c the force constant need to be scaled down to its original value.
c i.e. two torsions 1-2-3-4 & 5-2-3-6 have half the amplitude of
c the case that only 1-2-3-4 exists. We cannot just count central
c bonds since LES torsions considered as individual torsions alrady
c normalised. We ASSUME that this problem was taken care in the
c setup of the parameter list. This may create difficulties in
c transferring parameters from other force fields.

c If nangl.le.1 no torsions possible

        if (nangl.le.1) go to 12
        do 11 i=1,nangl
         do 11 j=i+1,nangl
          if ( (iangl1(i).eq.iangl2(j)) .and. 
     1          (iangl2(i).eq.iangl3(j)) .and.
     2          (iangl1(j).ne.iangl3(i))) then
                i1 = iangl1(j)
                j1 = iangl2(j)
                k1 = iangl3(j)
                l1 = iangl3(i)
                if (l1.lt.i1) then
                 kkk = l1
                 l1  = i1
                 i1  = kkk
                 kkk = j1
                 j1  = k1
                 k1  = kkk
                end if
                idtmp(1) = ptid(i1)
                idtmp(2) = ptid(j1)
                idtmp(3) = ptid(k1)
                idtmp(4) = ptid(l1)
                idtmp(5) = search(idtmp(1),idtmp(2),idtmp(3),idtmp(4),
     1                  torsp(1,1),torsp(1,2),torsp(1,3),torsp(1,4),
     2                  totutrs,4)
 1001           ntors = ntors + 1
                itor1(ntors) = i1
                itor2(ntors) = j1
                itor3(ntors) = k1
                itor4(ntors) = l1
                if(idtmp(5).gt.0) then
                phase1(ntors)  = basephs1(idtmp(5))
                phase2(ntors)  = basephs2(idtmp(5))
                phase3(ntors)  = basephs3(idtmp(5))
                period(ntors) = baseper(idtmp(5))
                ktors1(ntors) = basamp1(idtmp(5))
                ktors2(ntors) = basamp2(idtmp(5))
                ktors3(ntors) = basamp3(idtmp(5))

               
                else
                phase1(ntors) = 0
                phase2(ntors) = 0
                phase3(ntors) = 0
                period(ntors)= 0
                ktors1(ntors)= 0
                ktors2(ntors) = basamp2(idtmp(5))
                ktors3(ntors) = basamp3(idtmp(5))
                end if

c     ileana
                if(arith) then
                if(baseper(idtmp(5)).lt.0) then
c                   baseper(idtmp(5)) = -baseper(idtmp(5))
                    period(ntors) = -baseper(idtmp(5))
                   idtmp(5) = idtmp(5) +1
                   
                   goto 1001
                endif
                endif  
 

          else if ( (iangl2(i).eq.iangl1(j)) .and. 
     1          (iangl3(i).eq.iangl2(j)) .and.
     2          (iangl1(i).ne.iangl3(j))) then
                i1 = iangl1(i)
                j1 = iangl2(i)
                k1 = iangl3(i)
                l1 = iangl3(j)
                if (l1.lt.i1) then
                 kkk = l1
                 l1  = i1
                 i1  = kkk
                 kkk = j1
                 j1  = k1
                 k1  = kkk
                end if
                idtmp(1) = ptid(i1)
                idtmp(2) = ptid(j1)
                idtmp(3) = ptid(k1)
                idtmp(4) = ptid(l1)
                
                idtmp(5) = search(idtmp(1),idtmp(2),idtmp(3),idtmp(4),
     1                  torsp(1,1),torsp(1,2),torsp(1,3),torsp(1,4),
     2                  totutrs,4)
c                 write(stdo,*)'idtmp5 is ',idtmp(5)

 2001           ntors = ntors + 1
                itor1(ntors) = i1
                itor2(ntors) = j1
                itor3(ntors) = k1
                itor4(ntors) = l1
                if (idtmp(5).gt.0) then
                period(ntors) = baseper(idtmp(5))
                ktors1(ntors) = basamp1(idtmp(5))
                ktors2(ntors) = basamp2(idtmp(5))
                ktors3(ntors) = basamp3(idtmp(5))
                phase1(ntors)  = basephs1(idtmp(5))
                phase2(ntors)  = basephs2(idtmp(5))
                phase3(ntors)  = basephs3(idtmp(5))
                else
                period(ntors) = 0
                ktors1(ntors) = 0
                ktors2(ntors) = 0
                ktors3(ntors) = 0
                phase1(ntors)  = 0
                phase2(ntors)  = 0
                phase3(ntors)  = 0
                end if
                if (debug) then
                 write (stdo,*)' id phase1 ',idtmp(5),phase1(ntors)
                 write (stdo,*)' id phase2 ',idtmp(5),phase2(ntors)
                 write (stdo,*)' id phase3 ',idtmp(5),phase3(ntors)
                end if

c     ileana
                if(arith) then
                if(baseper(idtmp(5)).lt.0) then

c                    write(stdo,*)'found one negative periodicity'
c                   baseper(idtmp(5)) = -baseper(idtmp(5))
                    period(ntors) = -baseper(idtmp(5))
                   idtmp(5) = idtmp(5) +1
                  
                   goto 2001
                endif
                endif

          else if ( (iangl1(i).eq.iangl2(j)) .and. 
     1          (iangl2(i).eq.iangl1(j)) .and.
     2          (iangl3(i).ne.iangl3(j))) then
                i1 = iangl3(i)
                j1 = iangl2(i)
                k1 = iangl1(i)
                l1 = iangl3(j)
                if (l1.lt.i1) then
                 kkk = l1
                 l1  = i1
                 i1  = kkk
                 kkk = j1
                 j1  = k1
                 k1  = kkk
                end if
                idtmp(1) = ptid(i1)
                idtmp(2) = ptid(j1)
                idtmp(3) = ptid(k1)
                idtmp(4) = ptid(l1)
                idtmp(5) = search(idtmp(1),idtmp(2),idtmp(3),idtmp(4),
     1                  torsp(1,1),torsp(1,2),torsp(1,3),torsp(1,4),
     2                  totutrs,4)
c                 write(stdo,*)'idtmp5 is ',idtmp(5)

 3001           ntors = ntors + 1
                itor1(ntors) = i1
                itor2(ntors) = j1
                itor3(ntors) = k1
                itor4(ntors) = l1
                if (idtmp(5).gt.0) then
                period(ntors) = baseper(idtmp(5))
                ktors1(ntors) = basamp1(idtmp(5))
                ktors2(ntors) = basamp2(idtmp(5))
                ktors3(ntors) = basamp3(idtmp(5))
                phase1(ntors)  = basephs1(idtmp(5))
                phase2(ntors)  = basephs2(idtmp(5))
                phase3(ntors)  = basephs3(idtmp(5))
                else
                period(ntors) = 0
                ktors1(ntors) = 0
                ktors2(ntors) = 0
                ktors3(ntors) = 0
                phase1(ntors)  = 0
                phase2(ntors)  = 0
                phase3(ntors)  = 0
                end if
                if (debug) then
                 write (stdo,*)' id phase1 ',idtmp(5),phase1(ntors)
                 write (stdo,*)' id phase2 ',idtmp(5),phase2(ntors)
                 write (stdo,*)' id phase3 ',idtmp(5),phase3(ntors)
                end if

c     ileana

                if(arith) then
                if(baseper(idtmp(5)).lt.0) then

c                    write(stdo,*)'found one negative periodicity'
c                   baseper(idtmp(5)) = -baseper(idtmp(5))
                    period(ntors) = -baseper(idtmp(5))
                   idtmp(5) = idtmp(5) +1
                  
                   goto 3001
                endif
                endif

          else if ( (iangl2(i).eq.iangl3(j)) .and. 
     1          (iangl3(i).eq.iangl2(j)) .and.
     2          (iangl1(i).ne.iangl1(j))) then
                i1 = iangl1(i)
                j1 = iangl2(i)
                k1 = iangl3(i)
                l1 = iangl1(j)
                if (l1.lt.i1) then
                 kkk = l1
                 l1  = i1
                 i1  = kkk
                 kkk = j1
                 j1  = k1
                 k1  = kkk
                end if
                idtmp(1) = ptid(i1)
                idtmp(2) = ptid(j1)
                idtmp(3) = ptid(k1)
                idtmp(4) = ptid(l1)
                idtmp(5) = search(idtmp(1),idtmp(2),idtmp(3),idtmp(4),
     1                  torsp(1,1),torsp(1,2),torsp(1,3),torsp(1,4),
     2                  totutrs,4)
c                 write(stdo,*)'idtmp5 is ',idtmp(5)


 4001           ntors = ntors + 1
                itor1(ntors) = i1
                itor2(ntors) = j1
                itor3(ntors) = k1
                itor4(ntors) = l1
                if (idtmp(5).gt.0) then
                period(ntors) = baseper(idtmp(5))
                ktors1(ntors)  = basamp1(idtmp(5))
                ktors2(ntors)  = basamp2(idtmp(5))
                ktors3(ntors)  = basamp3(idtmp(5))
                phase1(ntors)  = basephs1(idtmp(5))
                phase2(ntors)  = basephs2(idtmp(5))
                phase3(ntors)  = basephs3(idtmp(5))
                else
                period(ntors) = 0
                ktors1(ntors)  =0
                ktors2(ntors)  =0
                ktors3(ntors)  =0
                phase1(ntors)  = 0
                phase2(ntors)  = 0
                phase3(ntors)  = 0
                end if
                if (debug) then
                 write (stdo,*)' id phase1 ',idtmp(5),phase3(ntors)
                 write (stdo,*)' id phase2 ',idtmp(5),phase2(ntors)
                 write (stdo,*)' id phase3 ',idtmp(5),phase1(ntors)
                end if

c     ileana

                if(arith) then
                if(baseper(idtmp(5)).lt.0) then

c                    write(stdo,*)'found one negative periodicity'
c                   baseper(idtmp(5)) = -baseper(idtmp(5))
                   period(ntors) = -baseper(idtmp(5))
                   idtmp(5) = idtmp(5) +1
                  
                   goto 4001
                endif
                endif

           end if
11      continue
12      continue

C Finally, set up the exclusion list. We assume no repetition
C of bonds or angles.
C
C VARIATION - FEB 7, 1995, The exclusion is extended to include both -
C               the 1,2,3 lists and the 1-4 list
C
C Do a loop over all particles, bonds and angles
C This is rather inefficient setup. Probably much fatser to sort 
C first the bond and angle list. Kept this way to keep the
C clarity of the code and noticing too, that this is never time
C bottleneck in the calculations.
C ix - the number of exclusion found for current atom

C the exclusions are set with the assumption that the bonds
C are sorted in the first Columm and the angle are sorted in the
c second and then in the first columm
        ix = 0
        do 15 i=1,npt
         do 13 k=1,nb
          if (ib1(k).eq.i) then
           ix = ix + 1
           exc2(ix) = ib2(k)
          end if
13       continue
         do 14 k=1,nangl
          if (iangl1(k).eq.i) then
           ix = ix + 1
           exc2(ix) = iangl3(k)
          end if
14       continue
c
c HERE THE TORSIONS ARE ADDED TO THE EXCLUSION LIST
c

         if(.not.arith) then

         do 141 itors=1,ntors 
                if (itor1(itors).eq.i) then
                        if (i.lt.itor4(itors)) then
                         ix = ix + 1
                         exc2(ix) = itor4(itors)
                        end if
                else if (itor4(itors).eq.i) then
                        if (i.lt.itor1(itors)) then
                         ix = ix + 1
                         exc2(ix) = itor1(itors)
                        end if
              endif            
141      continue  
c         write(stdo,*)'current ix is ', ix
         




        else

c  ileana  -need to eliminate 4 atom-sets that appear 
C  twice due to the Amber torsion list
          
           
        if (ntors .gt. 0) then
           tmp_ntors = 1
           tmp_itor1(tmp_ntors) = itor1(1)
           tmp_itor2(tmp_ntors) = itor2(1)
           tmp_itor3(tmp_ntors) = itor3(1)
           tmp_itor4(tmp_ntors) = itor4(1)

           do 1141 itors=2,ntors

         
            if(.not. (itor1(itors).eq.itor1(itors-1).and.
     1                  itor2(itors).eq.itor2(itors-1).and.
     2                  itor3(itors).eq.itor3(itors-1).and.
     3                  itor4(itors).eq.itor4(itors-1)) ) then
                tmp_ntors = tmp_ntors + 1
                tmp_itor1(tmp_ntors) = itor1(itors)
                tmp_itor2(tmp_ntors) = itor2(itors)
                tmp_itor3(tmp_ntors) = itor3(itors)
                tmp_itor4(tmp_ntors) = itor4(itors)
            end if

 1141       continue
        end if
            




c     -now generate exclusions from the temporary itor vector

             do 1241 itors=1,tmp_ntors 
                if (tmp_itor1(itors).eq.i) then
                        if (i.lt.tmp_itor4(itors)) then
                         ix = ix + 1
                         exc2(ix) = tmp_itor4(itors)
                        end if
                else if (tmp_itor4(itors).eq.i) then
                        if (i.lt.tmp_itor1(itors)) then
                         ix = ix + 1
                         exc2(ix) = tmp_itor1(itors)
                        end if
              endif            
1241     continue  

      endif
   
      
         exc1(i) = ix
15      continue
        totex = ix
        write(stdo,*)'total number of exclusions is ',totex

         

        if (totex.gt.maxex) then
         write(stdo,*)' intermediate totex maxex ',totex,maxex
         call alert(name,namel,'Max excl exceeded',17,1)
         endif
c
C check that there are no duplications in the exclusion list
c
        do 16 i=1,npt
         if ( exc1(i-1)+1.ge.exc1(i) ) go to 16
C@
          j = exc1(i-1)
151       j = j + 1
          k  = j
152       k = k + 1
C@
C       if (i.eq.45) then
C         write(*,*)' j k  ',j,k,
C     1    ' exc2(j) exc2(k) ',exc2(j),exc2(k)
C       write(*,*)(exc2(i1),i1=exc1(i-1)+1,exc1(i))
C       end if
          if (k.le.exc1(i)) then
                if (exc2(j).eq.exc2(k)) then
C@
C       write(*,*)' HIT. removing k = ',k
C       write(*,*)' j k i exc2(k) ',j,k,i,exc2(k)
C       write(*,*)' before '
C       write(*,*)(exc2(i1),i1=exc1(i-1)+1,exc1(i))
                  call rm_elemi(exc2,k,ix)
                  ix = ix - 1
                  do 153 i1=i,npt
                   exc1(i1) = exc1(i1)-1
153               continue
C       write(*,*)' AFTER '
C       write(*,*)(exc2(i1),i1=exc1(i-1)+1,exc1(i))
                  j  = j  - 1
                  go to 151
                end if
                go to 152
          end if
          if (j.lt.exc1(i)-1) go to 151
16       continue

        totex = ix
C@
C       write(*,*)' totex = ',totex
C Do similar thing for the special 1-4 interactions

        if (debug) write(stdo,*)' before 1-4 list '
C if *** NO TORSIONS *** MOST WHAT IS BELOW CAN BE ELIMINATED
c
        if (ntors.eq.0) then
          totspe = 0

C some torsions exist and therefore special list is required
        else


            if(.not.arith) then

               totspe = 0
               do 17 i=1,npt
                  do 17 j=1,ntors
                     if (itor1(j).eq.i ) then
                        if (.not. (
     1      (epsgm12(i).lt.1.d-3 .or. epsgm12(itor4(j)).lt.1.d-3)
     2     .and.
     3      (dabs(ptchg(i)*ptchg(itor4(j))).lt.1.d-5)
     4           ) ) then
            totspe = totspe + 1
            if (totspe.gt.maxspec) then
               write(stdo,*) ' Current 1-4 interactions ',ix
               write(stdo,*) ' Maximum allowed ',maxspec
               call alert(name,namel,'Maxspec exceeded',16,1)
            end if
            spec1(totspe) = i
            spec2(totspe) = itor4(j)
          end if
         end if
17      continue



       

             else

                
                
c     ileana-- arith flag is true here

                totspe = 0
               do 1117 i=1,npt
                  do 1117 j=1,tmp_ntors
                     if (tmp_itor1(j).eq.i ) then
                        if (.not. (
     1      (epsgm6(i).lt.0.25d0.or.epsgm6(tmp_itor4(j)).lt.0.25d0)
     2     .and.
     3      (dabs(ptchg(i)*ptchg(tmp_itor4(j))).lt.1.d-5)
     4           ) ) then
            totspe = totspe + 1
            if (totspe.gt.maxspec) then
               write(stdo,*) ' Current 1-4 interactions ',ix
               write(stdo,*) ' Maximum allowed ',maxspec
               call alert(name,namel,'Maxspec exceeded',16,1)
            end if
            spec1(totspe) = i
            spec2(totspe) = tmp_itor4(j)
          end if
         end if
1117    continue


        endif





c compress list of special interactions to avoid duplications
c if a 1-4 appears twice. There are two options: For OPLS BOTH are removed.
c For AMBER only one duplicate is removed. 
c RE. Aug 28 02
       
        if (totspe.eq.0) go to 1998
        i = 0
188     i = i + 1
        j = i
198     j = j + 1
          if (spec1(i).eq.spec1(j) .and. spec2(i).eq.spec2(j)) then
                if (.not.arith)then
                  call rm_elemi(spec1,j,totspe)
                  call rm_elemi(spec2,j,totspe)
                  j = j - 1
                  totspe = totspe - 1
                end if
                  call rm_elemi(spec1,i,totspe)
                  call rm_elemi(spec2,i,totspe)
                  i = i - 1
                  totspe = totspe - 1
          end if

          
          if (totspe.eq.0) go to 1998
          if (i.eq.0) go to 188
          
          
          
          if (j.lt.totspe) go to 198
          if (i.lt.totspe-1) go to 188
         
1998      continue
c
c if 1-4 interaction appears also as 1-3 in the angle list it
c corresponds to a 5 member ring and should be eliminated
c
        do 199 i=1,nangl
         j = 0
200      j = j + 1
         if (j.gt.totspe) go to 199
         if (iangl1(i).eq.spec1(j).and.iangl3(i).eq.spec2(j)) then
                call rm_elemi(spec1,j,totspe)
                call rm_elemi(spec2,j,totspe)
                totspe = totspe - 1
                j      = j - 1
         end if
         go to 200
199     continue
        

       

           
        

c     compressing torsion list; torsions with 
C      k1=k2=k3=0 will be eliminated
        
        if (ntors.ne.0) then
        i = 1
20      continue
        if ( (.not.arith.and. ktors1(i).eq.0.0 .and. ktors2(i).eq.0.0
     1        .and. ktors3(i).eq.0.0)
     2     .or. (arith.and.ktors2(i).eq.0.0) ) then
           ntors = ntors - 1
           if (ntors.eq.0) go to 23
           do 21 j=i,ntors
            itor1(j)  = itor1(j+1)
            itor2(j)  = itor2(j+1)
            itor3(j)  = itor3(j+1)
            itor4(j)  = itor4(j+1)
            ktors1(j) = ktors1(j+1)
            ktors2(j) = ktors2(j+1)
            ktors3(j) = ktors3(j+1)
            phase1(j)  = phase1(j+1)
            phase2(j)  = phase2(j+1)
            phase3(j)  = phase3(j+1)
            period(j) = period(j+1)
21         continue
           i = i - 1
          end if
          i = i + 1
         if  (i.le.ntors) go to 20

         if (debug) then
          write(stdo,*) ' List of remaining torsions '
          do 22 j=1,ntors
           write(stdo,*)itor1(j),itor2(j),itor3(j),itor4(j),
     1          ktors1(j),ktors2(j),ktors3(j),period(j),
     2          phase1(j),phase2(j),phase3(j)
22        continue
         end if
        end if

c END IF FOR *** NO TORSIONS *** case
        end if

23      continue


c   compresing improper torsions
c   improper torsion with kimp=0 will be eliminated
        if (nimp.ne.0) then
        i = 1
120      continue
        if ( kimp(i).eq.0.0) then
           nimp = nimp - 1
           if (nimp.eq.0) go to 123
           do j=i,nimp
            iimp1(j)  = iimp1(j+1)
            iimp2(j)  = iimp2(j+1)
            iimp3(j)  = iimp3(j+1)
            iimp4(j)  = iimp4(j+1)
            kimp(j)   = kimp(j+1)
            impeq(j)  = impeq(j+1)
           end do
           i = i - 1
          end if
          i = i + 1
         if  (i.le.nimp) go to 120

        end if  ! nimp==0
        
123      continue         

c also compress the angle list. TIP3 molecules should not have any
c angle. This is done only if shakm is planned (shakm keyword
c as an input to con is required).
c
        if (.not. shakm) return

        if (nangl.ne.0) then
        i = 1
24      continue
         if (moname(poimon(iangl1(i))).eq.'TIP3' .or.
     1      moname(poimon(iangl1(i)))(1:3).eq.'SPC') then
          nangl = nangl - 1
          if (nangl.eq.0) go to 26
          do 25 j=i,nangl
                iangl1(j) = iangl1(j+1)
                iangl2(j) = iangl2(j+1)
                iangl3(j) = iangl3(j+1)
                kangl(j)  = kangl(j+1)
                angleq(j) = angleq(j+1)
25        continue
          i = i-1
         end if
         i = i + 1
         if (i.le.nangl) go to 24
        end if
26      continue

c compress also bonds in TIP3 monomers
c
        if (nb.ne.0) then
         i = 1
27       continue
         if (moname(poimon(ib1(i))).eq.'TIP3' .or.
     1         moname(poimon(ib1(i)))(1:3).eq.'SPC') then
          nb = nb -1
          if (nb.eq.0) return
          do 28 j=i,nb
                ib1(j) = ib1(j+1)
                ib2(j) = ib2(j+1)
                kbond(j) = kbond(j+1)
                req(j)   = req(j+1)
28        continue
          i = i-1
         end if
         i = i + 1
         if (i.le.nb) go to 27
        end if

      return
      end


C********************************************
      subroutine pickup2(n, atom, arr1, arr2, out)
C********************************************
      integer n, i, atom, arr1(n), arr2(n), out(n)
C starting from arr1(1), finds all consecutive entries in that match
C atom, and places the corresponding entries of arr2() into out ()
C returns the number moved into out in the parameter "n"
C
      do 10 i=1,n
          if (arr1(i) .eq. atom) then
              out(i) = arr2(i)
          else
              n = i-1
              goto 20
          endif
10    continue
      return
20    continue
      return
      end
 
 
C********************************************
      subroutine sort1(n, ra)
C********************************************
C Sorts an array RA of length N into ascending numerical order using
C the Heapsort algorithm, while makeing the corresponding rearrangement
C of the array RB.  From Numerical Recipes, Press et al.
      integer n
      integer ra(n)
      integer rra,l,ir,i,j

      l = n/2 + 1
      ir = n
10    continue
          if (l .gt. 1) then
              l = l-1
              rra = ra(l)
          else
              rra = ra(ir)
              ra(ir) = ra(1)
              ir = ir-1
              if (ir .eq. 1) then
                  ra(1) = rra
                  return
               endif
           endif
           i = l
           j = l+l
20         if (j .le. ir) then
               if (j .lt. ir) then
                   if (ra(j) .lt. ra(j+1)) j = j+1
               endif
               if (rra .lt. ra(j)) then
                   ra(i) = ra(j)
                   i=j
                   j=j+j
               else
                   j = ir+1
               endif
           goto 20
           endif
           ra(i) = rra
       goto 10
       end            

        subroutine sort4(n,ra,rb,rc,rd)
        integer n,i,j
        integer ra(*),rb(*)
        double precision rc(*),rd(*)
        integer rra,rrb
        double precision rrc,rrd
        if (n.eq.1) return
        do 12 j=2,n
         rra = ra(j)
         rrb = rb(j)
         rrc = rc(j)
         rrd = rd(j)
         do 11 i=j-1,1,-1
          if(ra(i).le.rra) go to 10
          ra(i+1) = ra(i)
          rb(i+1) = rb(i)
          rc(i+1) = rc(i)
          rd(i+1) = rd(i)
11       continue
         i=0
10       ra(i+1) = rra
         rb(i+1) = rrb
         rc(i+1) = rrc
         rd(i+1) = rrd
12       continue
        return
        end

        subroutine sort2(n,ra,rb)
        integer n
        integer ra(*),rb(*)
        integer rra,rrb
        integer i,j
        if (n.eq.1) return
        do 12 j=2,n
         rra = ra(j)
         rrb = rb(j)
         do 11 i=j-1,1,-1
          if(ra(i).le.rra) go to 10
          ra(i+1) = ra(i)
          rb(i+1) = rb(i)
11       continue
         i=0
10       ra(i+1) = rra
         rb(i+1) = rrb
12       continue
        return
        end
 
 
C********************************************
      subroutine sort2a(n, ra, rb)
C********************************************
C Sorts an array RA of length N into ascending numerical order using
C the Heapsort algorithm, while makeing the corresponding rearrangement
C of the array RB.  From Numerical Recipes, Press et al.
      integer n,l,ir,i,j
      integer ra(n), rb(n)
      integer rra, rrb

      l = n/2 + 1
      ir = n
10    continue
          if (l .gt. 1) then
              l = l-1
              rra = ra(l)
              rrb = rb(l)
          else
              rra = ra(ir)
              rrb = rb(ir)
              ra(ir) = ra(1)
              rb(ir) = rb(1)
              ir = ir-1
              if (ir .eq. 1) then
                  ra(1) = rra
                  rb(1) = rrb
                  return
               endif
           endif
           i = l
           j = l+l
20         if (j .le. ir) then
               if (j .lt. ir) then
                   if (ra(j) .lt. ra(j+1)) j = j+1
               endif
               if (rra .lt. ra(j)) then
                   ra(i) = ra(j)
                   rb(i) = rb(j)
                   i=j
                   j=j+j
               else
                   j = ir+1
               endif
           goto 20
           endif
           ra(i) = rra
           rb(i) = rrb
       goto 10
       end            

 
C********************************************
      subroutine sort4a(n, ra, rb, rc, rd)
C********************************************
C Sorts an array RA of length N into ascending numerical order using
C the Heapsort algorithm, while makeing the corresponding rearrangement
C of the arrays RB,RC,RD.  From Numerical Recipes, Press et al.
      integer n,l,ir,i,j
      integer ra(n), rb(n)
      double precision rc(n),rd(n)
      integer rra, rrb 
      double precision rrc,rrd

      if (n.eq.2) then
       if (ra(1).lt.ra(2)) return
       rra = ra(1)
       ra(1) = ra(2)
       ra(2) = rra
       rrb = rb(1)
       rb(1) = rb(2)
       rb(2) = rrb
       rrc = rc(1)
       rc(1) = rc(2)
       rc(2) = rrc
       rrd = rd(1)
       rd(1) =rd(2)
       rd(2) = rrd
       return
      end if
      l = n/2 + 1
      ir = n
10    continue
          if (l .gt. 1) then
              l = l-1
              rra = ra(l)
              rrb = rb(l)
              rrc = rc(l)
              rrd = rd(l)
          else
              rra = ra(ir)
              rrb = rb(ir)
              rrc = rc(ir)
              rrd = rd(ir)
              ra(ir) = ra(1)
              rb(ir) = rb(1)
              rc(ir) = rc(1)
              rd(ir) = rd(1)
              ir = ir-1
              if (ir .eq. 1) then
                  ra(1) = rra
                  rb(1) = rrb
                  rc(1) = rrc
                  rd(1) = rrd
                  return
               endif
           endif
           i = l
           j = l+l
20         if (j .le. ir) then
               if (j .lt. ir) then
                   if (ra(j) .lt. ra(j+1)) j = j+1
               endif
               if (rra .lt. ra(j)) then
                   ra(i) = ra(j)
                   rb(i) = rb(j)
                   rc(i) = rc(j)
                   rd(i) = rd(j)
                   i=j
                   j=j+j
               else
                   j = ir+1
               endif
           goto 20
           endif
           ra(i) = rra
           rb(i) = rrb
           rc(i) = rrc
           rd(i) = rrd
       goto 10
       end            

C********************************************
      subroutine piksrt(n, arr)
C********************************************
C Sorts an array "arr" of length n into ascending numerical order,
C by straight insertion.  n is input; arr is replaced on output
C by its sorted rearrangement.  From Numerical Recipes, Press et al.
      integer n, arr(n)
      integer a,i,j

      do 12 j=2,n
          a = arr(j)
          do 11 i= j-1, 1, -1
              if (arr(i) .le. a) goto 10
              arr(i+1) = arr(i)
11        continue
          i=0
10        arr(i+1) = a
12    continue
      return
      end










