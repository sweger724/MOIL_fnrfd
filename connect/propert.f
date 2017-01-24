        subroutine prprt(uprm,basenm,basems,basechg,
     1          baseeps,basesgm,maxunq,totunq,
     2          bondp,basekb,basereq,maxubnd,totubnd,
     3          anglp,basekt,baseaeq,maxuagl,totuagl,
     4          torsp,basamp1,basamp2,basamp3,baseper,
     5          basephs1,basephs2,basephs3,maxutrs,totutrs,
     5          improp,baseimp,baseieq,maxuimp,totuimp,arith,el14,v14)
      implicit none
c
c reading particle properties
c
c basenm  : char*4  - particle names
c basems  : double  - particle masses
c basechg : double  - particle charges
c baseeps : double  - particle van der Waals (vdw) binding energy
c basesgm : double  - particle vdw radius
c maxunq  : integer - maximum number of unique particles (input)
c totunq  : integer - actual number of unique atom found (output)
c bondp   : integer - pointer to bonds between unique particles
c       bondp(i,1) bondp(i,2) - first and second atoms of i-th bond
c        1-2
c basekb  : double  - force constant of bonds
c basereq : double  - equilibrium distances of bonds
c maxubnd : integer - maximum number of unique bonds (input)
c totubnd : integer - actual number of unique bonds found (output)
c anglp   : integer - pointer to angles between unique particles
c       anglp(i,1-3) first to third atoms in the i-th angle
c        1-2-3
c basekt  : double  - force constants for angles
c baseaeq : double  - equilibrium positions for angles
c maxuagl : integer - maximum number of unique angles (input)
c totuagl : integer - actual number of unique angles found (output)
c torsp   : integer - pointer to torsion between unique particles
c       torsp(i,1-4) first to fourth atoms in the i-th torsion
c        1-2-3-4
c basamp[1-3] : double  - amplitude for unique torsions
c baseper : integer - period of torsional angles
c basephs[1-3]: double  - phase for unique torsions
c maxutrs : integer - maximum number of unique torsions (input)
c totutrs : integer - actual number of unique torsions found (output)
c improp  : integer - pointer to imp. torsions between unq prtcles
c       improp(i,1-4) first to fourth atom for the i-th imp.torsion
c          2
c          |
c        4-1-3
c baseimp : double  - amplitude for improper torsions
c baseieq : double  - equilibrium position of imp. torsions
c maxuimp : integer - maximum number on unique improper torsion (input)
c totuimp : integer - actual number of unq. imp. torsions (output)
c
        integer uprm,maxunq,totunq,maxubnd,totubnd,maxuagl,totuagl
        integer maxutrs,totutrs,maxuimp,totuimp
        integer bondp(maxubnd,2),anglp(maxuagl,3),torsp(maxutrs,4)
        integer improp(maxuimp,4),baseper(maxutrs)

        character*4 basenm(maxunq)

        double precision basems(maxunq),basechg(maxunq),baseeps(maxunq)
        double precision basesgm(maxunq)
        double precision basekb(maxubnd),basereq(maxubnd)
        double precision basekt(maxuagl),baseaeq(maxuagl)
        double precision basamp1(maxutrs),basephs1(maxutrs)
        double precision basephs2(maxutrs),basephs3(maxutrs)
        double precision basamp2(maxutrs),basamp3(maxutrs)
        double precision baseimp(maxuimp),baseieq(maxuimp)
        double precision el14,v14

c
c local
c
        character*7 name
        character*5 getchar
        character*4 char1,char2,char3,char4
        integer intgr
        integer namel,coun,i,j,lc,level
        double precision getd
        logical find, success, empty,arith
        double precision number

        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/LINE.BLOCK'

        name = 'propert'
        namel = 7
        lc    = 4

        call init_var()
        el14 = -1.d0
        v14 = -1.d0

c Read current input line. Line can be extended into more than one fortran
c read by adding "-" to the end of the line.
c
        silent = .true.
c@
c@	silent = .false.
c@
        call rline(name,namel,uprm)

c first exe line must include PRTC
c properties of individual particles
c
        if (.not.find('PRTC')) then
          level = 1
          call alert(name,namel,'Missing PRTC header',19,level)
        end if


c
c loop for reading particles properties
c
        do 1 coun=1,maxunq
         call rline(name,namel,uprm)
         if (find('DONE')) then
          totunq = coun - 1
          write(stdo,101)
101       format(//,1x,' *** DONE WITH PARTICLE PROPERTIES *** ',//)
          go to 2
         end if
         basenm(coun)  = getchar('PNAM','????',lc)
         basems(coun)  = (getd('PMAS',-1.d0))
         basechg(coun) = (getd('PCHG',999.d0))
         baseeps(coun) = (getd('PEPS',-1.d0))
         basesgm(coun) = (getd('PSGM',-1.d0))
         if (debug) then
         write(stdo,*)' coun = ',coun
         write(stdo,*)' name mass chg eps sgm '
         write(stdo,*)basenm(coun),basems(coun),basechg(coun)
     1   ,baseeps(coun),basesgm(coun)
         end if
         if (basenm(coun).eq.'????') then
                level = 1
                call alert(name,namel,'No prtc name',12,level)
         else if (basems(coun).lt.0) then
                level = 1
                call alert(name,namel,'No prtc mass',12,level)
         else if (basechg(coun).gt.998.) then
                level = 1
                call alert(name,namel,'No prtc charge',14,level)
         else if (baseeps(coun).lt.0) then
                level = 1
                call alert(name,namel,'No prtc epsilon',15,level)
         else if (basesgm(coun).lt.0) then
                level = 1
                call alert(name,namel,'No prtc sigma',13,level)
         end if
1       continue
        call rline(name,namel,uprm)
        if (.not.find('DONE')) then
         level = 1
         call alert(name,namel,'Missing DONE after PRTC',23,level)
        end if
        totunq = maxunq
        write(stdo,101)
2       continue

c
c DONE with reading atom properties
c

c
c Check if end of file
c
        call rline(name,namel,uprm)
        if (find('*EOD')) then
         level = 0
         call alert(name,namel,'No 1-4 interaction parameters',29,level)
         return
        end if

c Reading potential energy parameters
c Start with 1-4 interactions:
c

        if (.not.find('1-4P')) then
                level = 1
                call alert(name,namel,'1-4 header missing',18,level)
        end if
        
        do coun=1,100
          call rline(name,namel,uprm)
          if (find('DONE')) then
            write(stdo,1022)
1022         format(//,1x,' *** DONE WITH 1-4 PROPERTIES *** ',//)
            go to 66
          end if
          
          el14  = getd('el14',el14)
          v14  = getd('v14f',v14)

        end do
66      continue

c
c Check if end of file
c
        call rline(name,namel,uprm)
        if (find('*EOD')) then
         level = 0
         call alert(name,namel,'No covalent structure',21,level)
         return
        end if

c Continue with bonds:
        if (.not.find('BOND')) then
                level = 1
                call alert(name,namel,'Bond header missing',19,level)
        end if

        do 5 coun=1,maxubnd
         call rline(name,namel,uprm)
         if (find('DONE')) then
          totubnd = coun - 1
          write(stdo,102)
102       format(//,1x,' *** DONE WITH BOND PROPERTIES *** ',//)
          go to 6
         end if
         success = .false.
         do 4 i=1,totunq
          if (debug) write(stdo,*)'basenm(i) = ',basenm(i)
          if (find(basenm(i))) then
           do 3 j=i,totunq
            if (debug) write(stdo,*)'basenm(j) = ',basenm(j)
            if (find(basenm(j))) then
                bondp(coun,1) = i
                bondp(coun,2) = j
                basekb(coun)  = number()
                basereq(coun) = number()
                success = .true.
                go to 5
            end if
3          continue
          end if
4        continue
         if (.not.success) then
          level = 1
          call alert(name,namel,'Unable to revive bond',21,level)
         end if
5       continue
         if (find('DONE')) then
          totubnd = maxubnd
          write(stdo,102)
         else
          level = 1
          call alert(name,namel,'NO DONE AFTER BOND LIST',23,level)
         end if

6       continue
c
c Check if end of file
c
        call rline(name,namel,uprm)
        if (find('*EOD')) then
         level = 0
         call alert(name,namel,'No angles and below ',20,level)
         return
        end if

c Read angle parameters
c

        if (.not.find('ANGL')) then
                level = 1
                call alert(name,namel,'Angle header missing',20,level)
        end if

        do 10 coun=1,maxuagl
         call rline(name,namel,uprm)
         if (find('DONE')) then
          totuagl = coun - 1
          write(stdo,103)
103       format(//,1x,' *** DONE WITH ANGLE PROPERTIES *** ',//)
          go to 11
         end if
         call get4c(char1,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char2,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char3,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         success = .false.
         anglp(coun,1) = 0
         anglp(coun,2) = 0
         anglp(coun,3) = 0
         do 9 i=1,totunq
          if (char1.eq.basenm(i)) anglp(coun,1) = i
          if (char2.eq.basenm(i)) anglp(coun,2) = i
          if (char3.eq.basenm(i)) anglp(coun,3) = i
9        continue
         if (anglp(coun,1).gt.0 .and. anglp(coun,2).gt.0 .and.
     1    anglp(coun,3).gt.0) then
                basekt(coun)  = number()
                baseaeq(coun) = number()
                baseaeq(coun)  = baseaeq(coun)/pi180
                totuagl = totuagl + 1
         else
          level = 1
          call alert(name,namel,'Unable to revive angle',22,level)
         end if
                
10      continue

        call rline(name,namel,uprm)
        if (find('DONE')) then
                totubnd = maxuagl
                write(stdo,103)
        else
                level = 1
                write(*,*)"ERROR: ",coun,maxuagl
        call alert(name,namel,'NO DONE AFTER ANGL LIST',23,level)
        end if
11      continue
c
c check if end of file
c
        call rline(name,namel,uprm)
        if (find('*EOD')) then
         level = 0
      call alert(name,namel,'No torsion parameters and below',31,level)
         return
        end if

c read torsion angle parameters
c
        if  (.not.find('TORS')) then
                level = 1
        call alert(name,namel,'Torsion header missing',22,level)
        end if

        do 16 coun=1,maxutrs
         call rline(name,namel,uprm)
         if (find('DONE')) then
          totutrs = coun - 1
          write(stdo,104)
104       format(//,1x,' *** DONE WITH TORSION PROPERTIES *** ',//)
          go to 17
         end if
         call get4c(char1,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char2,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char3,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char4,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         torsp(coun,1) = -1000
         torsp(coun,2) = -1000
         torsp(coun,3) = -1000
         torsp(coun,4) = -1000
         success = .false.
         do 15 i=1,totunq
          if (char1.eq.basenm(i)) torsp(coun,1) = i
          if (char2.eq.basenm(i)) torsp(coun,2) = i
          if (char3.eq.basenm(i)) torsp(coun,3) = i
          if (char4.eq.basenm(i)) torsp(coun,4) = i

c "X" is a generic particle. Global selection is allowed
c  only for particle 1 and particle 4 of a torsion

          if (char1.eq.'X') torsp(coun,1) = -999
          if (char4.eq.'X') torsp(coun,4) = -999

15       continue
         if (torsp(coun,1).gt.-1000 .and. torsp(coun,2).gt.0
     1    .and. torsp(coun,3).gt.0 .and. torsp(coun,4).gt.-1000) then
c       ileana
c       --testing the Amber force field implementation
c           write(stdo,*)'count in propert is ',coun


                  basamp1(coun) = number()
                  basamp2(coun) = number()
                  basamp3(coun) = number()
                  baseper(coun) = intgr()
                  basephs1(coun)= number()
                  basephs2(coun)= number()
                  basephs3(coun)= number()
         else
          level = 1
          call alert(name,namel,'Unable to revive torsion',24,level)
         end if
16      continue

c       write(*,*)' before a try to find done +++++++++++++'
        call rline(name,namel,uprm)
        if (find('DONE')) then
c        write(*,*)' found done +++++++++++++'
         totutrs = maxutrs
         write(stdo,104)
        else
         level = 1
         call alert(name,namel,'No done after angl list',23,level)
        end if
17      continue

c check if end of file
c
        call rline(name,namel,uprm)
        if (find('*EOD')) then
         level = 0
         call alert(name,namel,'No improper torsions',20,level)
         return
        end if


c read improper torsions
c
        if (.not.find('IMPR')) then
         level = 1
         call alert(name,namel,'Improper header missing',23,level)
        end if

        do 22 coun=1,maxuimp
         call rline(name,namel,uprm)
         if (find('DONE')) then
          totuimp = coun - 1
          write(stdo,105)
105       format(//,1x,' *** DONE WITH IMPROPER TORSIONS *** ',//)
          go to 23
         end if

         call get4c(char1,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char2,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char3,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         call get4c(char4,empty)
         if (empty) call alert(name,namel,'Empty atom name',15,1)
         improp(coun,1) = 0
         improp(coun,2) = 0
         improp(coun,3) = 0
         improp(coun,4) = 0
         if (arith) then
             if (char1.eq.'X') then
                improp(coun,1) = -999
             endif
             if (char2.eq.'X') then 
                improp(coun,2) = -999
             endif
             if (char3.eq.'X') then
                improp(coun,3) = -999
             endif
             if (char4.eq.'X') then
                improp(coun,4) = -999
             endif
          endif
         do 21 i=1,totunq
          if (char1.eq.basenm(i)) improp(coun,1) =i
          if (char2.eq.basenm(i)) improp(coun,2) =i
          if (char3.eq.basenm(i)) improp(coun,3) =i
          if (char4.eq.basenm(i)) improp(coun,4) =i
21       continue
c       ileana -- in improper list one can also have "X" atoms
c       also, baseieq is now cos(gamma), so there is no need to be converted

         if(.not.arith) then
         if (improp(coun,1).gt.0 .and. improp(coun,2).gt.0 .and.
     1          improp(coun,3).gt.0 .and. improp(coun,4).gt.0) then
                 baseimp(coun) = number()
                 baseieq(coun) = number()
                 baseieq(coun) = baseieq(coun)/pi180
         else
          level = 1
          call alert(name,namel,'Unable to revive imp. trs.',29,level)
         end if

         else

           
      if (improp(coun,1).gt.-1000.and. improp(coun,2).gt.-1000.and.
     1  improp(coun,3).gt.-1000.and. improp(coun,4).gt.-1000) then
                 baseimp(coun) = number()
c       here baseieq has no meaning with the Amber prop file
c                baseieq(coun) = number()
         else
          level = 1
          call alert(name,namel,'Unable to revive imp. trs.',29,level)
         end if
         endif   
22      continue

        call rline(name,namel,uprm)
        if (find('DONE')) then
         totuimp = maxuimp
         write(stdo,105)
        else
         level = 1
         call alert(name,namel,'Missing DONE after impropers',28,level)
        end if
23      continue
        write(stdo,106)totunq,totubnd,totuagl,totutrs,totuimp
106     format(//,1x,' *** Number of parameters read: ',/,1x,
     1   'Unique particles, bonds, angles, torsions, improper torsions',
     2   /,8x,5(3x,i7),//)
        return
        end
