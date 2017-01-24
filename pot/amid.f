        subroutine amid()

        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'

c kamid  - amide plane torsion constraint value
c ncnsttmp - temporary ncnst, used to tell how many amide constraints added

        double precision kamid,getd
        integer ncnsttmp,level,namel,i,j
        character*4 name


        name = 'amid'
        namel = 4

        kamid = getd('kamd',100.d0)
        ecnyes = .true.
c
c pull them from torsion list
c checking for O-C-N-H or H-N-C-O
c
        ncnsttmp=ncnst
        do 112 i=1,ntors
                if (ptnm(itor1(i))(1:1).eq.'O') then
                  if (ptnm(itor2(i))(1:1).eq.'C') then
                    if (ptnm(itor3(i))(1:1).eq.'N') then
                      if ((ptnm(itor4(i))(1:1).eq.'H'))then 
c michele oplsaal
c                    .or.
c     1                    ((ptnm(itor4(i))(2:2).eq.'H') .and.
c     1                     (ptnm(itor4(i))(1:1) .eq. '1'))) then
                 
c  correction for OPLSAAL ... Peter                 
                if (ptnm(itor4(i)).eq.'HE22') goto 112
c michele
                if (ptnm(itor4(i)).eq.'H2') goto 112
                if (ptnm(itor4(i)).eq.'HD22') goto 112
c end of michele


c
c check for things like ASN where potentially two amides could
c be found for same plane (two H on N)
c look through constraints already assigned for match of first 3 atoms
c (but not last one)
                        if (ncnst.ne.0 .and. ncnsttmp.ne.0) then
                        do 10 j=ncnsttmp,ncnst
                          if (itor1(i).eq.icnst1(j).and.
     1                     itor2(i).eq.icnst2(j).and.
     2                     itor3(i).eq.icnst3(j)) go to 112
10                     continue
                        end if

                        ncnst = ncnst + 1
                        if (ncnst.gt.maxcnst) then
                          level = 1
                          call alert(name,namel,'Maxcnst exceeded',
     1                         16,level)
                        endif
                        icnst1(ncnst) = itor1(i)
                        icnst2(ncnst) = itor2(i)
                        icnst3(ncnst) = itor3(i)
                        icnst4(ncnst) = itor4(i)
                        kcns(ncnst)   = kamid
                        cnseq(ncnst)  = 180.d0/pi180
                      !write(6,*)"YY",itor1(i),itor2(i),itor3(i),itor4(i)


                        if (debug) then
                                write(stdo,*)' kcns cnseq ',
     &                          kcns(ncnst),cnseq(ncnst)*pi180
                        end if
                      endif
                    endif
                  endif
                  else if ((ptnm(itor1(i))(1:1).eq.'H')) then

c michele oplsaal
c                  .or.
c     1                     ((ptnm(itor1(i))(2:2).eq.'H') .and.
c     1                      (ptnm(itor1(i))(1:1) .eq. '1'))) then



                if (ptnm(itor1(i)).eq.'HE22') goto 112
                if (ptnm(itor1(i)).eq.'H2') goto 112
                if (ptnm(itor1(i)).eq.'HD22') goto 112
c end of michele



                     if (ptnm(itor2(i))(1:1).eq.'N') then
                       if (ptnm(itor3(i))(1:1).eq.'C') then
                         if (ptnm(itor4(i))(1:1).eq.'O') then

                     write(6,*)"XX",itor1(i),itor2(i),itor3(i),itor4(i)
c check for things like ASN where potentially two amides could
c be found for same plane (two H on the N)
c look through constraints already assigned for match of last 3 atoms
c (but not first one)
                        if (ncnst.ne.0 .and. ncnsttmp.ne.0) then
                        do 11 j=ncnsttmp,ncnst
                          if (itor4(i).eq.icnst4(j).and.
     1                     itor2(i).eq.icnst2(j).and.
     2                     itor3(i).eq.icnst3(j)) go to 112
11                     continue
                        end if
c

                        ncnst = ncnst + 1
                        if (ncnst.gt.maxcnst) then
                          level = 1
                          call alert(name,namel,'Maxcnst exceeded',
     1                         16,level)
                        endif
                        icnst1(ncnst) = itor1(i)
                        icnst2(ncnst) = itor2(i)
                        icnst3(ncnst) = itor3(i)
                        icnst4(ncnst) = itor4(i)
                        kcns(ncnst)   = kamid
                        cnseq(ncnst)  = 180.d0/pi180

c michele 
          write(132,*)"secondo",itor1(i),itor2(i),itor3(i),itor4(i)



                        if (debug) then
                                write(stdo,*)' kcns cnseq ',
     &                          kcns(ncnst),cnseq(ncnst)*pi180
                        end if
                      endif
                    endif
                  endif
                endif
112         continue
           write (stdo,100) ncnst-ncnsttmp,kamid
100        format(1x,'Constrained ',i7,
     1      ' amide planes with f cnst ',f8.3)

        return
        end
