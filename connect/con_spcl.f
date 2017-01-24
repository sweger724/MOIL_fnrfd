        program con_xite
c CREATE SPECIAL CONNECTIVITY FILES FOR FORMATION AND DESTRUCTION OF BONDS
c
c----------------------------------------------------------------
c LENGTH.BLOCK  - includes defintion for vector lengths
c CONNECT.BLOCK - includes the topology and energy term for the
c                       molecule
c PROPERT.BLOCK - includes definition for atom bond angle torsion
c                       and improper torsion properties
c MONOMERS.BLOCK- includes information on connectivity of individual
c                       monomers
c LINE.BLOCK    - used for line interpreter
c UNITS.BLOCK   - standard units used
c DEBUG.BLOCK   - To initiate debugging printout
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/PROPERT.BLOCK'
        include 'COMMON/MONOMERS.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SPECL.BLOCK'

        integer rvnewpoit(maxpt)
c----------------------------------------------------------------


c----------------------------------------------------------------
c of    - integer function return a unit number of file to be opened
c UNITS that we need
c uwcon - unit to write conn. and properties of the molecule (output)
        integer of 
        integer uwcon,urcon
c----------------------------------------------------------------
c ******************************************************
        integer ii,jj,kk,ipick(maxpt),npick,l,ll,ls
        integer mm,tmpexc1,npt1,npt2,npt3
c **************************************************** 

c----------------------------------------------------------------
c local variables
c
        character*4 name
        logical find
        integer unit ,imon,j,i,k,m
        integer namel,tmp(100)
        integer level
        
        double precision getd


        data uwcon,urcon/2*-1/
        data el14,v14,eps/2.0d0,8.0d0,1.d0/
c--------------------------------------------------------------

c-------------------------------------------------------------------
c initialize file units for standard i/o
c Initialize some molecular properties (see CONNECT.BLOCK for details)
c and unit numbers (to -1)
c
        stdi = 5
        stdo = 6
        stderr = 0
        totmon = 0
        totex  = 0
        npt    = 0
        nb     = 0
        nmb    = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
c--------------------------------------------------------------------
        specl=.false.
        emyes0 = .false.
        debug = .false.
        name  = 'CONS'
        namel = 4
        j = 0
        k = 0
        m = 0
        imon = 0
        ls = 0
        l  = 0
c----------------------------------------------------------------
c open junk file unit (25) and
c Input to connect is set to standard input (stdi)
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
        unit = stdi
c----------------------------------------------------------------
c intention: to open files and read variables
c
1       continue
        call rline(name,namel,unit)
        v14  = getd('v14f',v14)
        el14 = getd('el14',el14)
        eps  = getd('epsi',eps)
        if (find('debu')) debug = .true.
        if (find('file')) then
          
         if (find('rcon')) then
          urcon  = of()
          call rconn(urcon)
         end if
         if (find('wcon')) uwcon = of()
         else if (find('spec')) then
          specl=.true.

ccccccccccccccccccccccccccccccccccccccccccccc
         else if (find('mos1')) then
            call pick(ipick,npick)
             do 44 i=1,npt
              if(ipick(i) .gt. 0 ) then 
               j = j + 1
               if (j.gt.maxspcl) then
                 level = 1
                 call alert(name,namel,'maxspcl exceeded!',17,level)
               end if
               ntyp(j) = npick  
               newpoit(j) = i
               rvnewpoit(i) = j
               if (poimon(i) .ne. poimon(i-1)) then
                imon=imon+1
                moname(imon)=moname(poimon(i))
                tmp(imon) = poipt(poimon(i))
               end if
               poimon(j) = imon
             end if
44         continue
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
           npt1 = j
           do 59 mm=1,j
            i = newpoit(mm)
               tmpexc1=exc1(i)-exc1(i-1)
               exc1(mm) = tmpexc1 + l
               do 53 ll = exc1(i-1)+1,exc1(i)
                if(rvnewpoit(exc2(ll)) .eq. 0 ) then
                 exc1(mm) = exc1(mm)-1
                 goto 53
                end if
                exc2(l+1) = rvnewpoit(exc2(ll))
                l=l+1
53             continue
59         continue

           do 54 ll = 1,totspe
            ii = spec1(ll)
            jj = spec2(ll)
            if ((ipick(ii).gt.0) .and. (ipick(jj).gt.0)) then
                ls = ls + 1
                spec1(ls) = spec1(ll)
                spec2(ls) = spec2(ll)
                p14(1,ls) = p14(1,ll)
                p14(2,ls) = p14(2,ll)
                p14(3,ls) = p14(3,ll)
            end if
54         continue
           lz14(1) = ls
          
c++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
              do 11 i=1,imon
                poipt(i)=rvnewpoit(tmp(i))
11            continue
             poitype(1) = j
            do 2000 i=1,nb
            ii=ib1(i)
            jj=ib2(i)
            if(ipick(ii) .gt.0 .and. ipick(jj) .gt.0) then
             k = k + 1
             ib1(k) = rvnewpoit(ib1(i))
             ib2(k) = rvnewpoit(ib2(i))
             req(k) = req(i) 
             kbond(k) = kbond(i)
            end if
2000       continue
           snb(1) = k
            do 3000 i=1,nangl     
             ii=iangl1(i)
             kk=iangl3(i)
             jj=iangl2(i)
           if (ipick(ii).gt.0 .and. ipick(jj).gt.0) then
            if (ipick(kk).gt.0)then 
             m=m+1
             iangl1(m) = rvnewpoit(iangl1(i))
             iangl2(m) = rvnewpoit(iangl2(i))
             iangl3(m) = rvnewpoit(iangl3(i))
             kangl(m) = kangl(i)
             angleq(m) = angleq(i)
            end if
           end if
3000      continue
          snang(1) = m
c  ccccccccccccccccccccccccccccccccccccccccccccc
         else if(find('mos2')) then
          call pick(ipick,npick)
             do 144 i=1,npt
              if(ipick(i) .gt. 0 ) then 
               j = j + 1
               ntyp(j) = npick + 1 
               newpoit(j) = i
               rvnewpoit(i) = j
               if (poimon(i) .ne. poimon(i-1)) then
                imon=imon+1
                moname(imon)=moname(poimon(i))
                tmp(imon) = poipt(poimon(i))
               end if
               poimon(j) = imon
             end if
144         continue
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
           npt2 = j
           do 69 mm=npt1+1,j
            i = newpoit(mm)
               tmpexc1=exc1(i)-exc1(i-1)
               exc1(mm) = tmpexc1 + l
               do 63 ll = exc1(i-1)+1,exc1(i)
                if(rvnewpoit(exc2(ll)) .eq. 0 ) then
                 exc1(mm) = exc1(mm)-1
                 goto 63
                end if
                exc2(l+1)=rvnewpoit(exc2(ll))
                l=l+1
63             continue
69          continue
               do 154 ll = 1,totspe
                        ii = spec1(ll)
                        jj = spec2(ll)
                        if (ipick(ii)*ipick(jj).gt.0) then
                         ls = ls + 1
                         spec1(ls) = spec1(ll)
                         spec2(ls) = spec2(ll)
                         p14(1,ls) = p14(1,ll)
                         p14(2,ls) = p14(2,ll)
                         p14(3,ls) = p14(3,ll)
                        end if
154            continue
               lz14(2) = ls
          
c++++++++++++++++++++++++++++++++++++++++++++++++++++++
              do 51 i=1,imon
                poipt(i)=rvnewpoit(tmp(i))
51            continue
             poitype(2) = j
            do 2010 i=1,nb
            ii=ib1(i)
            jj=ib2(i)
            if(ipick(ii) .gt.0 .and. ipick(jj) .gt.0) then
             k = k + 1
             ib1(k) = rvnewpoit(ib1(i))
             ib2(k) = rvnewpoit(ib2(i))
             req(k) = req(i) 
             kbond(k) = kbond(i)
            end if
2010       continue
           snb(2) = k
            do 3001 i=1,nangl     
             ii=iangl1(i)
             kk=iangl3(i)
             jj=iangl2(i)
           if (ipick(ii).gt.0 .and. ipick(jj).gt.0) then
            if (ipick(kk).gt.0)then 
              m=m+1
              iangl1(m) = rvnewpoit(iangl1(i))
              iangl2(m) = rvnewpoit(iangl2(i))
              iangl3(m) = rvnewpoit(iangl3(i))
              kangl(m) = kangl(i)
              angleq(m) = angleq(i)
             end if
            end if
3001       continue
          snang(2) = m
c  ccccccccccccccccccccccccccccccccccccccccccccc
        else if(find('mos3')) then
         call pick(ipick,npick)
             do 244 i=1,npt
              if(ipick(i) .gt. 0 ) then 
               j = j + 1
               ntyp(j) = npick + 2 
               newpoit(j) = i
               rvnewpoit(i) = j
               if (poimon(i) .ne. poimon(i-1)) then
                imon=imon+1
                moname(imon)=moname(poimon(i))
                tmp(imon) = poipt(poimon(i))
               end if
               poimon(j) = imon
             end if
             poitype(3) = j
244         continue
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
           npt3 = j
           do 79 mm=npt2+1,j
            i = newpoit(mm)
               tmpexc1=exc1(i)-exc1(i-1)
               exc1(mm) = tmpexc1 + l
               do 73 ll = exc1(i-1)+1,exc1(i)
                if(rvnewpoit(exc2(ll)) .eq. 0 ) then
                 exc1(mm) = exc1(mm)-1
                 goto 73
                end if
                exc2(l+1)=rvnewpoit(exc2(ll))
                l=l+1
73             continue
79          continue
               do 254 ll = 1,totspe
                ii = spec1(ll)
                jj = spec2(ll)
                if (ipick(ii)*ipick(jj).gt.0) then
                        ls = ls + 1
                         spec1(ls) = spec1(ll)
                         spec2(ls) = spec2(ll)
                         p14(1,ls) = p14(1,ll)
                         p14(2,ls) = p14(2,ll)
                         p14(3,ls) = p14(3,ll)
                end if
254             continue
                lz14(3) = ls
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
              do 71 i=1,imon
                poipt(i)=rvnewpoit(tmp(i))
71            continue
            do 2210 i=1,nb
            ii=ib1(i)
            jj=ib2(i)
            if(ipick(ii) .gt.0 .and. ipick(jj) .gt.0) then
             k = k + 1
             ib1(k) = rvnewpoit(ib1(i))
             ib2(k) = rvnewpoit(ib2(i))
             req(k) = req(i) 
             kbond(k) = kbond(i)
            end if
2210       continue
           snb(3) = k
            do 3201 i=1,nangl     
             ii=iangl1(i)
             kk=iangl3(i)
             jj=iangl2(i)
           if (ipick(ii).gt.0 .and. ipick(jj).gt.0) then
            if (ipick(kk).gt.0)then 
             m=m+1
             iangl1(m) = rvnewpoit(iangl1(i))
             iangl2(m) = rvnewpoit(iangl2(i))
             iangl3(m) = rvnewpoit(iangl3(i))
             kangl(m) = kangl(i)
             angleq(m) = angleq(i)
            end if
           end if
3201       continue
          snang(3) = m
c  ccccccccccccccccccccccccccccccccccccccccccccc
        else if(find('mos4')) then
         call pick(ipick,npick)
             do 294 i=1,npt
              if(ipick(i) .gt. 0 ) then 
               j = j + 1
               ntyp(j) = npick + 3 
               newpoit(j) = i
               rvnewpoit(i) = j
               if (poimon(i) .ne. poimon(i-1)) then
                imon=imon+1
                moname(imon)=moname(poimon(i))
                tmp(imon) = poipt(poimon(i))
               end if
               poimon(j) = imon
            end if
294         continue
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
           do 89 mm=npt3+1,j
            i = newpoit(mm)
               tmpexc1=exc1(i)-exc1(i-1)
               exc1(mm) = tmpexc1 + l
               do 83 ll = exc1(i-1)+1,exc1(i)
                if(rvnewpoit(exc2(ll)) .eq. 0 ) then
                 exc1(mm) = exc1(mm)-1
                 goto 83
                end if
                exc2(l+1)=rvnewpoit(exc2(ll))
                l=l+1
83             continue
89          continue
               do 354 ll = 1,totspe
                ii = spec1(ll)
                jj = spec2(ll)
                if(ipick(ii)*ipick(jj).gt.0) then
                        ls = ls + 1
                         spec1(ls) = spec1(ll)
                         spec2(ls) = spec2(ll)
                         p14(1,ls) = p14(1,ll)
                         p14(2,ls) = p14(2,ll)
                         p14(3,ls) = p14(3,ll)
                end if
354             continue
                lz14(4) = ls
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
              do 91 i=1,imon
                poipt(i)=rvnewpoit(tmp(i))
91            continue
             poitype(4) = j

            do 2219 i=1,nb
            ii=ib1(i)
            jj=ib2(i)
            if(ipick(ii) .gt.0 .and. ipick(jj) .gt.0) then
             k = k + 1
             ib1(k) = rvnewpoit(ib1(i))
             ib2(k) = rvnewpoit(ib2(i))
             req(k) = req(i) 
             kbond(k) = kbond(i)
            end if
2219       continue
           snb(4) = k
            do 9201 i=1,nangl     
             ii=iangl1(i)
             kk=iangl3(i)
             jj=iangl2(i)
           if (ipick(ii).gt.0 .and. ipick(jj).gt.0) then
            if (ipick(kk).gt.0)then 
             m=m+1
             iangl1(m) = rvnewpoit(iangl1(i))
             iangl2(m) = rvnewpoit(iangl2(i))
             iangl3(m) = rvnewpoit(iangl3(i))
             kangl(m) = kangl(i)
             angleq(m) = angleq(i)
            end if
           end if
9201       continue
         snang(4) = m
ccccccccccccccccccccccccccccccccccccccccccccccccc


         else if (find('acti')) then
          go to 2
        end if
        go to 1
2       continue
           totmon=imon
           npt = j
           print*,'npt=',npt
           print*,'m=',m
           nangl = m
           nb    = k
            totex = l
            totspe = ls
            print*,'totex=',totex
            print*,'totspe=',totspe
                write(*,*)' NBULK = ',NBULK
         pbulk(NBULK) = npt
        

        if (debug) then
         write(stdo,*)' After opening files '
        end if
c
c
c ******************************************************************
c *****************************************************************
c Call wconn: write down current topolopgy and properies of molecular
c system
        call wconn(uwcon,v14,el14,0)
        do 3 i=1,NBULK
         write(stdo,100)BULK(i)
100      format(/,1x,'*** MOLECULE ',a4, ' ASSEMBLED ***')
3       continue
        stop
        end
