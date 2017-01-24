        subroutine rconn(urcon)
        implicit none
c
c Read connectivity data, the one stored in
c connectivity common block: CONNECT.BLCOK
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/MUTA.BLOCK'

        integer urcon,namel
        integer i,ii,j,k,k1,l,ll,level
        double precision tmp

        character*4 version,vrs_cmp
        character*5 name
        character*1 TILDA

        integer istart,iend,total_n
        logical trash

c hydrogen bonds
        integer donor_index,acceptor_index
        character*1 hb
c end hydrogen bonds

        parameter (version='12.1')
        data name/'rconn'/
        data namel/5/

        
        rewind urcon

        call init_var()

c hydrogen bonds initialization
        donor_index = 1
        acceptor_index = 1

        read(urcon,99)vrs_cmp
99      format(a4)
        if (vrs_cmp.ne.version) then
                level = 1
                write(*,*)' *** Version on con file ',vrs_cmp
                write(*,*)' *** but expecting... ', version
                write(*,*)' *** RECREATE YOUR CONNECTIVITY FILE *** '
                call alert(name,namel,' con file too old',17,level)
        end if

        read(urcon,100)TILDA
100     format(a1)
        read(urcon,100)TILDA
c totmon npt nb nangl ntors nimp totex totspe lestyp NBULK
        read(urcon,101)totmon,npt,nb,nangl,ntors,nimp,totex,totspe,
     1lestyp,NBULK,nmb,nchg,nbeta,nwaters,arith,trash,mdivyes
101     format(1x,14i8,3l2)
c@@@@@@@@
c       write(*,101)totmon,npt,nb,nangl,ntors,nimp,totex,totspe,
c     1  lestyp,NBULK,nmb,nchg,nbeta,nwaters,arith,prll_on_off,mdivyes
c@@@@@@@@@

        if (NBULK.EQ.0) then
         level = 0
         call alert(name,namel,'No molecules assembled',22,level)
         return
        else if (totmon.eq.0) then
         level = 0
         call alert(name,namel,'No monomers ?!',14,level)
         return
        else if (npt.eq.0) then
         level = 0
         call alert(name,namel,'No particles ?!',15,level)
         return
        end if


        poipt(0) = 0
        do 1 i=1,NBULK
         read(urcon,102)BULK(I)
102      format(a)
1       continue
        read(urcon,103)(pbulk(i),i=1,nbulk)
103     format(/,10(1x,i7))
        read(urcon,104)(moname(i),i=1,totmon)
104     format(/,10(1x,a4))
        read(urcon,105)(poipt(i),i=1,totmon)
105     format(/,10(1x,i7))
        read(urcon,100)TILDA
        read(urcon,100)TILDA
        do 2 i=1,npt
         read(urcon,107)j,poimon(i),ptid(i),lesid(i),mutaid(i),ptnm(i)
     1       ,ptms(i),ptchg(i),epsgm6(i),epsgm12(i),ptwei(i),ptsaid(i)
     2       ,hb
107      format(1x,i7,1x,i7,1x,i3,1x,i3,1x,i3,1x,a4,1x,f7.3,2x,f9.5,
     1          1x,e12.5,e12.5,e10.3,1x,i3,1x,a1)
        if (debug) then
         write(stdo,*) ' at 107'
         write(stdo,107)j,poimon(i),ptid(i),lesid(i),mutaid(i),ptnm(i)
     1       ,ptms(i),ptchg(i),epsgm6(i),epsgm12(i),ptwei(i),ptsaid(i)
     2       ,hb
        end if
c hydrogen bonds
        if (hb.eq.'A')then
            poi_acc(acceptor_index)=i
            acceptor_index=acceptor_index + 1
        elseif (hb.eq.'D')then
            poi_don(donor_index)=i
            donor_index=donor_index + 1
        endif
c end hydrogen bonds
2       continue
c hydrogen bonds
        n_donors=donor_index - 1
        n_acceptors= acceptor_index - 1
c end hydrogen bonds
        tmp = 0.d0
        do i=1,npt
          tmp = tmp + ptchg(i)
        end do
c@        write(*,*)' The total charge of the system is ',tmp

        if (nb.gt.0) then

         total_n = nb
         if (prll_on_off) then
          call load_balance(nb,my_pe,num_pes,istart,iend)
          bonds_start=istart
          bonds_end=iend
         else
          istart = 1
          iend   = nb
          bonds_start=istart
          bonds_end=iend
         end if



          read(urcon,100)TILDA
          read(urcon,100)TILDA

           do i=1,total_n
            j = i
            read(urcon,109)ib1(j),ib2(j),kbond(j),req(j)
109         format(2(1x,i7),2(1x,f10.4))
          end do

          nb_all = total_n

        end if

c@
C       write(*,*)' Last bond is ',ib1(nb),ib2(nb),kbond(nb),req(nb)
C       write(*,*)' PI180 = ',pi180
C       write(*,*)' sin (PI) = ',sin(pi180*180)
c------------------------------------------------------
        if (nmb.gt.0) then

         emyes0 = .true.
         total_n = nmb
         if (prll_on_off) then
          call load_balance(nmb,my_pe,num_pes,istart,iend)
         else
          istart = 1
          iend   = nmb
         end if

         read(urcon,100)TILDA
         read(urcon,100)TILDA

         do 198 i=1,istart-1
          read(urcon,100)TILDA
198      continue

         do 199 i=istart,iend
          j = i - istart + 1
          read(urcon,189)imb1(j),imb2(j),rmeq(j),D(j),alpha(j)
189       format(2(1x,i7),3(1x,f10.4))
          emyes(i)  = .true.
199      continue

          do 200 i=iend+1,total_n
           read(urcon,100)TILDA
200       continue

        end if
c----------------------------------------------------

        if (nangl.gt.0) then

         total_n = nangl
         if (prll_on_off) then
          call load_balance(nangl,my_pe,num_pes,istart,iend)
         else
          istart = 1
          iend   = nangl
         end if
c@
C       write(*,*)'Angle: my_pe num_pes istart iend ',my_pe,num_pes,istart,iend
          read(urcon,100)TILDA
          read(urcon,100)TILDA

          do i=1,istart-1
             !read(urcon,100)TILDA
             j = i+iend-istart+1
             read(urcon,111)iangl1(j),iangl2(j),iangl3(j),
     1            kangl(j),angleq(j)
             angleq(j) = angleq(j)/pi180
          end do

          do 42 i=istart,iend
           j = i - istart + 1
           read(urcon,111)iangl1(j),iangl2(j),iangl3(j),
     1          kangl(j),angleq(j)
C@
C          write(*,*)i,iangl1(j),iangl2(j),iangl3(j),kangl(j),angleq(j)
111        format(3(1x,i7),2(1x,f10.5))
           angleq(j) = angleq(j)/pi180
42        continue

          do 43 i=iend+1,total_n
            !read(urcon,100)TILDA
            read(urcon,111)iangl1(i),iangl2(i),iangl3(i),
     1            kangl(i),angleq(i)
            angleq(i) = angleq(i)/pi180
43        continue

          nangl_all = total_n

        end if

        if (ntors.gt.0) then

         total_n = ntors
         if (prll_on_off) then
          call load_balance(ntors,my_pe,num_pes,istart,iend)
         else
          istart = 1
          iend   = ntors
         end if
c@
C         write(*,*)' my-pe ntors istart iend ',my_pe,ntors,istart,iend

          read(urcon,100)TILDA
          read(urcon,100)TILDA

          do 51 i=1,istart-1
           read(urcon,100)TILDA
51        continue

          do 52 i=istart,iend
           j = i - istart +1
           read(urcon,113)itor1(j),itor2(j),itor3(j),itor4(j),
     1          period(j),ktors1(j),ktors2(j),ktors3(j),
     2          phase1(j),phase2(j),phase3(j)
113        format(5(1x,i7),3(1x,f10.4),1x,3(f6.2))
52        continue

          do 53 i=iend+1,total_n
           read(urcon,100)TILDA
53        continue

        end if

        if (nimp.gt.0) then

        total_n = nimp
         if (prll_on_off) then
          call load_balance(nimp,my_pe,num_pes,istart,iend)
         else
          istart = 1
          iend   = nimp
         end if
C@
C       write(*,*)' nimp istart iend ',nimp,istart,iend
C       write(*,*)' total_n = ',total_n

          read(urcon,100)TILDA
          read(urcon,100)TILDA

          do 61 i=1,istart-1
           read(urcon,100)TILDA
61        continue

          do 6 i=istart,iend
           j = i - istart + 1
           read(urcon,115)iimp1(j),iimp2(j),iimp3(j),iimp4(j),
     1kimp(j),impeq(j)
115        format(4(1x,i7),2(1x,e15.8))
           impeq(j) = impeq(j)/pi180
6         continue

          do 62 i=iend+1,total_n
           read(urcon,100)TILDA
62        continue
        end if
        
        if (totex.gt.0) then
c the exclusion list is not split between the processors
c since the list for searching of neighbours varies in the simulation
c
         read(urcon,100)TILDA
c@
C       write(*,100)TILDA
         read(urcon,100)TILDA
c@
C       write(*,100)TILDA
         k = 0
         i = 0
         exc1(0) = 0
7        continue
         read(urcon,117)j,k1
71       continue
         i = i + 1
         if (debug) write(stdo,*)' j k1 = ',j,k1

C@
C       write(*,*) ' reading totex...j k1 = ',j,k1

117        format(2(1x,i7))
           if (i.lt.j) then
            exc1(i) = k
            go to 71
         else if (i.gt.j) then
            level = 1
            call alert(name,namel,'Fishy exclusion list!',21,level)
         else
            exc1(j) = k + k1
            read(urcon,118)(exc2(l),l=k+1,exc1(j))
118         format(1x,10(i7,1x))
            k = k + k1
         end if
         if (k.ne.totex) go to 7
        end if

C reading special exclusion list for muta
        if (muta) then
        read (urcon,100) TILDA
         k = 0
         j = 0
         l = 0
         ii = 1
         do i=1,npt
           if (mutaid(i).gt.0) then
           read(urcon,117)j,k1
           exm1(ii) = k+k1
           read(urcon,118)(exm2(l),l=k+1,exm1(ii))
           write(stdo,*)(exm2(l),l=k+1,exm1(ii))
           k = k+k1
           ii = ii + 1
           write (stdo,*) j,exm1(j)
           write (stdo,*) k+1,exm2(k+1),exm2(exm1(j))
          endif

         enddo
         endif
        
        if (totspe.gt.0) then

         total_n = totspe
         if (prll_on_off) then
          call load_balance(totspe,my_pe,num_pes,istart,iend)
         else
          istart = 1
          iend   = totspe
         end if

         read(urcon,100)TILDA
         read(urcon,100)TILDA

         do 81 i=1,istart-1
          read(urcon,100)TILDA
81       continue
         do 8 i=istart,iend
          j = i - istart + 1
         read(urcon,120)spec1(j),spec2(j),p14(1,j),p14(2,j),p14(3,j)
c120      format(2(1x,i5),3(1x,f15.5))
120       format(2(1x,i7),1x,f30.5,2(1x,f15.5))
8        continue

         do 82 i=iend+1,total_n
          read(urcon,100)TILDA
82       continue

        end if

        read(urcon,100)TILDA
        read(urcon,123)(flagchr(i),i=1,npt)
123     format(1x,15l2)
        read(urcon,100)TILDA
        read(urcon,124)(cplbl(i),i=1,npt)
124     format(1x,15i6)
        if (nbeta.gt.0) then
                
                total_n = nbeta
                if (prll_on_off) then
                 call load_balance(nbeta,my_pe,num_pes,istart,iend)
                else
                 istart = 1
                 iend   = nbeta
                end if

                read(urcon,100)TILDA
                read(urcon,100)TILDA

                read(urcon,125)(j,i=1,istart-1)
125             format(1x,10i8)
                read(urcon,125)(betap(i-istart+1),i=istart,iend)
                read(urcon,125)(j,i=iend+1,total_n)

                read(urcon,126)(tmp,i=1,istart-1)
126             format(1x,10f6.3)
                read(urcon,126)(cbeta(i-istart+1),i=istart,iend)
                read(urcon,126)(tmp,i=iend+1,total_n)
                read(urcon,127)avg_hydro
127             format(1x,2f9.4)
        end if

         if (totmon.gt.maxmono) then
         level = 1
         write(stdo,*)' totmon = ',totmon,'  maxmono = ',maxmono
         call alert(name,namel,' number of mono exceeded ',25,level)
        else if (npt.gt.maxpt) then
         level = 1
         write(stdo,*)' npt  = ',npt,'  maxpt = ',maxpt
         call alert(name,namel,' number of prtc exceeded ',25,level)
        else if (nb.gt.maxbond) then
         level = 1
         write(stdo,*)' nb = ',nb,'  maxbond = ',maxbond
         call alert(name,namel,' number of bond exceeded ',25,level)
        else if(nangl.gt.maxangl) then
         level = 1
         write(stdo,*)' nangl = ',nangl,'  maxangl = ',maxangl
         call alert(name,namel,' number of angl exceeded ',25,level)
        else if(ntors.gt.maxtors) then
         level = 1
         write(stdo,*)' ntors = ',ntors,'  maxtors = ',maxtors
         call alert(name,namel,' number of tors exceeded ',25,level)
        else if(nimp.gt.maximp) then
         level = 1
         write(stdo,*)' nimp = ',nimp,'  maximp =',maximp
         call alert(name,namel,' number of impr exceeded ',25,level)
        else if(totex.gt.maxex) then
         level = 1
         write(stdo,*)' totex = ',totex,'  maxex =',maxex
         call alert(name,namel,' number of excl exceeded ',25,level)
        else if(totspe.gt.maxspec) then
         level = 1
         write(stdo,*)' totspe = ',totspe,' maxspec = ',maxspec
         call alert(name,namel,' number of spec exceeded ',25,level)
        else if (nmb.gt.maxmorsb) then
         level = 1
         write(stdo,*)' nmb = ',nmb,' maxmorsb = ',maxmorsb
         call alert(name,namel,' Number of morseb exceeded ',27,level)
        else if (NBULK.gt.MAXBULK) then
         level = 1
         write(stdo,*)' NBULK = ',NBULK, '  MAXBULK = ',MAXBULK
         call alert(name,namel,' number of BULK exceeded ',25,level)
        end if
c yael
c read the monomers division data

      if (mdivyes) then
        read(urcon,100)TILDA
        do 240 i=1,totmon
             read(urcon,210)mdivlist(0,i)
           if (mdivlist(0,i).ne.-1) then
                read(urcon,210)(mdivlist(j,i),j=1,mdivlist(0,i)+1)
           else
               read(urcon,210)(mdivlist(j,i),j=1,2)
           endif
240     continue
210     format(1x,7i9)
      else
         do 250 i=1,totmon
            mdivlist(0,i) = 1
            mdivlist(1,i) = poipt(i-1)
            mdivlist(2,i) = poipt(i)
250      continue
      endif

      call filldpt(mdivlist,dpoipt,poidmon,poipt,realmono,totdmon,
     1          totmon,npt)
c
c if number of TIP3 waters is not zero read their monomer # below
c
         if (nwaters.gt.0) then
          read(urcon,100)TILDA
          read(urcon,211)(idxtip3(i),i=1,nwaters)
211       format(10i8)
         end if
        
      return
      end

c---------------------------------------------------

c yael.
c create the realmono and the dpoipt lists using the mdivlist
c data .dpoipt(i) points to the last particle belong to the i'th
c monomer part. realmono(i) points the monomer the i'th part belongs too.
c totdmon is the total number of monomer parts.
      subroutine filldpt(mdivlist,dpoipt,poidmon,poipt,realmono,
     1          totdmon,totmon,npt)
      include 'COMMON/LENGTH.BLOCK'
      integer mdivlist(0:maxdiv+1,maxmono)
      integer dpoipt(0:maxmono*(maxdiv+1)),realmono(maxmono*maxdiv)
      integer poipt(0:maxmono)
      integer poidmon(maxpt)
      integer totdmon,totmon,npt
c local
      integer i,j,pt
      totdmon = 0
      do 100 i=1,totmon
         if (mdivlist(0,i).eq.-1) then
            dpoipt(totdmon) = poipt(i-1)
            totdmon = totdmon + 1
            realmono(totdmon) = i
         endif
         do 200 j=1,mdivlist(0,i)
            dpoipt(totdmon) = mdivlist(j,i)
            totdmon = totdmon + 1
            realmono(totdmon) = i
200      continue
100   continue

        dpoipt(totdmon) = npt
        realmono(totdmon) = totmon

c       do 400 i=1,totdmon
c         do 300 pt=dpoipt(i-1)+1,dpoipt(i)
c            poidmon(pt) = i
c300      continue
c400   continue

        return
        end
