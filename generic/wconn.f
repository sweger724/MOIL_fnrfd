c hk change start
c originally, subroutine wconn(uwcon,v14,el14)
        subroutine wconn(uwcon,v14,el14,flag_pair)
c       hk change end

c
c Write down connectivity data, mostly the one stored in
c connectivity common block: CONNECT.BLCOK
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'


        integer uwcon,namel
        integer level,i,j,k,ll
        character*4 vrsn
        character*5 name
        character*1 TILDA
        double precision tmp,v14,el14

cccccccccccccccccccc huhnkie start ccccccccccccccccccccccc
c type_arr is an array that non-redundantly saves all the existing pair types.
c at maximum, there can be "#types*#types" pair types.
c (maxunq is the number of possible types.)
c (type_arr(1, i), type_arr(2, i)) is the pair type stored 
c at ith position in the array. 
        integer type_arr(2,maxunq*maxunq)
        integer begin

c temporary variables to store Aij, Bij, Qij.
        double precision temp_A, temp_B, temp_Q

c flag_pair is 0 if wconn does not write out pair type info
c flag_pair is 1 if wconn does not write out pair type info
        integer flag_pair
cccccccccccccccccccc huhnkie end ccccccccccccccccccccccccc

c hydrogen bonds
        character*1 hb
c end hydrogen bonds

        parameter (vrsn='12.1')

        data TILDA,name/'~','wconn'/
        data namel/5/

        call init_var()

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

        write(uwcon,99)vrsn
99      format(a4)
        write(uwcon,100)TILDA
100     format(a,' CONNECTIVITY FILE FOR MOLECULES: ')

c------------------------------------------------------
        if (specl) then
       write(uwcon,1001)
1001   format('~ totmon npt nb nangl totex totspe lestyp NBULK  ')
       write(uwcon,1002)totmon,npt,nb,nangl,totex,totspe,lestyp,NBULK
1002   format(1x,8i8)
         write(uwcon,51)('snb(',i,') ',i=1,nmb),('sang(',i,') '
     1,i=1,nmb)
51       format('~',8(a5,i1,a2))
       write(uwcon,5550)(snb(i),i=1,nmb),(snang(i),i=1,nmb)
5550     format(1x,8i8)

        lz14(0) = 0
        write(uwcon,1003)(lz14(i),i=0,nmb)
1003    format('~ Pointers to landau-Zener 14 ',/,5(1x,i7))
        
      write(uwcon,136)(poitype(i),i=1,nmb)
136   format('~ Pointers to last particle of each special unit'
     1   ,/,10(1x,i7))
c-----------------------------------------------------

      else
         write(uwcon,101)
101      format('~ totmon  npt    nb nangl  ntors  nimp totex totspe ',
     1 'lestyp NBULK nmb nchg nbeta nwaters arith prll_on_off mdivyes')
         write(uwcon,102)totmon,npt,nb,nangl,ntors,nimp,totex,totspe,
     1   lestyp,NBULK,nmb,nchg,nbeta,nwaters,arith,prll_on_off,mdivyes
102      format(1x,14i8,3(l2))
       end if


        do 1 i=1,NBULK
         write(uwcon,103)BULK(I)
103      format(a)
1       continue
        write(uwcon,104)(pbulk(i),i=1,nbulk)
104     format('~ Pointers to last particle of BULK',/,10(1x,i7))
        write(uwcon,105)(moname(i),i=1,totmon)
105     format('~ Monomer names ',/,10(1x,a4))
        write(uwcon,106)(poipt(i),i=1,totmon)
106     format('~ Pointers to last particle of monomer',/,10(1x,i7))
         if(specl) then
        write(uwcon,710)
710     format('~ Properties of particles list : ',/,
     1'~ pt mono grupid ptid lesid mutaid ptnm ptms ptchg epsgm6 epsgm12
     2     ptwei ptsaid')
        else
        write(uwcon,107)
107     format('~ Properties of particles list : ',/,
     1'~pt mono ptid lesid mutaid ptnm   ptms   ptchg   epsgm6 epsgm12 
     2   ptwei ptsaid hbond')
        end if


        do 2 i=1,npt

c---------------------------------------------------------------
       if (specl) then
        k = newpoit(i)
        write(uwcon,977)i,k,poimon(i),ntyp(i),ptid(k),lesid(k),mutaid(k)
     1          ,ptnm(k),ptms(k),
     1                ptchg(k),epsgm6(k),epsgm12(k),ptwei(k),ptsaid(k)
977     format(1x,i3,1x,i7,1x,i3,1x,i5,1x,i3,1x,i3,1x,i3,1x,a4,1x,f7.3,
     1                              2x,f9.5,1x,e12.5,e12.5,e12.5,1x,i3)
         else
c hydrogen bonds
        hb=' '
        if (ptnm(i).eq.'O   ') hb='A'
        if (ptnm(i).eq.'H   ') hb='D'        
c end hydgrogen bonds
        write(uwcon,108)i,poimon(i),ptid(i),lesid(i),mutaid(i),ptnm(i)
     1                ,ptms(i),ptchg(i),epsgm6(i),epsgm12(i),ptwei(i),
     2                ptsaid(i),hb
108     format(1x,i7,1x,i7,1x,i3,1x,i3,1x,i3,1x,a4,1x,f7.3,2x,f9.5,1x,
     1          e12.5,
     1                e12.5,e10.3,1x,i3,1x,a1)
         end if
c--------------------------------------------------------------
2       continue

        if (nb.gt.0) then
         write(uwcon,109)
109      format('~ Bonds list: ',/,'~ ib1 ib2 kbond req ')
         do 3 i=1,nb
          write(uwcon,110)ib1(i),ib2(i),kbond(i),req(i)
110       format(2(1x,i7),2(1x,f10.4))
3        continue
        end if
         
c----------------------------------------------
        if (nmb.gt.0 .and. (.not.specl)) then
         write(uwcon,1000)
1000     format('~ morse bond list:',/,'~ imb1 imb2 rmeq ')
         do 2000 i=1,nmb
          write(uwcon,1500)imb1(i),imb2(i),rmeq(i),D(i),alpha(i)
1500      format(2(1x,i7),3(1x,f10.4))
2000     continue
        end if
c------------------------------------------------

        if (nangl.gt.0) then
          write(uwcon,111)
111       format('~ Angles list: ',/,
     1'~ iangl1 iangl2 iangl3 kangl angleq ')
          do 4 i=1,nangl
           write(uwcon,112)iangl1(i),iangl2(i),iangl3(i),
     1kangl(i),angleq(i)*pi180
112        format(3(1x,i7),2(1x,f10.5))
4         continue
        end if

       if (.not.specl) then

        if (ntors.gt.0) then
          write(uwcon,113)
113       format('~ Torsions list: ',/,
     1'~ itor1 itor2 itor3 itor4 period ktors1 ktors2
     2 ktors3 phase1 phase2 phase3 ')
          do 5 i=1,ntors
           write(uwcon,114)itor1(i),itor2(i),itor3(i),itor4(i),period(i)
     1      ,ktors1(i),ktors2(i),ktors3(i),phase1(i),phase2(i),phase3(i)
114        format(5(1x,i7),3(1x,f10.4),1x,3(f6.2))
5         continue
        end if


c       ileana

        if(.not.arith)then

        if (nimp.gt.0) then
          write(uwcon,115)
115       format('~ Improper torsion properties: ',/
     1'~iimp1 iimp2 iimp3 iimp4 kimp impeq ')
          do 6 i=1,nimp
           write(uwcon,116)iimp1(i),iimp2(i),iimp3(i),iimp4(i),
     1kimp(i),impeq(i)*pi180
116        format(4(1x,i7),2(1x,e15.8))
6         continue
        end if

        else

           if (nimp.gt.0) then
          write(uwcon,1115)
1115      format('~ Improper torsion properties: ',/
     1'~iimp1 iimp2 iimp3 iimp4 kimp ')
          do 1126 i=1,nimp
           write(uwcon,1116)iimp1(i),iimp2(i),iimp3(i),iimp4(i),
     1kimp(i)
1116       format(4(1x,i7),1x,e15.8)
1126      continue
        end if
           
        endif
       end if
        
        if (totex.gt.0) then
         write(uwcon,117)
117      format('~ Exclusion list 1-2 1-3 1-4, set as followed:',/
     1'~ atom number, number of exclusions and list')
         k = 0
         do 7 i=1,npt
          if (exc1(i)-k.gt.0) then
           write(uwcon,118)i,exc1(i)-k
118        format(2(1x,i7))
           write(uwcon,119)(exc2(j),j=k+1,exc1(i))
119        format(1x,10(i7,1x))
           k = exc1(i)
          end if
7        continue
        end if

        if (muta) then
         k = 0
         j = 0
         l = 0
         write(uwcon,230)
230      format('~ Special exclusion list for mutants')
         do 237 i=1,npt

          if (mutaid(i).gt. 0) then
           l=l+1
           write(uwcon,118)i,exm1(l)-k
           write(uwcon,119)(exm2(j),j=k+1,exm1(l))
           k = exm1(l)
          endif

237      enddo
        endif
 
        if (totspe.gt.0) then

         if (specl) then

          write(uwcon,120)
          do 71 i=1,totspe
            j = spec1(i)
            k = spec2(i)
            write(uwcon,121)j,k,p14(1,i),p14(2,i),p14(3,i)
71        continue

         else
c NOT special

c       ileana

          
           write(uwcon,120)
120        format('~ Special list 1-4  set as followed:',/
     1'~ atom number, atom number, esp12*v14 eps6*v14 qiqj*el14 ')

         
C@
C        write(*,*)' v14 el14 ',v14,el14
           do 8 i=1,totspe
            j = spec1(i)
            k = spec2(i)

c       ileana

            if(.not.arith) then

           if ((epsgm12(j).lt.1.d-3 .or. epsgm12(k).lt.1.d-3)
     1          .and. hvdw0 ) then
                p14(1,i) = 0.d0
                p14(2,i) = 0.d0
           else

                p14(1,i) = epsgm12(j)*epsgm12(k)*v14
                p14(2,i) = epsgm6(j)*epsgm6(k)*v14

                
                
        endif   

        else
c               if ((epsgm6(j).lt.0.25 .or. epsgm6(k).lt.0.25)
c     1         .and. hvdw0 ) then
c               p14(1,i) = 0.d0
c               p14(2,i) = 0.d0
c               else
           
                p14(2,i) = 0.5d0*(epsgm6(j)+epsgm6(k))
                p14(2,i) = p14(2,i)*p14(2,i)*p14(2,i)
                p14(2,i) = p14(2,i)*p14(2,i)
                p14(1,i) = (epsgm12(j)*epsgm12(k))*p14(2,i)*p14(2,i)*v14
                p14(2,i) = (epsgm12(j)*epsgm12(k))*p14(2,i)*v14
c               endif 
           end if
           p14(3,i) = ptchg(j)*ptchg(k)*eps*el14*kofdie
c@
c          if (specl) then
c                       write(*,*)' j k p14(1,i) p14(3,i) '
c                       write(*,*)j,k,p14(1,i),p14(3,i)
c                       write(*,*)' epsgm12(j) epsgm12(k) '
c                       write(*,*)epsgm12(j),epsgm12(k)
c          end if
           if ( (lesid(j).ne.0) .and. (lesid(k).eq.lesid(j)) ) then
                tmp = 1.d0/ptwei(j)
                p14(1,i) = p14(1,i)*tmp
                p14(2,i) = p14(2,i)*tmp
                p14(3,i) = p14(3,i)*tmp
           end if
         



           
           write(uwcon,121)j,k,p14(1,i),p14(2,i),p14(3,i)
c121       format(2(1x,i5),3(1x,f15.5))
121       format(2(1x,i7),1x,f30.5,2(1x,f15.5))

        

             

          

8        continue

        endif
        endif


c write down charge flags:
c flagchr is a logical vector.
c if flagchr(i) is true the atom is charged.
       if ( .not. specl ) then
         write(uwcon,123)
123      format('~ Charge flags. T = charged atom , F = uncharged')
        write(uwcon,124)(flagchr(i),i=1,npt)
124     format(1x,15l2)

c chen
c write down les copy lables:
c cplbl is an integer vector.
         write(uwcon,125)
125      format('~ LES copy labels ')
         write(uwcon,126)(cplbl(i),i=1,npt)
126     format(1x,15i6)

c
c write down hydrophobic forces arrays
        if (nbeta.ne.0) then
         write(uwcon,127)
127      format('~ Hydrophobic forces parameters ',/
     1          '~ first pointers to Cb, then constants ')
         write(uwcon,128)(betap(i),i=1,nbeta)
128      format(1x,10i8)
         write(uwcon,129)(cbeta(i),i=1,nbeta)
129      format(1x,10f6.3)
         write(uwcon,130)avg_hydro
130      format(1x,2f9.4)
        end if
       end if
c yael
c write the monomers division data
         if (mdivyes) then
           write(uwcon,*)'~ the monomers division data'
           do 200 i=1,totmon
           write(uwcon,210)(mdivlist(0,i))
           if (mdivlist(0,i).eq.-1) then
                write(uwcon,210)(mdivlist(j,i),j=1,2)
           else
                 write(uwcon,210)(mdivlist(j,i),j=1,mdivlist(0,i)+1)
           endif

200       continue
210       format(1x,7i9)
         endif

c write down monomer numbers of TIP3 waters
c
         if (nwaters.gt.0) then
          write(uwcon,*)
     1'~ The monomers divisions below are TIP3/SPC/SPCE waters'
          write(uwcon,211)(idxtip3(i),i=1,nwaters)
211       format(10i8)
         end if


c hk add start
         if (flag_pair .eq. 1) then
c hk add end

ccccccccccccccccccc huhnkie start cccccccccccccccccccccccccccccccccccc
c
c i and j are particle types (ptid).
c ii and jj are the particle identifiers (pt)

            index=0
            begin=0
            do 10 ii=1,npt   
               do 20 jj=ii+1,npt
                            
c       see if (ii,jj) is in the exclusion list. 
c       if so, do not store the pair type.     
c       begin=0
                  
                  do 30 k=begin+1,exc1(ii)
                     if(jj .eq. exc2(k)) then
                        go to 20
                     end if    
 30               continue

c       see if the pair type is created already in the type_arr.
c       if so, do not store the pair type.

              do 40 kk = 1, index
                 if( ( (type_arr(1,kk).eq.ptid(ii)).and.
     6                 (type_arr(2,kk).eq.ptid(jj))      )
     6                           .or.
     6               ( (type_arr(1,kk).eq.ptid(jj)).and.
     6                 (type_arr(2,kk).eq.ptid(ii))      )

     6             )           then
                    go to 20
                 end if
 40                         continue

              index = index + 1
              type_arr(1,index)=ptid(ii)
              type_arr(2,index)=ptid(jj)

              i = ptid(ii)
              j = ptid(jj)

              A_arr(i,j) = 0.0d0
              B_arr(i,j) = 0.0d0
              Q_arr(i,j) = 0.0d0

              A_arr(i,j) = epsgm12(ii)*epsgm12(jj)
              B_arr(i,j) = epsgm6(ii)*epsgm6(jj)
              Q_arr(i,j) = ptchg(ii)*ptchg(jj)*kofdie/eps

 20                   continue

           begin = exc1(i)

 10             continue

c num_pairs is a global variable declared in CONNECT.BLOCK
c print out num_pairs

        num_pairs = index
        write(uwcon,183)
 183        format('~ the total number of pair types')
        write(uwcon, 74) num_pairs
 74          format(1x,i5)

c print out i, j, Aij, Bij, Qij, where i and j are particle types.
        write(uwcon,180)
 180        format('~ Pair type list : ',/,
     6  '~ pt1_type  pt2_type  Aij Bij Qij')

        do 75 mm=1, num_pairs
           i = type_arr(1,mm)
           j = type_arr(2,mm)
           write(uwcon,280) i,j,
     6          A_arr(i,j),B_arr(i,j),Q_arr(i,j)
 280              format(1x,i5,1x,i5,1x,e14.7,1x,e14.7,1x,e14.7)

 75                       continue

ccccccccccccccccccccc huhnkie end cccccccccccccccccccccccccccccccccc

c hk add start
        endif
c hk add end


        return
        end


