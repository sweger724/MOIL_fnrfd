      subroutine d2vanlyt(d2v,eigen)
c     
c copy from the special packing of the non bonded second derivative 
c matrix to the regular 3N*3N second derivative matrix
c     
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c     
      character*6 name
      integer namel,i,nloops,iloop
      integer pcounter
      double precision tmp,a1,b1,q
      double precision pick 
      double precision epstmp
      integer j,k,jbeg,jend,symidx,istart,ireal,inter
      logical eigen,sepfast
c
      integer kk,jj,ii,ll,mm,nn
      integer ith,iith,iphi,iiphi,k1
c      
      double precision sqtmp
      double precision d2v(3*maxpt2d,3*maxpt2d)
      
c vectors needed for the digonalization of the matrix
c
      double precision eigenv(3*maxpt),work(3*maxpt)
c
      double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg

      integer jbeg1,jend1,jbeg2,jend2
      integer io,ih1,ih2,jo,jh1,jh2,iwater

      data oosgm6,oosgm12,ohchrg,oochrg,hhchrg/
     1		595.05497170775496395324,
     2		582002.66166028445200793331,
     3		-115.4871969048,230.9743938096,
     4		57.7435984524/
c
      sepfast=.false.
c
c.....the diagonal values are acumulated and used only at the end of the
c.....calculation.
      do 10 i=1,6*npt
         diag(i) = 0.d0
 10   continue
c
      call nbondm()
      if (esymyes) call syminit()

c
c
c.....bond contributions................................................
      if (ebyes) then
c
         if (debug) write(*,*)' Calling bond2 '
c
c........nbond uses 6*dimension=6*9*maxptoffd
         nloops=nb/(9*maxptoffd)+1
c
         do 20 iloop=1,nloops
c
            ifirst=(iloop-1)*9*maxptoffd+1
            ilast=iloop*9*maxptoffd
            ilast=min(ilast,nb)
c
            call bond2_pcs()
c
            do 120 k=ifirst,ilast
               ii=3*(ib1(k)-1)+1
               jj=3*(ib2(k)-1)+1
               kk = 6*(k-ifirst)+1
c
               call convpair(d2v,ii,jj,kk)
c
 120        continue
c
 20      continue
c
      end if
c
c.....bond-angle contributions..........................................
      if (ethyes) then
c
         if (debug) write(*,*)' calling theta2 '
c
c........ntheta uses 27*dimension=27*2*maxptoffd
         nloops=nangl/(2*maxptoffd)+1
c
         do 30 iloop=1,nloops
c
            ifirst=(iloop-1)*2*maxptoffd+1
            ilast=iloop*2*maxptoffd
            ilast=min(ilast,nangl)
c
            call theta2_pcs()
c
            do 130 ith=ifirst,ilast
c
               iith= 27*(ith-ifirst)+1
               ii  = 3*(iangl1(ith)-1)+1
               jj  = 3*(iangl2(ith)-1)+1
               kk  = 3*(iangl3(ith)-1)+1
c pair i,j
               call convnotpair(d2v,ii,jj,iith)
c pair i,k
               iith=iith+9
               call convnotpair(d2v,ii,kk,iith)
c pair j,k
               iith=iith+9
               call convnotpair(d2v,jj,kk,iith)
 130        continue
c
 30      continue
c
      end if
c
c.....torsion contributions.............................................
      if (etoyes) then
c
         if (debug) write(*,*)' calling etors2 '
c
c........ntors uses 54*dimension=54*maxptoffd
         nloops=ntors/(maxptoffd)+1
c
         do 40 iloop=1,nloops
c
            ifirst=(iloop-1)*maxptoffd+1
            ilast=iloop*maxptoffd
            ilast=min(ilast,ntors)
c
            call etors2_pcs()

            do 140 iphi=ifirst,ilast
c
               iiphi = 54*(iphi-ifirst)+1
               kk    = 3*(itor1(iphi)-1)+1
               ll    = 3*(itor2(iphi)-1)+1
               mm    = 3*(itor3(iphi)-1)+1
               nn    = 3*(itor4(iphi)-1)+1
c..............pair k,l (note the opposite order in the call)
               call convnotpair(d2v,ll,kk,iiphi)
c..............pair k,m
               iiphi=iiphi+9
               call convnotpair(d2v,mm,kk,iiphi)
c..............pair k,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,kk,iiphi)
c..............pair l,m
               iiphi=iiphi+9
               call convnotpair(d2v,mm,ll,iiphi)
c..............pair l,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,ll,iiphi)
c..............pair m,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,mm,iiphi)
c
 140        continue
c
 40      continue

      end if
c
c.....improper torsion contributions....................................
      if (eimyes) then 
c
         if (debug) write(*,*)' calling eimphi2 '
c
c........nimptors uses 54*dimension=54*maxptoffd
         nloops=nimp/(maxptoffd)+1
c
         do 50 iloop=1,nloops
c
            ifirst=(iloop-1)*maxptoffd+1
            ilast=iloop*maxptoffd
            ilast=min(ilast,nimp)
c    
            call eimphi2_pcs()
c
            do 150 iphi=ifirst,ilast
               
               kk = 3*(iimp1(iphi)-1)+1
               ll = 3*(iimp2(iphi)-1)+1
               mm = 3*(iimp3(iphi)-1)+1
               nn = 3*(iimp4(iphi)-1)+1
               iiphi = 54*(iphi-ifirst)+1
c..............pair k,l (note the opposite order in the call)
               call convnotpair(d2v,ll,kk,iiphi)
c..............pair k,m
               iiphi=iiphi+9
               call convnotpair(d2v,mm,kk,iiphi)
c..............pair k,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,kk,iiphi)
c..............pair l,m
               iiphi=iiphi+9
               call convnotpair(d2v,mm,ll,iiphi)
c..............pair l,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,ll,iiphi)
c..............pair m,n
               iiphi=iiphi+9
               call convnotpair(d2v,nn,mm,iiphi)
c
 150        continue
 50      continue
c
      end if

c
c.....non-bonded contributions (divided in 3 lists).....................
      if (evdyes.or.eelyes) then
c
         if (eelyes) then
            epstmp = 332.0716d0/eps
          else
            epstmp =  0.0d0
         end if
c
         if (evdyes) then
            pick = 1.0d0
          else
            pick = 0.0d0
         end if

c
c....... contributions from the first list (Elec. and vdW)..............
         if (debug) write(*,*)' calling cdie2_pcs1'
c
         jbeg=1
         pcounter=0
c
         if (debug) write(*,*) 'first non-bonded list'

         do 60 i=1,npt-1
c
            jend = point1(i)
c
            if (jbeg.le.jend) then
c
                tmp = 1.d0/ptwei(i)
                a1   = epsgm12(i)*pick
                b1   = epsgm6(i)*pick
                q   = ptchg(i)*epstmp
c
                do 600 k=jbeg,jend
c
                    j=list1(k)
c
                    pcounter=pcounter+1
                    pairs1(pcounter)=i
                    pairs2(pcounter)=j
                    if (debug) write(*,*) i,j
                    a_pair(pcounter)= a1*epsgm12(j)
                    b_pair(pcounter)= b1*epsgm6(j)
                    q_pair(pcounter)= q*ptchg(j)
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
                       call cdie2_pcs1(sepfast)
c
                       do 160 k1=ifirst,ilast
                          kk = 6*(k1-ifirst)+1
                          ii = 3*(pairs1(k1)-1)+1
                          jj = 3*(pairs2(k1)-1)+1 
                          call convpair(d2v,ii,jj,kk)
 160                   continue
c
                       pcounter=0
c
                    end if
c
 600            continue
c
            end if 
c
            jbeg = jend + 1
c
 60      continue
c
         if (pcounter.lt.(9*maxptoffd)) then
c
            ifirst=1
            ilast=pcounter
            call cdie2_pcs1(sepfast)
c
            do 260 k1=ifirst,ilast
               kk = 6*(k1-ifirst)+1
               ii = 3*(pairs1(k1)-1)+1
               jj = 3*(pairs2(k1)-1)+1 
               call convpair(d2v,ii,jj,kk)
 260        continue
c
         end if

c....... contributions from the second list (vdW only)..................
         if (debug) write(*,*)' calling cdie2_pcs2'

         jbeg=1
         pcounter=0
c
         if (debug) write(*,*) 'second non-bonded list'

         do 70 i=1,npt-1
c
            jend = point2(i)
c
            if (jbeg.le.jend) then
c
                tmp = 1.d0/ptwei(i)
                a1   = epsgm12(i)*pick
                b1   = epsgm6(i)*pick
c
                do 700 k=jbeg,jend
c
                    j=list2(k)
c
                    pcounter=pcounter+1
                    pairs1(pcounter)=i
                    pairs2(pcounter)=j
                    if (debug) write(*,*) i,j
                    a_pair(pcounter)= a1*epsgm12(j)
                    b_pair(pcounter)= b1*epsgm6(j)
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
                       call cdie2_pcs2(sepfast)
c
                       do 170 k1=ifirst,ilast
                          kk = 6*(k1-ifirst)+1
                          ii = 3*(pairs1(k1)-1)+1
                          jj = 3*(pairs2(k1)-1)+1 
                          call convpair(d2v,ii,jj,kk)
 170                   continue
c
                       pcounter=0
c
                    end if
c
 700            continue
c
            end if 
c
            jbeg = jend + 1
c
 70      continue
c
         if (pcounter.lt.(9*maxptoffd)) then
c
            ifirst=1
            ilast=pcounter
            call cdie2_pcs2(sepfast)
c
            do 270 k1=ifirst,ilast
               kk = 6*(k1-ifirst)+1
               ii = 3*(pairs1(k1)-1)+1
               jj = 3*(pairs2(k1)-1)+1 
               call convpair(d2v,ii,jj,kk)
 270        continue
c
         end if
c
c
c
c....... contributions from the third list (Elet. only).................
         if (debug) write(*,*)' calling cdie2_pcs3'
c
         jbeg=1
         pcounter=0
c
         if (debug) write(*,*) 'third non-bonded list'

         do 80 i=1,npt-1
c
            jend = point3(i)
c
            if (jbeg.le.jend) then
c
                tmp = 1.d0/ptwei(i)
                q   = ptchg(i)*epstmp
c
                do 800 k=jbeg,jend
c
                    j=list3(k)
c
                    pcounter=pcounter+1
                    pairs1(pcounter)=i
                    pairs2(pcounter)=j
                    if (debug) write(*,*) i,j
                    q_pair(pcounter)= q*ptchg(j)
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
                       call cdie2_pcs3(sepfast)
c
                       do 180 k1=ifirst,ilast
                          kk = 6*(k1-ifirst)+1
                          ii = 3*(pairs1(k1)-1)+1
                          jj = 3*(pairs2(k1)-1)+1 
                          call convpair(d2v,ii,jj,kk)
 180                   continue
c
                       pcounter=0
c
                    end if
c
 800            continue
c
            end if 
c
            jbeg = jend + 1
c
 80      continue
c
         if (pcounter.lt.(9*maxptoffd)) then
c
            ifirst=1
            ilast=pcounter
            call cdie2_pcs3(sepfast)
c
            do 280 k1=ifirst,ilast
               kk = 6*(k1-ifirst)+1
               ii = 3*(pairs1(k1)-1)+1
               jj = 3*(pairs2(k1)-1)+1 
               call convpair(d2v,ii,jj,kk)
 280        continue
c
         end if

c
c        start wat-wat real contributions
c
	 if (nwaters.gt.1) then

            oosgm12=oosgm12*pick
            oosgm6=oosgm6*pick
            if (.not.(eelyes)) then
               oochrg= epstmp
               ohchrg= epstmp
               hhchrg= epstmp
            end if
            pcounter=0

            do 405 i=1,nwaters-1

		jbeg1  = poinwat1(i)+1
		jbeg2  = poinwat2(i)+1
		jend1  = poinwat1(i+1)
		jend2  = poinwat2(i+1)

		io = dpoipt(indxwat(i))-2
		ih1=io+1
		ih2=io+2

		if (jbeg1.le.jend1) then
                write(6,*)'first ww list jbeg1 jend1 ',jbeg1,jend1
		do 205 k=jbeg1,jend1
                   write(6,*)'pcounter ',pcounter
c
                    if (pcounter+9.ge.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=pcounter
c
                       call cdie2_wat1()
c
                       do 175 k1=ifirst,ilast
                          kk = 6*(k1-ifirst)+1
                          ii = 3*(pairs1(k1)-1)+1
                          jj = 3*(pairs2(k1)-1)+1 
                          call convpair(d2v,ii,jj,kk)
 175                   continue
c
                       pcounter=0
c
                    end if
c
                    j  = listwt1(k)

                    jo = dpoipt(j)-2
                    jh1= jo + 1
                    jh2= jo + 2

c Oxygen - Oxygen part
                    if (debug) write(*,*)'o_wat-o_wat ', io,jo
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jo
                    a_pair(pcounter)= oosgm12
                    b_pair(pcounter)= oosgm6
                    q_pair(pcounter)= oochrg
                    tmp_pair(pcounter)= 1

c electrostatic interaction of Oi1-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0
		
c electrostatic interaction of Oi1-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Oj
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jo
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Oj
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jo
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0


 205            continue
		end if

c
                if (pcounter.lt.(9*maxptoffd)) then
c
                   ifirst=1
                   ilast=pcounter
                   call cdie2_wat1()
c
                   do 176 k1=ifirst,ilast
                      kk = 6*(k1-ifirst)+1
                      write(6,*)'conv i j ',pairs1(k1),pairs2(k1)
                      ii = 3*(pairs1(k1)-1)+1
                      jj = 3*(pairs2(k1)-1)+1 
                      call convpair(d2v,ii,jj,kk)
 176               continue  
                end if

c start second loop including particles with upper cutoff
c cutele2 -  includes ONLY electrostic forces

c if distance of O-O larger than buffer cutoff
c exclude the whole water-water interaction
                
                pcounter=0
                write(6,*)'scnd ww list jbeg2 jend2 ',jbeg2,jend2
                if (jbeg2.le.jend2) then

                do 215 k=jbeg2,jend2
                   write(6,*)'pcounter ',pcounter
c
                    if (pcounter+9.ge.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=pcounter
c
                       call cdie2_wat2()
c
                       do 185 k1=ifirst,ilast
                          kk = 6*(k1-ifirst)+1
                          ii = 3*(pairs1(k1)-1)+1
                          jj = 3*(pairs2(k1)-1)+1 
                          call convpair(d2v,ii,jj,kk)
 185                   continue
c
                       pcounter=0
c
                    end if
c
                    j  = listwt2(k)
                    jo = dpoipt(j)-2
                    jh1= jo + 1
                    jh2= jo + 2

c Oxygen - Oxygen part
                    if (debug) write(*,*)'o_wat-o_wat ', io,jo
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jo
                    q_pair(pcounter)= oochrg
                    tmp_pair(pcounter)= 1

c electrostatic interaction of Oi1-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0
		
c electrostatic interaction of Oi1-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=io
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Oj
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jo
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Oj
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jo
                    q_pair(pcounter)= ohchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih1
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj1
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jh1
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj2
                    pcounter=pcounter+1
                    pairs1(pcounter)=ih2
                    pairs2(pcounter)=jh2
                    q_pair(pcounter)= hhchrg
                    tmp_pair(pcounter)= 0


 215            continue
                end if

c
                if (pcounter.lt.(9*maxptoffd)) then
c
                   ifirst=1
                   ilast=pcounter
                   call cdie2_wat2()
c
                   do 186 k1=ifirst,ilast
                      kk = 6*(k1-ifirst)+1
                      ii = 3*(pairs1(k1)-1)+1
                      jj = 3*(pairs2(k1)-1)+1 
                      call convpair(d2v,ii,jj,kk)
 186               continue  
                end if



 405       continue

         end if
c        end wat-wat

c
c        end real (including wat-wat) contributions

c
c        start symmetry contributiuons

         if (esymyes) then
      write(6,*) ' symmetry ',iblock1(symnum),symnum
      write(6,*) ' symmetry ',iblock2(symnum),iblock3(symnum)
            if (.not.(iblock1(symnum).gt.0 .or.
     1      iblock2(symnum).gt.0 .or. iblock3(symnum).gt.0))
     2      go to 110

            if (debug) write(*,*) 'symmetry contributions'

            tmp = 1.d0
            usesym = .true.
c
c           This first loop is up to the short (vdW) radius

            istart = 1
            jbeg   = point1(npt-1) + 1
            pcounter = 0

            if (debug) write(*,*) 'first non-bonded list'

c           103 is a loop on all symmetry operations
            do 103 symidx = 1,symnum
c              102 is a loop on all atoms that belong to the current symmetry
c              operation
               do 102 i=istart,iblock1(symidx)
                  jend = psym1(i)
                  ireal = symreal1(i)
                  if (jbeg.le.jend) then

                     a1=epsgm12(ireal)*pick
                     b1=epsgm6(ireal)*pick
                     q=ptchg(ireal)*epstmp

                     do 101 k=jbeg,jend

                        j=list1(k)
                        if (j.eq.ireal) go to 101
c
                        pcounter=pcounter+1
                        pairs1(pcounter)=ireal
                        pairs1s(pcounter)=symidx
                        pairs2(pcounter)=j
                        a_pair(pcounter)= a1*epsgm12(j)
                        b_pair(pcounter)= b1*epsgm6(j)
                        q_pair(pcounter)= q*ptchg(j)
                        if (debug) write(*,*) ireal,j
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs1(sepfast)

                           do 161 k1=ifirst,ilast
                              kk = 6*(k1-ifirst)+1
                              ii = 3*(pairs1(k1)-1)+1
                              jj = 3*(pairs2(k1)-1)+1 
                              call convpair(d2v,ii,jj,kk)
 161                       continue
c                           call matvec_die(vector,vectmp)
c
                           pcounter=0
c
                        end if
c
 101                 continue

                  end if

                  jbeg   = jend + 1

 102           continue

               istart = iblock1(symidx) + 1

 103        continue

         if (pcounter.lt.(9*maxptoffd)) then
c
            ifirst=1
            ilast=pcounter
            call cdie2_pcs1(sepfast)
c
            do 261 k1=ifirst,ilast
               kk = 6*(k1-ifirst)+1
               ii = 3*(pairs1(k1)-1)+1
               jj = 3*(pairs2(k1)-1)+1
               call convpair(d2v,ii,jj,kk)
 261        continue

c
         end if

c
c This second loop deals with vdW of uncharged particles

            istart = 1
            jbeg   = point2(npt-1) + 1
            pcounter = 0

            if (debug) write(*,*) 'second non-bonded list'

c           203 is a loop on all symmetry operation
            do 203 symidx = 1,symnum
c              202 is a loop on all atoms that belong to the current symmetry
c              operation
               do 202 i=istart,iblock2(symidx)
                  jend = psym2(i)
                  ireal = symreal2(i)
                  if (jbeg.le.jend) then

                     a1=epsgm12(ireal)*pick
                     b1=epsgm6(ireal)*pick

                     do 201 k=jbeg,jend

                        j=list2(k)
                        if (j.eq.ireal) go to 201
c
                        pcounter=pcounter+1
                        pairs1(pcounter)=ireal
                        pairs1s(pcounter)=symidx
                        pairs2(pcounter)=j
                        if (debug) write(*,*) ireal,j
                        a_pair(pcounter)= a1*epsgm12(j)
                        b_pair(pcounter)= b1*epsgm6(j)
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs2(sepfast)

                           do 171 k1=ifirst,ilast
                              kk = 6*(k1-ifirst)+1
                              ii = 3*(pairs1(k1)-1)+1
                              jj = 3*(pairs2(k1)-1)+1
                              call convpair(d2v,ii,jj,kk)
 171                       continue
c                           call matvec_die(vector,vectmp)
c
                           pcounter=0
c
                        end if
c
 201                 continue

                  end if

                  jbeg   = jend + 1

 202           continue

               istart = iblock2(symidx) + 1

 203        continue
c
         if (pcounter.lt.(9*maxptoffd)) then
c
            ifirst=1
            ilast=pcounter
            call cdie2_pcs2(sepfast)
c
            do 271 k1=ifirst,ilast
               kk = 6*(k1-ifirst)+1
               ii = 3*(pairs1(k1)-1)+1
               jj = 3*(pairs2(k1)-1)+1
               call convpair(d2v,ii,jj,kk)
 271        continue
c
         end if
c

c This third loop is on longer range electrostatic

            istart = 1
            jbeg   = point3(npt-1) + 1
            pcounter = 0

            if (debug) write(*,*) 'third non-bonded list'
c           303 is a loop on all symmetry operation
            do 303 symidx = 1,symnum
c              302 is a loop on all atoms that belong to the current symmetry 
c              operation
               do 302 i=istart,iblock3(symidx)
                  jend = psym3(i)
                  ireal = symreal3(i)
                  if (jbeg.le.jend) then

                     q=ptchg(ireal)*epstmp

                     do 301 k=jbeg,jend

                        j=list3(k)
                        if (j.eq.ireal) go to 301
c
                        pcounter=pcounter+1
                        pairs1(pcounter)=ireal
                        pairs1s(pcounter)=symidx
                        pairs2(pcounter)=j
                        if (debug) write(*,*) ireal,j
                        q_pair(pcounter)= q*ptchg(j)
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs3(sepfast)

                           do 181 k1=ifirst,ilast
                              kk = 6*(k1-ifirst)+1
                              ii = 3*(pairs1(k1)-1)+1
                              jj = 3*(pairs2(k1)-1)+1 
                              call convpair(d2v,ii,jj,kk)
 181                       continue
c                          call matvec_die(vector,vectmp)
c
                           pcounter=0
c
                        end if
c
 301                 continue

                  end if

                  jbeg   = jend + 1

 302           continue

               istart = iblock3(symidx) + 1

 303        continue


            if (pcounter.lt.(9*maxptoffd)) then
c
               ifirst=1
               ilast=pcounter
               call cdie2_pcs3(sepfast)

               do 281 k1=ifirst,ilast
                  kk = 6*(k1-ifirst)+1
                  ii = 3*(pairs1(k1)-1)+1
                  jj = 3*(pairs2(k1)-1)+1 
                  call convpair(d2v,ii,jj,kk)
 281           continue
c               call matvec_die(vector,vectmp)
c
            end if
c

c
c start wat-wat symmetry terms
c          
            if (.not.(iblckwt1(symnum).gt.0 .or. iblckwt2(symnum).gt.0)) 
     &      go to 110

            istart= 1
            jbeg1 = poinwat1(nwaters)+1
            pcounter = 0

            if (debug) write(6,*)' starting wat-wat symmetry terms '
c loop on symmetry operations
c
            do 407 symidx = 1,isym
               do 307 iwater=istart,iblckwt1(symidx)
                  jend1  = psymwt1(iwater)
                  i      = symmwat1(iwater)

                  io = dpoipt(i)-2
                  ih1=io+1
                  ih2=io+2

                  if (jbeg1.le.jend1) then

                  do 207 k=jbeg1,jend1
c
                     if (pcounter+9.ge.(9*maxptoffd)) then
c
                        ifirst=1
                        ilast=pcounter
c
                        call cdie2_wat1()
c
                        do 276 k1=ifirst,ilast
                           kk = 6*(k1-ifirst)+1
                           ii = 3*(pairs1(k1)-1)+1
                           jj = 3*(pairs2(k1)-1)+1 
                           call convpair(d2v,ii,jj,kk)
 276                    continue
c
                        pcounter=0
c
                     end if
c
                     j  = listwt1(k)
                     jo = dpoipt(j)-2
                     jh1= jo + 1
                     jh2= jo + 2

c Oxygen - Oxygen part
                     if (debug) write(*,*)'o_wat-o_wat ', io,jo
c skip the contribution from water own image
                     if (io.eq.jo) go to 667
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     a_pair(pcounter)= oosgm12
                     b_pair(pcounter)= oosgm6
                     q_pair(pcounter)= oochrg
                     tmp_pair(pcounter)= 1

c electrostatic interaction of Hi1-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

 667                 continue
c electrostatic interaction of Oi1-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0
		
c electrostatic interaction of Oi1-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Oj
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Oj
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

 207              continue
                  end if
                  jbeg1 = jend1 + 1
 307           continue
               istart = iblckwt1(symidx) + 1
 407        continue
c
            if (pcounter.lt.(9*maxptoffd)) then
c
               ifirst=1
               ilast=pcounter
               call cdie2_wat1()
c
               do 177 k1=ifirst,ilast
                  kk = 6*(k1-ifirst)+1
                  ii = 3*(pairs1(k1)-1)+1
                  jj = 3*(pairs2(k1)-1)+1 
                  call convpair(d2v,ii,jj,kk)
 177           continue  
            end if


c start second loop including particles with upper cutoff
c cutele2 - includes ONLY electrostic forces
c

            istart = 1
            jbeg2 = poinwat2(nwaters)+1
            pcounter = 0
            if (debug) write(6,*)' starting second symm wat-wat list'
c loop on symmetry operations
c
            do 707 symidx = 1,isym
c@	 write(*,*)' iblckwt2(symidx) = ',iblckwt2(symidx)
               do 607 iwater=istart,iblckwt2(symidx)
                  jend2  = psymwt2(iwater)
c@		write(*,*)' jbeg2 jend2 ',jbeg2,jend2
                  i      = symmwat2(iwater)

                  io = dpoipt(i)-2
                  ih1=io+1
                  ih2=io+2

                  if (jbeg2.le.jend2) then

                  do 507 k=jbeg2,jend2
c
                     if (pcounter+9.ge.(9*maxptoffd)) then
c
                        ifirst=1
                        ilast=pcounter
c
                        call cdie2_wat2()
c
                        do 178 k1=ifirst,ilast
                           kk = 6*(k1-ifirst)+1
                           ii = 3*(pairs1(k1)-1)+1
                           jj = 3*(pairs2(k1)-1)+1 
                           call convpair(d2v,ii,jj,kk)
 178                    continue
c
                        pcounter=0
c
                     end if
c
                     j  = listwt2(k)
c@			write(*,*)' #2 i j symidx ',i,j,symidx
                     jo = dpoipt(j)-2
                     jh1= jo + 1
                     jh2= jo + 2

c Oxygen - Oxygen part
                     if (debug) write(*,*)'o_wat-o_wat ', io,jo
                     if (io.eq.jo) go to 567
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     a_pair(pcounter)= oosgm12
                     b_pair(pcounter)= oosgm6
                     q_pair(pcounter)= oochrg
                     tmp_pair(pcounter)= 1

c electrostatic interaction of Hi1-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

 567                 continue
c electrostatic interaction of Oi1-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0
		
c electrostatic interaction of Oi1-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=io
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Oj
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Oj
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jo
                     q_pair(pcounter)= ohchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi1-Hj2
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih1
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh2
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

c electrostatic interaction of Hi2-Hj1
                     pcounter=pcounter+1
                     pairs1(pcounter)=ih2
                     pairs1s(pcounter)=symidx
                     pairs2(pcounter)=jh1
                     q_pair(pcounter)= hhchrg
                     tmp_pair(pcounter)= 0

 507              continue
                  end if
                  jbeg2 = jend2 + 1
 607           continue
               istart = iblckwt2(symidx) + 1
 707        continue
c
            if (pcounter.lt.(9*maxptoffd)) then
c
               ifirst=1
               ilast=pcounter
               call cdie2_wat2()
c
               do 179 k1=ifirst,ilast
                  kk = 6*(k1-ifirst)+1
                  ii = 3*(pairs1(k1)-1)+1
                  jj = 3*(pairs2(k1)-1)+1 
                  call convpair(d2v,ii,jj,kk)
 179           continue  
            end if

c           end of wat-wat symmetry contributions
            
            usesym = .false.

         end if  
c........end of the symmetry contributions
c
 110     continue
c

      end if
c     end of nonbonded elec. and vdw contributions

c
c
c
c.....special 1-4 bonded contributions..................................
      if (evdyes.or.eelyes) then
c
         if (debug) write(*,*)' calling cdie2_14'
c
c........special-14 uses 6*dimension=6*9*maxptoffd
         nloops=totspe/(9*maxptoffd)+1
c
         if (debug) write(*,*) 'special non-bonded list'

         do 90 iloop=1,nloops
c
            ifirst=(iloop-1)*9*maxptoffd+1
            ilast=iloop*9*maxptoffd
            ilast=min(ilast,totspe)
c
            call cdie2_14()
            do 190 k=ifirst,ilast

               kk = 6*(k-ifirst)+1
               ii = 3*(spec1(k)-1)+1
               jj = 3*(spec2(k)-1)+1 
               if (debug) write(*,*)spec1(k),spec2(k)
               call convpair(d2v,ii,jj,kk)
 190        continue
c
 90      continue
c
      end if
c
c
c.....commented out (left as in d2vdr2), because it is not used yet.
cout     if (nbeta.gt.0) then
cout         if (debug) write(*,*)' calling hyd2'
cout         call dhyd2()
cout---Added 12.01.94 by M. K-I for Hydrophobic Interactions--------
cout-----------New block for hydrophobicity-------------------------
cout---- ehyes is the flag for calculating hydrophobic interactions-
cout
cout         if (ehyes) then
cout            jj=0
cout            do 111 l=1,nbeta-1
cout               i=betap(l)
coutc              j = 6*(i-1)+1
cout               k = 3*(i-1)+1
cout               do 121 m=l+1,nbeta
cout                  kk=betap(m)
cout                  jj=jj+1
cout                  kkk=  3*(kk-1) + 1
cout                  jjj=  6*(jj-1) + 1
cout.element xi,xj which is the same as xj,xi (i,j, particle indices)   
cout                  d2v(k,kkk)     = d2hyd(jjj)  +  d2v(k,kkk)
cout                  d2v(kkk,k)     = d2hyd(jjj)  +  d2v(kkk,k)
coutc element xi,yj which is the same as yj,xi , xj,yi and yi,xj
coutc the last two are the same only for function of distances
cout                  d2v(k+1,kkk)   = d2hyd(jjj+1)  +   d2v(k+1,kkk)
cout                  d2v(kkk+1,k)   = d2hyd(jjj+1)  +   d2v(kkk+1,k)
cout                  d2v(k,kkk+1)   = d2hyd(jjj+1)  +   d2v(k,kkk+1)
cout                  d2v(kkk,k+1)   = d2hyd(jjj+1)  +   d2v(kkk,k+1)
coutc elements zi,xj xj,zi zj,xi xi,zj
cout                  d2v(k+2,kkk)   = d2hyd(jjj+2)  +   d2v(k+2,kkk)
cout                  d2v(kkk,k+2)   = d2hyd(jjj+2)  +   d2v(kkk,k+2)
cout                  d2v(kkk+2,k)   = d2hyd(jjj+2)  +   d2v(kkk+2,k)
cout                  d2v(k,kkk+2)   = d2hyd(jjj+2)  +   d2v(k,kkk+2)
coutc element yi,yj and yj,yi
cout                  d2v(k+1,kkk+1) = d2hyd(jjj+3)  +   d2v(k+1,kkk+1)
cout                  d2v(kkk+1,k+1) = d2hyd(jjj+3)  +   d2v(kkk+1,k+1)
coutc elements yi,zj zj,yi yj,zi zi,yj
cout                  d2v(k+1,kkk+2) = d2hyd(jjj+4)  +   d2v(k+1,kkk+2)
cout                  d2v(kkk+2,k+1) = d2hyd(jjj+4)  +   d2v(kkk+2,k+1)
cout                  d2v(kkk+1,k+2) = d2hyd(jjj+4)  +   d2v(kkk+1,k+2)
cout                  d2v(k+2,kkk+1) = d2hyd(jjj+4)  +   d2v(k+2,kkk+1)
cout elements zi,zj zj,zi
cout                  d2v(k+2,kkk+2) = d2hyd(jjj+5)  +   d2v(k+2,kkk+2)
cout                  d2v(kkk+2,k+2) = d2hyd(jjj+5)  +   d2v(kkk+2,k+2)
cout 121           continue
cout               
cout 111        continue
cout         endif
cout------------end of New block-------------------------------
cout
cout      end if
c
c      if (debug) call prd2v(stdo)
c
c
c
c    
      do 41 k=1,npt
         kk = 6*(k-1)+1
         ii = 3*(k-1)+1
         d2v(ii,ii)     = diag(kk)
         d2v(ii+1,ii)   = diag(kk+1)
         d2v(ii,ii+1)   = diag(kk+1)
         d2v(ii+2,ii)   = diag(kk+2)
         d2v(ii,ii+2)   = diag(kk+2)
         d2v(ii+1,ii+1) = diag(kk+3)
         d2v(ii+1,ii+2) = diag(kk+4)
         d2v(ii+2,ii+1) = diag(kk+4)
         d2v(ii+2,ii+2) = diag(kk+5)
 41   continue
c
c
      if (eigen) then

         do 31 i=0,3*npt-1
            do 311 j=0,3*npt-1
               ii = int(i/3) + 1 
               jj = int(j/3) + 1
               sqtmp = dsqrt(ptms(ii)*ptms(jj))
               d2v(i+1,j+1) = d2v(i+1,j+1)/sqtmp
 311        continue
 31      continue
c
         call house(d2v,3*npt,3*npt,eigenv,work,i)
         write(stdo,*)' Eigenvalues '
         write(stdo,100)(eigenv(j),j=1,3*npt)
 100     format(1x,3(f9.2,1x))
         write(stdo,*)' Error = ',i
      end if

      return
      end
