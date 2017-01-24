      subroutine getvec(vector,sepfast)
c
c calling specific second derivative routines, part 1("small matrices").
c
c LATER-invert the order of the calculations: first non-bonded. 
c This may be important in minimazing of roundoff errors.
c
c If sepfast=.true. skip wat-wat contributions and all offdiagonal terms
c of the type fast-slow (this allows shorter vectors! )
c It is simply assumed that slow modes are those associated with waters
c It is also assumed that all the particles, bonds, angles that belong
c to the slow part come first !!!!
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'

      include 'COMMON/EXTERNAL.BLOCK'
      include 'COMMON/SGB.BLOCK'
c
c     

      character*6 name
      integer namel,i,nloops,iloop
      double precision vector(3*maxpt),vectmp(3*maxpt)
c # this should be enough ??
c #      double precision vector(3*npts),vectmp(3*npts)
      integer pcounter,nbs,nangls
      double precision tmp,a1,b1,q,dpottemp(3,maxpt)
      double precision pick 
      double precision epstmp
      integer j,k,jbeg,jend,symidx,ireal,istart,inter,moildebug
      logical sepfast
c
c
c


      double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
      double precision e_gbsatmp

      integer jbeg1,jend1,jbeg2,jend2
      integer io,ih1,ih2,jo,jh1,jh2,iwater

      data oosgm6,oosgm12,ohchrg,oochrg,hhchrg/
     1		595.05497170775496395324,
     2		582002.66166028445200793331,
     3		-115.4871969048,230.9743938096,
     4		57.7435984524/
c

      data name/'getvec'/
      data namel/6/

c
c
c.....the diagonal values are acumulated and used only at the end of the
c.....calculation.
      if (.not.sepfast) npts = npt
      do 10 i=1,6*npts
         diag(i) = 0.d0
 10   continue
c
c.....the same with vectmp
      do 15 i=1,3*npts
         vectmp(i) = 0.d0
 15   continue

C     SGB 2nd derivatives

      
      if(gbsabool) then
         sndbool = 1


         moildebug = 0
         

         call egb_f(npt,ptchg,dielectric,sndbool
     $        ,vector,vectmp, debugsgb)

         
	sndbool = 0

	end if
        

      

c
c.....bond contributions................................................
      if (ebyes) then
c
         if (debug) write(*,*)'GETVEC: Calling bond2 '
c
c........nbond uses 6*dimension=6*9*maxptoffd
         if (sepfast) then
           nbs=nb-2*nwaters
           nloops=nbs/(9*maxptoffd)+1
c          in general case nbs should be defined on higher level
         else
           nloops=nb/(9*maxptoffd)+1
           nbs=nb
         end if
c
         do 20 iloop=1,nloops
c
            ifirst=(iloop-1)*9*maxptoffd+1
            ilast=iloop*9*maxptoffd
            ilast=min(ilast,nbs)
c
c            write (6,*) 'GETVEC: calling matvec bond stuff '
            call bond2_pcs()
            call matvec_bond(vector,vectmp)
c            write (6,*) 'GETVEC: calling bond' 

c
 20      continue
c
      end if


c
c.....bond-angle contributions..........................................
      if (ethyes) then
c
c         if (debug) write(*,*)'GETVEC: calling theta2 '
c
c........ntheta uses 27*dimension=27*2*maxptoffd
         if (sepfast) then
           nangls=nangl-nwaters
           nloops=nangls/(2*maxptoffd)+1
c          in general case nangls should be defined on higher level
         else
           nloops=nangl/(2*maxptoffd)+1
           nangls=nangl
         end if
c
         do 30 iloop=1,nloops
c
            ifirst=(iloop-1)*2*maxptoffd+1
            ilast=iloop*2*maxptoffd
            ilast=min(ilast,nangls)
c
c            write (6,*) 'GETVEC:calling angle bits '
            call theta2_pcs()
            call matvec_theta(vector,vectmp)
c
 30      continue
c
      end if

c     
c.....torsion contributions.............................................
      if (etoyes) then
c
c         if (debug) write(*,*) 'GETVEC:calling etors2 '
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

c            write (6,*) 'GETVEC:calling torsional bits '
            call etors2_pcs()
            call matvec_tors(vector,vectmp)
c
 40      continue

      end if


c
c.....improper torsion contributions....................................
      if (eimyes) then 
c
c         if (debug) write(*,*) 'GETVEC:calling eimphi2 '
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
c            write (6,*) 'GETVEC:calling improper bits '
            call eimphi2_pcs()
            call matvec_impto(vector,vectmp)
c
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

         usesym=.false.
c
c....... contributions from the first list (Elec. and vdW)..............
c         if (debug) write(*,*) 'GETVEC:calling cdie2_pcs1'
c
         jbeg=1
         pcounter=0
c
         if (sepfast) then
            nloops=npts
c           if there are no fast prts (wats) separation does not make sense
c           so npts is not the last prt by assumption
         else
            nloops=npt-1
         end if
c
         do 60 i=1,nloops
c
            jend = point1(i)
c
            if (jbeg.le.jend) then
c
                tmp = 1.d0/ptwei(i)
                a1   = epsgm12(i)*pick
                if (arith) then
                b1   = epsgm6(i)
                else
                b1   = epsgm6(i)*pick
                endif
                q   = ptchg(i)*epstmp
c
c
                do 600 k=jbeg,jend
c
                    j=list1(k)
c
                    pcounter=pcounter+1
                    pairs1(pcounter)=i
                    pairs2(pcounter)=j
                    a_pair(pcounter)= a1*epsgm12(j)
                    if (arith) then
                    b_pair(pcounter)= b1+epsgm6(j)
                    else
                    b_pair(pcounter)= b1*epsgm6(j)
                    endif
                    q_pair(pcounter)= q*ptchg(j)
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
c                       write (6,*) 'GETVEC:calling nonbond bits '
                       call cdie2_pcs1(sepfast)
                       call matvec_die(vector,vectmp,sepfast,npts)
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
c            write (6,*) 'GETVEC:nb2 '
            call cdie2_pcs1(sepfast)
            call matvec_die(vector,vectmp,sepfast,npts)
c
         end if


c....... contributions from the second list (vdW only)..................
         if (debug) write(*,*)' calling cdie2_pcs2'

         jbeg=1
         pcounter=0
c
         do 70 i=1,nloops
c
            jend = point2(i)
c
            if (jbeg.le.jend) then
c
                tmp = 1.d0/ptwei(i)
                a1   = epsgm12(i)*pick
                if (arith) then
                b1   = epsgm6(i)
                else
                b1   = epsgm6(i)*pick
                endif
c
                do 700 k=jbeg,jend
c
                    j=list2(k)
c
                    pcounter=pcounter+1
                    pairs1(pcounter)=i
                    pairs2(pcounter)=j
                    a_pair(pcounter)= a1*epsgm12(j)
                    if (arith) then
                    b_pair(pcounter)= b1+epsgm6(j)
                    else
                    b_pair(pcounter)= b1*epsgm6(j)
                    endif
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
                       call cdie2_pcs2(sepfast)
                       call matvec_die(vector,vectmp,sepfast,npts)
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
            call matvec_die(vector,vectmp,sepfast,npts)
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
         do 80 i=1,nloops
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
                    q_pair(pcounter)= q*ptchg(j)
                    tmp_pair(pcounter)= tmp
c
                    if (pcounter.eq.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=9*maxptoffd
c
                       call cdie2_pcs3(sepfast)
                       call matvec_die(vector,vectmp,sepfast,npts)
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
            call matvec_die(vector,vectmp,sepfast,npts)
c
         end if



c
c
c        start wat-wat real contributions
c
	 if ((nwaters.gt.1).and.(.not.(sepfast))) then

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

		do 205 k=jbeg1,jend1
c
                    if (pcounter+9.ge.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=pcounter
                       call cdie2_wat1()
                       call matvec_wat(vector,vectmp)
                       pcounter=0
c
                    end if
c
                    j  = listwt1(k)

                    jo = dpoipt(j)-2
                    jh1= jo + 1
                    jh2= jo + 2

c Oxygen - Oxygen part
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
                   call matvec_wat(vector,vectmp)
 176               continue  
                end if

c start second loop including particles with upper cutoff
c cutele2 -  includes ONLY electrostic forces

c if distance of O-O larger than buffer cutoff
c exclude the whole water-water interaction
                
                pcounter=0

                if (jbeg2.le.jend2) then

                do 215 k=jbeg2,jend2
c
                    if (pcounter+9.ge.(9*maxptoffd)) then
c
                       ifirst=1
                       ilast=pcounter
                       call cdie2_wat2()
                       call matvec_wat(vector,vectmp)
                       pcounter=0
c
                    end if
c
                    j  = listwt2(k)
                    jo = dpoipt(j)-2
                    jh1= jo + 1
                    jh2= jo + 2

c Oxygen - Oxygen part
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
                   ifirst=1
                   ilast=pcounter
                   call cdie2_wat2()
                   call matvec_wat(vector,vectmp) 
                end if



 405       continue

         end if
c        end wat-wat

c
c        end real (including wat-wat) contributions

c
c
c........symmetry contributions
c

         if (esymyes) then

            if (.not.(iblock1(symnum).gt.0 .or.
     1      iblock2(symnum).gt.0 .or. iblock3(symnum).gt.0))
     2      go to 110

            tmp = 1.d0
            usesym = .true.
c
c           This first loop is up to the short (vdW) radius

            istart = 1
            jbeg   = point1(npt-1) + 1
            pcounter = 0

c           103 is a loop on all symmetry operations
            do 103 symidx = 1,symnum
c              102 is a loop on all atoms that belong to the current symmetry
c              operation
               do 102 i=istart,iblock1(symidx)
                  jend = psym1(i)
                  ireal = symreal1(i)
                  if (jbeg.le.jend) then

                     a1=epsgm12(ireal)*pick
                     if (arith) then
                     b1=epsgm6(ireal)
                     else
                     b1=epsgm6(ireal)*pick
                     endif
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
                        if (arith) then
                        b_pair(pcounter)= b1+epsgm6(j)
                        else
                        b_pair(pcounter)= b1*epsgm6(j)
                        endif
                        q_pair(pcounter)= q*ptchg(j)
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs1(sepfast)
                           call matvec_die(vector,vectmp,
     &                          sepfast,npts)
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
             call matvec_die(vector,vectmp,sepfast,npts)
         end if

c
c This second loop deals with vdW of uncharged particles

            istart = 1
            jbeg   = point2(npt-1) + 1
            pcounter = 0

c           203 is a loop on all symmetry operation
            do 203 symidx = 1,symnum
c              202 is a loop on all atoms that belong to the current symmetry
c              operation
               do 202 i=istart,iblock2(symidx)
                  jend = psym2(i)
                  ireal = symreal2(i)
                  if (jbeg.le.jend) then

                     a1=epsgm12(ireal)*pick
                     if (arith) then
                     b1=epsgm6(ireal)
                     else
                     b1=epsgm6(ireal)*pick
                     endif

                     do 201 k=jbeg,jend

                        j=list2(k)
                        if (j.eq.ireal) go to 201
c
                        pcounter=pcounter+1
                        pairs1(pcounter)=ireal
                        pairs1s(pcounter)=symidx
                        pairs2(pcounter)=j
                        a_pair(pcounter)= a1*epsgm12(j)
                        if (arith) then
                        b_pair(pcounter)= b1+epsgm6(j)
                        else
                        b_pair(pcounter)= b1*epsgm6(j)
                        endif
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs2(sepfast)
                           call matvec_die(vector,vectmp,
     &                          sepfast,npts)
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
            call matvec_die(vector,vectmp,sepfast,npts)
         end if
c


c This third loop is on longer range electrostatic

            istart = 1
            jbeg   = point3(npt-1) + 1
            pcounter = 0

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
                        q_pair(pcounter)= q*ptchg(j)
                        tmp_pair(pcounter)= tmp
c
                        if (pcounter.eq.(9*maxptoffd)) then
c
                           ifirst=1
                           ilast=9*maxptoffd
c
                           call cdie2_pcs3(sepfast)
                           call matvec_die(vector,vectmp,
     &                          sepfast,npts)
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
               call matvec_die(vector,vectmp,sepfast,npts)
c
            end if
c
c
            if (sepfast) go to 110
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
                        call cdie2_wat1()
                        call matvec_wat(vector,vectmp)
                        pcounter=0
                     end if
c
                     j  = listwt1(k)
                     jo = dpoipt(j)-2
                     jh1= jo + 1
                     jh2= jo + 2

c Oxygen - Oxygen part
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
               call matvec_wat(vector,vectmp)
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
               do 607 iwater=istart,iblckwt2(symidx)
                  jend2  = psymwt2(iwater)
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
                        call cdie2_wat2()
                        call matvec_wat(vector,vectmp)
                        pcounter=0
                     end if
c
                     j  = listwt2(k)

                     jo = dpoipt(j)-2
                     jh1= jo + 1
                     jh2= jo + 2

c Oxygen - Oxygen part
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
               call matvec_wat(vector,vectmp)
            end if

c           end of wat-wat symmetry contributions
c            
            
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
c.....special 1-4 bonded contributions..................................
      if (evdyes.or.eelyes) then
c
         if (debug) write(*,*)' calling cdie2_14'
c
c........special-14 uses 6*dimension=6*9*maxptoffd
         nloops=totspe/(9*maxptoffd)+1
c
         do 90 iloop=1,nloops
c
            ifirst=(iloop-1)*9*maxptoffd+1
            ilast=iloop*9*maxptoffd
            ilast=min(ilast,totspe)
c
            call cdie2_14()
            call matvec_14(vector,vectmp)
c
 90      continue
c
      end if


c
c
c.....commented out (left as in d2vdr2), because it is not used yet.
cout      if (nbeta.gt.0) then
cout         if (debug) write(*,*)' calling hyd2'
cout         call dhyd2()
cout      end if
c
c      if (debug) call prd2v(stdo)
c
c
c
      call matvec_diag(vector,vectmp,npts)
c      write (6,*) 'GETVEC : DIAG '
c      do i = 1,npt
c         ii = (i-1)*3+1
c         write (6,*) vectmp(ii)-vectmpold(ii),
c     $        vectmp(ii+1)-vectmpold(ii+1),
c     $        vectmp(ii+2)-vectmpold(ii+2)
c         vectmpold(ii) = vectmp(ii)
c         vectmpold(ii+1) = vectmp(ii+1)
c         vectmpold(ii+2) = vectmp(ii+2)
c      end do
      
c
c  adding the contribution of the averaged second derivatives 
c
      if(sepfast) then
         do 1000 i=1,npts
            k=3*(i-1)
            vectmp(k+1)=vectmp(k+1)+vector(k+1)*d2pt_ave(1,i)+
     >           vector(k+2)*d2pt_ave(2,i)+vector(k+3)*d2pt_ave(3,i)
            vectmp(k+2)=vectmp(k+2)+vector(k+1)*d2pt_ave(4,i)+
     >           vector(k+2)*d2pt_ave(5,i)+vector(k+3)*d2pt_ave(6,i)
            vectmp(k+3)=vectmp(k+3)+vector(k+1)*d2pt_ave(7,i)+
     >           vector(k+2)*d2pt_ave(8,i)+vector(k+3)*d2pt_ave(9,i)
 1000    continue
      endif
c
c
      do 100 i=1,3*npts
         vector(i) = vectmp(i)
 100  continue
      return

      end

