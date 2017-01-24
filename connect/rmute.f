      subroutine rmute(ipick)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/MUTA.BLOCK'

      integer ipick(*)
      integer rmcnt,i,j,rmlink(maxtors),exch(maxtors)
      integer k,l,k1,kk,k2,k3
      integer mutcnt1,mutcnt2
      logical NEW

      NEW=.TRUE.
c maxtors is the bigger one between maxbond maxangl etc

      rmcnt=0

        if(nbd.gt.100) then
         write(*,*) 'nbd.gt.maxtors'
         stop
        endif

c -bonds

        do i=1,maxtors
         rmlink(i)=0
        enddo

        rmcnt=0

        if (nb.gt.0) then
         do 3 i=1,nb

c - case 1-2

          if ((ipick(ib1(i)).ne.ipick(ib2(i))) .and.
     1        (ipick(ib1(i)).ne.0) .and.
     2        (ipick(ib2(i)).ne.0) ) then
              rmcnt=rmcnt+1
              rmlink(rmcnt)=i
          endif

3        continue

         if (rmcnt.ne.0) then

          do i=1,rmcnt

           write(*,*) 'i,rmlink(i)',i,rmlink(i)
           if (rmlink(i).lt.nb) then

            do j=rmlink(i),nb-1
             ib1(j)=ib1(j+1)
             ib2(j)=ib2(j+1)
             kbond(j)=kbond(j+1)
             req(j)=req(j+1)
            enddo

            do j=i,rmcnt
             rmlink(j)=rmlink(j)-1
            enddo

            nb=nb-1

           else

            nb=nb-1

           endif
          enddo
       endif
      endif



c -angles

        do i=1,maxtors
         rmlink(i)=0
        enddo

        rmcnt=0

        if (nangl.gt.0) then
          do 4 i=1,nangl
c - case 1-2-?
           if     ((ipick(iangl1(i)).ne.ipick(iangl2(i))) .and.
     &              ipick(iangl1(i)).ne.0                 .and.
     &              ipick(iangl2(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case 1-?-2
           else if((ipick(iangl1(i)).ne.ipick(iangl3(i))) .and.
     &              ipick(iangl1(i)).ne.0                 .and.
     &              ipick(iangl3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-1-2
           else if((ipick(iangl2(i)).ne.ipick(iangl3(i))) .and.
     &              ipick(iangl2(i)).ne.0                 .and.
     &              ipick(iangl3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
          endif
4         continue
        if (rmcnt.gt.0) then
         do i=1,rmcnt
          write(*,*) 'i,rmlink(i)',i,rmlink(i)
          if (rmlink(i).lt.nangl) then
           do j=rmlink(i),nangl-1
            iangl1(j)=iangl1(j+1)
            iangl2(j)=iangl2(j+1)
            iangl3(j)=iangl3(j+1)
            kangl(j)=kangl(j+1)
            angleq(j)=angleq(j+1)
           enddo
           do j=i,rmcnt
            rmlink(j)=rmlink(j)-1
           enddo
           nangl=nangl-1 
          else
           nangl=nangl-1
          endif
         enddo
        end if
      endif

c -torsions

        do i=1,maxtors
         rmlink(i)=0
        enddo

        rmcnt=0

        if (ntors.gt.0) then
          do 5 i=1,ntors
c - case 1-2-?-?
           if     ((ipick(itor1(i)).ne.ipick(itor2(i))) .and.
     &              ipick(itor1(i)).ne.0                 .and.
     &              ipick(itor2(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case 1-?-2-?
           else if((ipick(itor1(i)).ne.ipick(itor3(i))) .and.
     &              ipick(itor1(i)).ne.0                 .and.
     &              ipick(itor3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case 1-?-?-2
           else if((ipick(itor1(i)).ne.ipick(itor4(i))) .and.
     &              ipick(itor1(i)).ne.0                 .and.
     &              ipick(itor4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-1-2-?
           else if((ipick(itor2(i)).ne.ipick(itor3(i))) .and.
     &              ipick(itor2(i)).ne.0                 .and.
     &              ipick(itor3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-1-?-2
           else if((ipick(itor2(i)).ne.ipick(itor4(i))) .and.
     &              ipick(itor2(i)).ne.0                 .and.
     &              ipick(itor4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-?-1-2
           else if((ipick(itor3(i)).ne.ipick(itor4(i))) .and.
     &              ipick(itor3(i)).ne.0                 .and.
     &              ipick(itor4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
          endif
5         continue
        if (rmcnt.gt.0) then
         do i=1,rmcnt
          if (rmlink(i).lt.ntors) then
           do j=rmlink(i),ntors-1
            itor1(j)=itor1(j+1)
            itor2(j)=itor2(j+1)
            itor3(j)=itor3(j+1)
            itor4(j)=itor4(j+1)
            period(j)=period(j+1)
            ktors1(j)=ktors1(j+1)
            ktors2(j)=ktors2(j+1)
            ktors3(j)=ktors3(j+1)
            phase1(j)=phase1(j+1)
            phase2(j)=phase2(j+1)
            phase3(j)=phase3(j+1)
           enddo
           do j=i,rmcnt
            rmlink(j)=rmlink(j)-1
           enddo
           ntors=ntors-1
          else
           ntors=ntors-1
          endif
         enddo
        end if
      endif

c - improper torsions

        do i=1,maxtors
         rmlink(i)=0
        enddo

        rmcnt=0

        if (nimp.gt.0) then
          do 6 i=1,nimp
c - case 1-2-?-?
           if     ((ipick(iimp1(i)).ne.ipick(iimp2(i))) .and.
     &              ipick(iimp1(i)).ne.0                 .and.
     &              ipick(iimp2(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case 1-?-2-?
           else if((ipick(iimp1(i)).ne.ipick(iimp3(i))) .and.
     &              ipick(iimp1(i)).ne.0                 .and.
     &              ipick(iimp3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case 1-?-?-2
           else if((ipick(iimp1(i)).ne.ipick(iimp4(i))) .and.
     &              ipick(iimp1(i)).ne.0                 .and.
     &              ipick(iimp4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-1-2-?
           else if((ipick(iimp2(i)).ne.ipick(iimp3(i))) .and.
     &              ipick(iimp2(i)).ne.0                 .and.
     &              ipick(iimp3(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-1-?-2
           else if((ipick(iimp2(i)).ne.ipick(iimp4(i))) .and.
     &              ipick(iimp2(i)).ne.0                 .and.
     &              ipick(iimp4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
c - case ?-?-1-2
           else if((ipick(iimp3(i)).ne.ipick(iimp4(i))) .and.
     &              ipick(iimp3(i)).ne.0                 .and.
     &              ipick(iimp4(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
          endif
6         continue
        if (rmcnt.gt.0) then
         do i=1,rmcnt
          if (rmlink(i).lt.nimp) then
           do j=rmlink(i),nimp-1
            iimp1(j)=iimp1(j+1)
            iimp2(j)=iimp2(j+1)
            iimp3(j)=iimp3(j+1)
            iimp4(j)=iimp4(j+1)
            period(j)=period(j+1)
            kimp(j)=kimp(j+1)
           enddo
           do j=i,rmcnt
            rmlink(j)=rmlink(j)-1
           enddo
           nimp=nimp-1
          else
           nimp=nimp-1
          endif
         enddo
        end if
      endif

c - add atoms of the mutant to the exclusion list of the original, and viceversa

c - 2/10 now a new version is created. All the atoms of the mutant and the native
c        sidechain are added to the exclusion list of all the other atoms, so are
c        the other atoms to the exlusion list of mutants.
c        A new exclusion list is created for the mutants and the native sidechain.
c        This new list is used in the calculation of non bonded interaction within
c        mutants and within natives and between those sidechains and the rest of the
c        protein.

        if (totex.gt.0) then
         k=0
         do 2001 i=1,npt
          if (ipick(i).ne.0) then
           do l=i+1,npt
            if ((ipick(l).ne.0).and.(ipick(l).ne.ipick(i))) then
             do j=k+1,exc1(i)
              if (exc2(j).eq.l) then
               go to 2323
              endif
             enddo
c take the next items in exc1 and copy it on a new array
             do j=k+1,totex
              exch(j+1) = exc2(j)
             enddo
c add the new member
             exch(k+1) = l
c increase by one the size of the exclusions list
             totex = totex+1
             DO j=i,npt
              exc1(j) = exc1(j)+1
             ENDDO
c recreate the new exclusion list
             do j=k+1,totex
              exc2(j) = exch(j)
             enddo
2323     continue
            endif
           enddo
          endif
          k=exc1(i)
2001     continue
        endif

C 2/10 here it starts the new version

C 2/10 first copy all the exlusion lists for the mutants and natives
C      involved in mutations on a new exclusion list.
!        if (1.eq.2) then
        mutcnt1 = 0

        if (totex.gt.0) then
         k=0
         l=0
         totexm=0
         do 2010 i=1,npt
          if (ipick(i).ne.0) then
             mutcnt1 = mutcnt1+1
             mutl(mutcnt1) = i
             exm1(mutcnt1) = exc1(i)-k+l
             do j=k+1,exc1(i)
             exm2(j-k+l)=exc2(j)
             totexm=totexm+1
             enddo
             l = exm1(mutcnt1)
          endif
          k=exc1(i)
2010     continue

       nmutant=mutcnt1

!       if (1.eq.1) then


C 2/10 Now, check on the previous particles if there are exclusions with the
C      native or the mutant and add it to the exclusion lists of these two
       k2=0
C      consider all natives and mutants
       DO 2100 i=1,nmutant
        k1=0
C      look at the particles that come before the mut/nat 
        DO 2200 j=1,mutl(i)-1
C      look at all the excluded prt for this prt
         DO 2300 k=k1+1,exc1(j)
C      is there our nat/mut there?
          if (exc2(k).eq.mutl(i)) then 
C      if yes, copy the list on a tmp vector, add the new excluded and copy all back
           DO l=k2+1,totexm
            exch(l)=exm2(l)
           ENDDO
           k2=k2+1
           DO kk=i,nmutant
           exm1(kk)=exm1(kk)+1
           ENDDO
           exm2(k2)=j
           totexm=totexm+1
           DO l=k2+1,totexm
            exm2(l)=exch(l-1)
           ENDDO
          endif
2300     CONTINUE
        k1=exc1(j)
2200    CONTINUE
        k2=exm1(i)
2100   CONTINUE

!       endif

C 3/10 remove the natives and mutants from all the exclusion list. In this way
C      they will not be count in the correction term in Ewald summation. Interaction
C      between mutants and natives will still not be counted since a controll that
C      removes them from nbond lists had been placed.

         k = 0

         DO 2020 i=1,npt

          if (ipick(i).ne.0) go to 2019 

2018      continue
          DO j=k+1,exc1(i)
           if (ipick(exc2(j)).ne.0) then          
           DO l=j+1,totex
            exch(l)=exc2(l)
           ENDDO
           DO l=i,npt
            exc1(l) = exc1(l) - 1
           ENDDO
           totex = totex-1
           DO l=j,totex
            exc2(l)=exch(l+1)
           ENDDO
           go to 2018
           endif
          ENDDO

2019     CONTINUE
         k = exc1(i)
2020     CONTINUE

C 3/10 now remove the mutants and natives from the original exclusion list

         k=0

         DO 2030 i=1,npt

          if (ipick(i).ne.0) then

           DO j=exc1(i)+1,totex
            exch(j)=exc2(j)
           ENDDO

           totex=totex-(exc1(i)-k)

           DO j=k+1,totex
            exc2(j) = exch(j+exc1(i)-k)
           ENDDO

           k1 = exc1(i) - k
       
           DO j=i,npt
            exc1(j) = exc1(j) - k1
           ENDDO

          endif

         k=exc1(i)
2030     CONTINUE

         endif

!         endif

c - special

        do i=1,maxtors
         rmlink(i)=0
        enddo

        if (totspe.gt.0) then
         do 2008 i=1,totspe
c - case 1-2
          if ((ipick(spec1(i)).ne.ipick(spec2(i))) .and.
     &         ipick(spec1(i)).ne.0 .and.
     &         ipick(spec2(i)).ne.0) then
           rmcnt=rmcnt+1
           rmlink(rmcnt)=i
          endif
2008     continue
        if (rmcnt.gt.0) then
         do i=1,rmcnt
          if (rmlink(i).lt.totspe) then
           do j=rmlink(i),totspe-1
            spec1(j)=spec1(j+1)
            spec2(j)=spec2(j+1)
           enddo
           do j=i,rmcnt
            rmlink(j)=rmlink(j)-1
           enddo
           totspe=totspe-1
          else
           totspe=totspe-1
          endif
         enddo
        end if
       end if



      return
      end
