      subroutine DEE_pairs1()
c
c Marks pair of rotamers as dead-end pairs. 
c In this routine we use the inequality for pairs rotamers i and j. 
c (Less restrictive Equation than the one used in DEE.f) 
c    (I.Lasters, M.DeMayer, J.Desmet, Protein Engineering 8, 815(1995))
c    e(g_i,h_j)-e(p_i,q_j)+
c       SUM_{k}{ min_{f'} {e[(g_i,h_j),f'_k]}-
c       SUM_{k}{ max_{f'} {e[(p_i,q_j),f'_k]}>0
c
c Added the fact that DEP imply DE-triples, and use this condition in 
c the first sum above (min), by analogy with the rotamer elimination 
c using DEP. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'

      
      character*10 name
      integer namel,level
c
      double precision EABmin(maxnpairslocal),EABmax(maxnpairslocal)
      double precision EABminlocal,EABmaxlocal
      double precision E1A1Bmin,E2A2Bmin,E1A1Bmax,E2A2Bmax
      integer iA,iB,indAB,indABlocal,ind1A1B,ind2A2B
      integer niterations,newdeadend,numberrotB
      integer pointEijA,pointEijB,pointEij1A,pointEij2A
      integer rotA,rotfirstA,rotlastA,rotAreal
      integer rot1A,rot2A,rot1Areal,rot2Areal
      integer rot1Atemp,rot2Atemp,rotAtemp,rot1Btemp,rot2Btemp,rotBtemp
      integer rotB,rotfirstB,rotlastB,rotBreal
      integer rot1B,rot2B,rot1Breal,rot2Breal
c
      name = 'DEE_pairs1'
      namel = 10
      level=1
c
c
      newdeadend=999
      niterations=0
c
      write (6,*) 
      write (6,*) 'DEE pairs 1:'
c
c
      do while (newdeadend.gt.0)
c
         newdeadend=0
         niterations=niterations+1
c
c........loop over all pair of enhanced positions (A & B)
         do iA=1,nposenh-1
c
            rotfirstA=poirotenhaux(iA-1)+1
            rotlastA=poirotenhaux(iA)
c
            do iB=iA+1,nposenh
c     
               rotfirstB=poirotenhaux(iB-1)+1
               rotlastB=poirotenhaux(iB)
               numberrotB=rotlastB-rotfirstB+1
c
c
c
c
c..............generate all possible pairs (AB) and calculate the 
c..............minimal and maximal energies with the other rotamers.
c
               do rotA=rotfirstA,rotlastA
c
                  rotAreal=rotaux(rotA)
                  rotAtemp=rotA-rotfirstA+1
                  pointEijA=pointEij(rotAreal-1)-poirotenh(iA)
c
                  do rotB=rotfirstB,rotlastB
c
                     rotBreal=rotaux(rotB)
                     rotBtemp=rotB-rotfirstB+1
                     pointEijB=pointEij(rotBreal-1)-poirotenh(iB)
c
                     indAB=pointEij(rotAreal-1)+rotBreal-
     &                    poirotenh(iA)
c
c....................if not a Dead-End-Pair:
                     if (.not.DEPij(indAB)) then
c
                        indABlocal=rotBtemp+(rotAtemp-1)*numberrotB
c
                        EABminlocal=Eiback(rotAreal)+
     &                       Eiback(rotBreal)+Eij(indAB)
                        EABmaxlocal=EABminlocal
c
                        call DEE_pairs1_innerloop(iA,iB,rotAreal,
     &                       rotBreal,pointEijA,pointEijB,EABminlocal,
     &                       EABmaxlocal)
c
                        EABmin(indABlocal)=EABminlocal
                        EABmax(indABlocal)=EABmaxlocal
c
                     end if
c
                  end do
c
               end do
c
c..............generate all possible combinations of pairs(1A1B & 2A2B):
c..............a)rot1A=rotfirstA,rotlastA-1; rot1B=rotfirstB,rotlastB
c..............  rot2A=rot1A+1,rotlastA;     rot2B=rotfirstB,rotlastB
c..............b)rot1A=rotfirstA,rotlastA;   rot1B=rotfirstB,rotlastB-1
c..............  rot2A=rot1A;                rot2B=rot1B+1,rotlastB
c
c
c..............a)
               do rot1A=rotfirstA,rotlastA-1
c
                  rot1Areal=rotaux(rot1A)
                  rot1Atemp=rot1A-rotfirstA+1
                  pointEij1A=pointEij(rot1Areal-1)-poirotenh(iA)
c
                  do rot1B=rotfirstB,rotlastB
c
                     rot1Breal=rotaux(rot1B)
                     rot1Btemp=rot1B-rotfirstB+1
c
                     ind1A1B=pointEij1A+rot1Breal
c
c....................if not a Dead-End-Pair:
                     if (.not.DEPij(ind1A1B)) then
c
                        indABlocal=rot1Btemp+(rot1Atemp-1)*numberrotB
                        E1A1Bmin=EABmin(indABlocal)
                        E1A1Bmax=EABmax(indABlocal)
c
c
                        do rot2A=rot1A+1,rotlastA
c
                           rot2Areal=rotaux(rot2A)
                           rot2Atemp=rot2A-rotfirstA+1
                           pointEij2A=pointEij(rot2Areal-1)-
     &                          poirotenh(iA)
c     
                           do rot2B=rotfirstB,rotlastB
c
                              rot2Breal=rotaux(rot2B)
                              rot2Btemp=rot2B-rotfirstB+1
c
                              ind2A2B=pointEij2A+rot2Breal
c
c.............................if not a Dead-End-Pair:
                              if (.not.DEPij(ind2A2B)) then
c
                                 indABlocal=rot2Btemp+
     &                                (rot2Atemp-1)*numberrotB
                                 E2A2Bmin=EABmin(indABlocal)
                                 E2A2Bmax=EABmax(indABlocal)
c
                                 if (E1A1Bmin.gt.E2A2Bmax) then
                                    DEPij(ind1A1B)=.true.
                                    newdeadend = newdeadend + 1
                                 else if (E2A2Bmin.gt.E1A1Bmax) then
                                    DEPij(ind2A2B)=.true.
                                    newdeadend = newdeadend + 1
                                 end if
c
c
                              end if
c
                           end do
c
                        end do
c
c
c
                     end if
c
                  end do
c
               end do
c
c
c
c
c..............b)
               do rot1A=rotfirstA,rotlastA
c
                  rot1Areal=rotaux(rot1A)
                  rot1Atemp=rot1A-rotfirstA+1
                  pointEij1A=pointEij(rot1Areal-1)-poirotenh(iA)
c
                  do rot1B=rotfirstB,rotlastB-1
c
                     rot1Breal=rotaux(rot1B)
                     rot1Btemp=rot1B-rotfirstB+1
c
                     ind1A1B=pointEij1A+rot1Breal
c     
c....................if not a Dead-End-Pair:
                     if (.not.DEPij(ind1A1B)) then
c     
                        indABlocal=rot1Btemp+(rot1Atemp-1)*numberrotB
                        E1A1Bmin=EABmin(indABlocal)
                        E1A1Bmax=EABmax(indABlocal)
c
                        rot2A=rot1A
c
                           rot2Areal=rotaux(rot2A)
                           rot2Atemp=rot2A-rotfirstA+1
                           pointEij2A=pointEij(rot2Areal-1)-
     &                          poirotenh(iA)
c
                           do rot2B=rot1B+1,rotlastB
c
                              rot2Breal=rotaux(rot2B)
                              rot2Btemp=rot2B-rotfirstB+1
c
                              ind2A2B=pointEij2A+rot2Breal
c
c.............................if not a Dead-End-Pair:
                              if (.not.DEPij(ind2A2B)) then
c
                                 indABlocal=rot2Btemp+
     &                                (rot2Atemp-1)*numberrotB
                                 E2A2Bmin=EABmin(indABlocal)
                                 E2A2Bmax=EABmax(indABlocal)
c
                                 if (E1A1Bmin.gt.E2A2Bmax) then
                                    DEPij(ind1A1B)=.true.
                                    newdeadend = newdeadend + 1
                                 else if (E2A2Bmin.gt.E1A1Bmax) then
                                    DEPij(ind2A2B)=.true.
                                    newdeadend = newdeadend + 1
                                 end if
c
c
                              end if
c
                           end do
c
                        continue
c
c
c
                     end if
c
                  end do
c
               end do
c
CDEB               write (6,*) 'iA,iB,newdeadend',iA,iB,newdeadend
CDEB               call flush(6)
c
c
            end do
c
         end do
c
         write (6,*) '  finished iteration',niterations
         write (6,*) '     new dead-end-pairs',newdeadend
         call flush(6)
c
c
      end do
c
c
      return
      end
c     
c
c-------------------------------------------------------------------
c
      subroutine DEE_pairs1_innerloop(iA,iB,rotAreal,rotBreal,pointEijA,
     &     pointEijB,EABminlocal,EABmaxlocal)
c
c Just the inner loop, splitted for clarity of the code. 
c
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'

      
      character*20 name
      integer namel,level
c
      double precision EABminlocal,EABmaxlocal,EijAB
      integer iA,iB,indA,indB,rotAreal,rotBreal
      integer pointEijA,pointEijB
c
      double precision emin,emax
      integer rotT,rotfirstT,rotlastT,rotTreal
      integer j
      integer pointEijT
c
c
      name = 'DEE_pairs1_innerloop'
      namel = 20
      level=1
c
c     
      do j=1,iA-1
c     
         emin=1.0d30
         emax=-1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
            pointEijT=pointEij(rotTreal-1)-poirotenh(j)
c
            indA=pointEijT+rotAreal
            indB=pointEijT+rotBreal
            EijAB=Eij(indA)+Eij(indB)
c
            emax=dmax1(emax,EijAB)
            if (.not.(DEPij(indA).or.DEPij(indB))) then
               emin=dmin1(emin,EijAB)
            end if
c
         end do
c     
         EABminlocal=EABminlocal+emin
         EABmaxlocal=EABmaxlocal+emax
c                     
      end do
c     
      do j=iA+1,iB-1
c
         emin=1.0d30
         emax=-1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
            pointEijT=pointEij(rotTreal-1)-poirotenh(j)
c
            indA=pointEijA+rotTreal
            indB=pointEijT+rotBreal
            EijAB=Eij(indA)+Eij(indB)
c
            emax=dmax1(emax,EijAB)
            if (.not.(DEPij(indA).or.DEPij(indB))) then
               emin=dmin1(emin,EijAB)
            end if
c
         end do
c     
         EABminlocal=EABminlocal+emin
         EABmaxlocal=EABmaxlocal+emax
c                                          
      end do
c     
c     
      do j=iB+1,nposenh
c
         emin=1.0d30
         emax=-1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
c
            indA=pointEijA+rotTreal
            indB=pointEijB+rotTreal
            EijAB=Eij(indA)+Eij(indB)
c
            emax=dmax1(emax,EijAB)
            if (.not.(DEPij(indA).or.DEPij(indB))) then
               emin=dmin1(emin,EijAB)
            end if
c
         end do
c     
         EABminlocal=EABminlocal+emin
         EABmaxlocal=EABmaxlocal+emax
c                     
      end do
c     
      return
      end







