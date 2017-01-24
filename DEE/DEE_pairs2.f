      subroutine DEE_pairs2()
c
c Marks pair of rotamers as dead-end pairs. 
c In this routine we use the inequality for pairs rotamers i and j. 
c (Equation analogous to the one used in DEE.f) 
c    (I.Lasters, M.DeMayer, J.Desmet, Protein Engineering 8, 815(1995))
c    e(g_i,h_j)-e(p_i,q_j)+
c       SUM_{k}{ min_{f'} {e[(g_i,h_j),f'_k]-e[(p_i,q_j),f'_k]}>0
c
c It is the first check, therefore it does not use the iterative search
c to save CPU time. 
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
      double precision E1A1B,E2A2B
      double precision dtotal1,dtotal2
      integer iA,iB,newdeadend,ind1A1B,ind2A2B,niterations
      integer rot1A,rot2A,rotfirstA,rotlastA,rot1Areal,rot2Areal
      integer rot1B,rot2B,rotfirstB,rotlastB,rot1Breal,rot2Breal
      integer pointEij1A,pointEij1B,pointEij2A,pointEij2B
c
      name = 'DEE_pairs2'
      namel = 10
      level=1
c
c
c
      write (6,*) 
      write (6,*) 'DEE pairs 2:'
c
c
c
      newdeadend=0
c
c.....loop over all pair of enhanced positions (A & B)
      do iA=1,nposenh-1
c
         rotfirstA=poirotenhaux(iA-1)+1
         rotlastA=poirotenhaux(iA)
c
         do iB=iA+1,nposenh
c     
            rotfirstB=poirotenhaux(iB-1)+1
            rotlastB=poirotenhaux(iB)
c
c...........generate all possible combinations of pairs(1A1B & 2A2B):
c...........a)rot1A=rotfirstA,rotlastA-1; rot1B=rotfirstB,rotlastB
c...........  rot2A=rot1A+1,rotlastA;     rot2B=rotfirstB,rotlastB
c...........b)rot1A=rotfirstA,rotlastA;   rot1B=rotfirstB,rotlastB-1
c...........  rot2A=rot1A;                rot2B=rot1B+1,rotlastB
c
c
c...........a)
            do rot1A=rotfirstA,rotlastA-1
c
               rot1Areal=rotaux(rot1A)
               pointEij1A=pointEij(rot1Areal-1)-poirotenh(iA)
c
               do rot1B=rotfirstB,rotlastB
c
                  rot1Breal=rotaux(rot1B)
                  pointEij1B=pointEij(rot1Breal-1)-poirotenh(iB)
c
                  ind1A1B=pointEij1A+rot1Breal
c
c.................if not a Dead-End-Pair:
                  if (.not.DEPij(ind1A1B)) then
c
                     E1A1B=Eiback(rot1Areal)+Eiback(rot1Breal)+
     &                    Eij(ind1A1B)
c
c
c
                     do rot2A=rot1A+1,rotlastA
c
                        rot2Areal=rotaux(rot2A)
                        pointEij2A=pointEij(rot2Areal-1)-
     &                       poirotenh(iA)
c
                        do rot2B=rotfirstB,rotlastB
c
                           rot2Breal=rotaux(rot2B)
                           pointEij2B=pointEij(rot2Breal-1)-
     &                          poirotenh(iB)
c
                           ind2A2B=pointEij2A+rot2Breal
c
c..........................if not a Dead-End-Pair:
                           if (.not.DEPij(ind2A2B)) then
c
                              E2A2B=Eiback(rot2Areal)+
     &                             Eiback(rot2Breal)+
     &                             Eij(ind2A2B)
                              dtotal1=E1A1B-E2A2B
                              dtotal2=-dtotal1
c
                              call DEE_pairs2_innerloop(iA,iB,
     &                             rot1Areal,rot1Breal,rot2Areal,
     &                             rot2Breal,pointEij1A,
     &                             pointEij1B,pointEij2A,
     &                             pointEij2B,dtotal1,dtotal2)
c
                              if (dtotal1.gt.0) then
                                 DEPij(ind1A1B)=.true.
                                 newdeadend = newdeadend + 1
                              else if (dtotal2.gt.0) then
                                 DEPij(ind2A2B)=.true.
                                 newdeadend = newdeadend + 1
                              end if
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
c...........b)
            do rot1A=rotfirstA,rotlastA
c
               rot1Areal=rotaux(rot1A)
               pointEij1A=pointEij(rot1Areal-1)-poirotenh(iA)
c     
               do rot1B=rotfirstB,rotlastB-1
c
                  rot1Breal=rotaux(rot1B)
                  pointEij1B=pointEij(rot1Breal-1)-poirotenh(iB)
c
                  ind1A1B=pointEij1A+rot1Breal
c
c.................if not a Dead-End-Pair:
                  if (.not.DEPij(ind1A1B)) then
c
                     E1A1B=Eiback(rot1Areal)+Eiback(rot1Breal)+
     &                    Eij(ind1A1B)
c
c
c
                     rot2A=rot1A
c
                        rot2Areal=rotaux(rot2A)
                        pointEij2A=pointEij(rot2Areal-1)-
     &                       poirotenh(iA)
c
                        do rot2B=rot1B+1,rotlastB
c     
                           rot2Breal=rotaux(rot2B)
                           pointEij2B=pointEij(rot2Breal-1)-
     &                          poirotenh(iB)
c
                           ind2A2B=pointEij2A+rot2Breal
c
c..........................if not a Dead-End-Pair:
                           if (.not.DEPij(ind2A2B)) then
c
                              E2A2B=Eiback(rot2Areal)+
     &                             Eiback(rot2Breal)+
     &                             Eij(ind2A2B)
                              dtotal1=E1A1B-E2A2B
                              dtotal2=-dtotal1
c
                              call DEE_pairs2_innerloop(iA,iB,
     &                             rot1Areal,rot1Breal,rot2Areal,
     &                             rot2Breal,pointEij1A,
     &                             pointEij1B,pointEij2A,
     &                             pointEij2B,dtotal1,dtotal2)
c
                              if (dtotal1.gt.0) then
                                 DEPij(ind1A1B)=.true.
                                 newdeadend = newdeadend + 1
                              else if (dtotal2.gt.0) then
                                 DEPij(ind2A2B)=.true.
                                 newdeadend = newdeadend + 1
                              end if
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
c
CDEB            write (6,*) 'iA,iB,newdeadend',iA,iB,newdeadend
CDEB            call flush(6)
c
c
         end do
c
      end do
c
      niterations = 0 
      write (6,*) '  finished iteration',niterations
      write (6,*) '     new dead-end-pairs',newdeadend
      call flush(6)
c
c
c
      return
      end
c     
c
c-------------------------------------------------------------------
c
      subroutine DEE_pairs2_innerloop(iA,iB,rot1Areal,rot1Breal,
     &     rot2Areal,rot2Breal,pointEij1A,pointEij1B,pointEij2A,
     &     pointEij2B,dtotal1,dtotal2)
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
      double precision dtotal1,dtotal2
      integer iA,iB,rot1Areal,rot1Breal,rot2Areal,rot2Breal
      integer pointEij1A,pointEij1B,pointEij2A,pointEij2B
c
      double precision ediffmin1,ediffmin2,ediff1,ediff2
      double precision Eij1A1B,Eij2A2B
      integer rotT,rotfirstT,rotlastT,rotTreal
      integer j,ind1A,ind1B,ind2A,ind2B
      integer pointEijT
c
c
      name = 'DEE_pairs2_innerloop'
      namel = 20
      level=1
c
c
c
      do j=1,iA-1
c     
         ediffmin1=1.0d30
         ediffmin2=1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
            pointEijT=pointEij(rotTreal-1)-poirotenh(j)
c
            ind1A=pointEijT+rot1Areal
            ind1B=pointEijT+rot1Breal
            ind2A=pointEijT+rot2Areal
            ind2B=pointEijT+rot2Breal
            Eij1A1B=Eij(ind1A)+Eij(ind1B)
            Eij2A2B=Eij(ind2A)+Eij(ind2B)
c
            ediff1=Eij1A1B-Eij2A2B
            ediff2=-ediff1
            ediffmin1=dmin1(ediffmin1,ediff1)
            ediffmin2=dmin1(ediffmin2,ediff2)
c
         end do
c     
         dtotal1=dtotal1+ediffmin1
         dtotal2=dtotal2+ediffmin2
c                     
      end do
c     
      do j=iA+1,iB-1
c
         ediffmin1=1.0d30
         ediffmin2=1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
            pointEijT=pointEij(rotTreal-1)-poirotenh(j)
c
            ind1A=pointEij1A+rotTreal
            ind1B=pointEijT +rot1Breal
            ind2A=pointEij2A+rotTreal
            ind2B=pointEijT +rot2Breal
            Eij1A1B=Eij(ind1A)+Eij(ind1B)
            Eij2A2B=Eij(ind2A)+Eij(ind2B)
c
            ediff1=Eij1A1B-Eij2A2B
            ediff2=-ediff1
            ediffmin1=dmin1(ediffmin1,ediff1)
            ediffmin2=dmin1(ediffmin2,ediff2)
c
         end do
c     
         dtotal1=dtotal1+ediffmin1
         dtotal2=dtotal2+ediffmin2
c                     
      end do
c     
c     
      do j=iB+1,nposenh
c
         ediffmin1=1.0d30
         ediffmin2=1.0d30
c     
         rotfirstT=poirotenhaux(j-1)+1
         rotlastT=poirotenhaux(j)
c
         do rotT=rotfirstT,rotlastT
            rotTreal=rotaux(rotT)
c
            ind1A=pointEij1A+rotTreal
            ind1B=pointEij1B+rotTreal
            ind2A=pointEij2A+rotTreal
            ind2B=pointEij2B+rotTreal
            Eij1A1B=Eij(ind1A)+Eij(ind1B)
            Eij2A2B=Eij(ind2A)+Eij(ind2B)
c
            ediff1=Eij1A1B-Eij2A2B
            ediff2=-ediff1
            ediffmin1=dmin1(ediffmin1,ediff1)
            ediffmin2=dmin1(ediffmin2,ediff2)
c
         end do
c     
         dtotal1=dtotal1+ediffmin1
         dtotal2=dtotal2+ediffmin2
c                     
      end do
c     
c     
      return
      end







