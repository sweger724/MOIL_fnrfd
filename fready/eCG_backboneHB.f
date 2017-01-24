        Subroutine eCG_backboneHB()
C      HB backbone potential according to
C      J. Phys. Chem. B 2004, 108, 9421 - 9438
C      Parametrization of Backbone-Electrostatic and Multibody Contributions to the
C      UNRES Force Field for Protein-Structure Prediction from Ab Initio Energy Surfaces
C      of Model Systems;  Liwo, Oldziej, Czaolewski, Kozlowska, Scheraga

C  04/17/2008
C  Significant changes introduced. Only best HB partners
C  are considered... not all partners
C
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision rrj(3),rrl(3),Ej,El
        double precision pi(3),di(3)
        double precision pj(3),dj(3), pl(3),dl(3)
        double precision si,sj,sij,sl,sil,Rj(3),Rl(3),sRl,sRj
        double precision cosAj,cosBj,cosCj,dEjdr,dEjdA,dEjdB,dEjdC
        double precision cosAl,cosBl,cosCl,dEldr,dEldA,dEldB,dEldC
        double precision fA, dfA, factor, d12di(3),d12dj(3),d12dl(3)
        double precision cos12
        integer i,j,l,T,jpt,ipt,j1pt,i1pt,lpt,l1pt
        integer m
        
        e_el = 0.d0

C   First precalculate energy of each posible
C   hydrogen bonding peptide pair
        do i=1,totmon
        
          if (HBbond1(i) .ne. 0) then

                ipt = poipt(i-1)+1
                i1pt = poipt(i)+1

                do m = 1,3
                  pi(m) = 0.5d0 * (coor(m,ipt) + coor(m,i1pt))
                  di(m) =coor(m,i1pt) - coor(m,ipt) 
                end do

                j=HBbond1(i)

                  jpt = poipt(j-1)+1
                  j1pt = poipt(j)+1

                  do m =1,3
                    pj(m) = 0.5d0 * (coor(m,jpt) + coor(m,j1pt))
                    dj(m) = coor(m,j1pt) - coor(m,jpt)
                    rrj(m) = pj(m) - pi(m)
                  end do

                  si = dsqrt(di(1)**2 + di(2)**2 + di(3)**2)

                  sj = dsqrt(dj(1)**2 + dj(2)**2 + dj(3)**2)
                  sij = dsqrt(rrj(1)**2 + rrj(2)**2 + rrj(3)**2)

                  sl = dsqrt(dl(1)**2 + dl(2)**2 + dl(3)**2)
                  sil = dsqrt(rrl(1)**2 + rrl(2)**2 + rrl(3)**2)

C  calculate cos Alpha, Beta, Gamma angles for partner j

               cosAj = (di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3))/(si*sj)
               cosBj = (di(1)*rrj(1) + di(2)*rrj(2) + di(3)*rrj(3))
     &               /(si*sij)
               cosCj = (dj(1)*rrj(1) + dj(2)*rrj(2) + dj(3)*rrj(3))
     &               /(sj*sij)
               T = HBtype1(i)

                call eCG_backboneNB(T,sij,cosAj,cosBj,cosCj,Ej,dEjdr,
     &                              dEjdA,dEjdB,dEjdC,.TRUE.)


             if (HBbond2(i) .ne. 0) then 
                
                  l = HBbond2(i)
                  lpt = poipt(l-1)+1
                  l1pt = poipt(l)+1

                  do m =1,3
                    pl(m) = 0.5d0 * (coor(m,lpt) + coor(m,l1pt))
                    dl(m) = coor(m,l1pt) - coor(m,lpt)
                    rrl(m) = pl(m) - pi(m)
                  end do
                
                  sl = dsqrt(dl(1)**2 + dl(2)**2 + dl(3)**2)
                  sil = dsqrt(rrl(1)**2 + rrl(2)**2 + rrl(3)**2)

C  calculate cos Alpha, Beta, Gamma angles for partner l


               cosAl = (di(1)*dl(1) + di(2)*dl(2) + di(3)*dl(3))/(si*sl)
               cosBl = (di(1)*rrl(1) + di(2)*rrl(2) + di(3)*rrl(3))
     &               /(si*sil)
               cosCl = (dl(1)*rrl(1) + dl(2)*rrl(2) + dl(3)*rrl(3))
     &               /(sl*sil)
               T = HBtype2(i)

               call eCG_backboneNB(T,sil,cosAl,cosBl,cosCl,El,dEldr,
     &                              dEldA,dEldB,dEldC,.TRUE.)

C  calculate the cosine of angle between l i j (cos12)
C and its derivative contribution to the potential
                  do m=1,3
                    Rj(m) = pj(m) - pi(m)
                    Rl(m) = - ( pl(m) - pi(m) )
                  end do

                  sRl = dsqrt(Rl(1)**2 + Rl(2)**2 + Rl(3)**2)
                  sRj = dsqrt(Rj(1)**2 + Rj(2)**2 + Rj(3)**2)

             cos12 = (Rl(1)*Rj(1) + Rl(2)*Rj(2) + Rl(3)*Rj(3))/(sRl*sRj)

                e_el=e_el+(Ej+El)*fA(cos12)
                factor = 0.5d0 * dfA(cos12) * (Ej+El)
                
                  do m =1,3
                    d12di(m) = (Rj(m)-Rl(m)) / (sRl*sRj)
     &                      - cos12 * (Rl(m)/sRl**2 - Rj(m)/sRj**2)
                    d12dj(m) = Rl(m)/(sRl*sRj) - cos12 * Rj(m)/sRj**2
                    d12dl(m) =-Rj(m)/(sRl*sRj) + cos12 * Rl(m)/sRl**2

                    d12di(m) = d12di(m) * factor
                    d12dj(m) = d12dj(m) * factor
                    d12dl(m) = d12dl(m) * factor

                    dpot(m,ipt)  = dpot(m,ipt)  + d12di(m)
                    dpot(m,i1pt) = dpot(m,i1pt) + d12di(m)
                    dpot(m,jpt)  = dpot(m,jpt)  + d12dj(m)
                    dpot(m,j1pt) = dpot(m,j1pt) + d12dj(m)
                    dpot(m,lpt)  = dpot(m,lpt)  + d12dl(m)
                    dpot(m,l1pt) = dpot(m,l1pt) + d12dl(m)
                  end do

          call BasicHBDerivative(dEjdr,dEjdA,dEjdB,dEjdC,rrj,di,dj,si,sj
     &        ,i,j,sij,cosAj,cosBj,cosCj,fA(cos12),ipt,i1pt,jpt,j1pt,Ej)
          call BasicHBDerivative(dEldr,dEldA,dEldB,dEldC,rrl,di,dl,si,sl
     &        ,i,l,sil,cosAl,cosBl,cosCl,fA(cos12),ipt,i1pt,lpt,l1pt,El)

          else   
       
          e_el=e_el+ Ej          
          
          call BasicHBDerivative(dEjdr,dEjdA,dEjdB,dEjdC,rrj,di,dj,si,sj
     &        ,i,j,sij,cosAj,cosBj,cosCj,1.d0,ipt,i1pt,jpt,j1pt,Ej)
        
          endif

               
          end if ! if (HBbond1(i) .ne. 0)
               
        end do   ! i= 1, totmon


        return
        end

C*************************************
C    Do basic HB derivative
C*************************************
      subroutine BasicHBDerivative(dEdr,dEdA,dEdB,dEdC,r,di,dj,si,sj,
     &              i,j,sij,cosA,cosB,cosC,cos12,ipt,i1pt,jpt,j1pt,E)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'      
        include 'COMMON/FREADY.BLOCK'

        double precision E,dEdr,dEdA,dEdB,dEdC
        double precision r(3),di(3),dj(3),si,sj,sij,cosA,cosB,cosC,cos12
        integer ipt, jpt,i1pt,j1pt

C   LOCAL

        double precision drdi(3),dAdi(3),dAdj(3)
        double precision dBdi(3),dBdi1(3),dBdj(3)
        double precision dCdj(3),dCdj1(3),dCdi(3)
        integer m,i,j

                if(E .gt. E_CG_max) then
                  write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~'
                  write(6,*)'i and j are ',i,j
                  write(6,*)'r,A,B,C: ',sij,cosA,cosB,cosC
                  write(6,'(a7,2(x,f8.3))')'eHB is ',
     &                  cos12,E
                  write(6,*)"dE is ",dEdr,dEdA,dEdB,dEdC
                  write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~~'
                 end if 
               
C                write(6,*)"HydrogenBond ",E,sij
               
                dEdr = dEdr * cos12
                dEdA = dEdA * cos12
                dEdB = dEdB * cos12
                dEdC = dEdC * cos12

              do m =1,3
                drdi(m) = -0.5d0 * r(m)/sij
C  drdj = -drdi
C  drdi1 = dzdi
C  drdj1 = dzdj

                dAdi(m) = -dj(m)/(si*sj) + cosA * di(m)/si**2
C dAdi1 = - dAdi

                dAdj(m) = -di(m)/(si*sj) + cosA * dj(m)/sj**2
C dAdj1 = -dAdj

                dBdi(m) = (-r(m)-0.5d0*di(m))/(si*sij) + 
     &                  cosB*(di(m)/si**2 + 0.5d0*r(m)/sij**2)

                dBdi1(m) = (r(m)-0.5d0*di(m))/(si*sij) - 
     &                  cosB*(di(m)/si**2 - 0.5d0*r(m)/sij**2)

                dBdj(m) = 0.5d0 *(di(m)/(si*sij) - cosB * r(m)/sij**2)
C  dBdj1(m) = dBdj(m)

                dCdj(m) = (-r(m)+0.5d0*dj(m))/(sj*sij) 
     &                + cosC*(dj(m)/sj**2 - 0.5d0*r(m)/sij**2)

                dCdj1(m) = (r(m)+0.5d0*dj(m))/(sj*sij) 
     &                 - cosC*(dj(m)/sj**2 + 0.5d0*r(m)/sij**2)

                dCdi(m) = -0.5d0 * (dj(m)/(sj*sij) - cosC * r(m)/sij**2)
C   dCdi1(m) = dCdi(m)
                
                dpot(m,ipt) = dpot(m,ipt) + dEdr*drdi(m) + dEdA*dAdi(m)
     &                                    + dEdB*dBdi(m) + dEdC*dCdi(m)

                dpot(m,i1pt)=dpot(m,i1pt) + dEdr*drdi(m) - dEdA*dAdi(m)
     &                                   + dEdB*dBdi1(m) + dEdC*dCdi(m)

                dpot(m,jpt) = dpot(m,jpt) - dEdr*drdi(m) + dEdA*dAdj(m)
     &                                    + dEdB*dBdj(m) + dEdC*dCdj(m)

                dpot(m,j1pt)=dpot(m,j1pt) - dEdr*drdi(m) - dEdA*dAdj(m)
     &                                    + dEdB*dBdj(m) + dEdC*dCdj1(m)
             end do  ! m

        return 
        end

C*************************************        
        function fA(c)
          implicit none
          double precision fA,c

          if (c .gt. 0.3 .and. c.lt.0.9) then
            fA = (c - 0.3)/0.6
          else if (c .ge. 0.9) then
            fA = 1.d0
          else
            fA = 0.d0
          end if

          return 
        end

C************************************
        function dfA(c)
          implicit none
          double precision dfA,c
        
          if (c .gt. 0.3 .and. c.lt.0.9) then
            dfA = 1.d0/0.6
          else
            dfA = 0.d0
          end if 
          
          return
        end



