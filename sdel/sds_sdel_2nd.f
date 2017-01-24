      subroutine sds_sdel(sener,dsall,r,d0,e0,nselec,pointr
     $     ,ipick,npt,gamma,fixend
     $     ,debug,ptms,time,clo)

      implicit none

c     a subroutine to calculate the action associated with SDEL.
C     It was developed from chmin routine so there are similarities
C     considering the treatment of the path as a polymer of system copies
c     The polymer energy (sener) is giving by
c     
c     S =     SUM  [ sqrt(2(E-U_i)) +sqrt(2(E-U_i+1)) ] dl_i,i+1
C             +   gamma SUM (d(i,i+1) - <d>)^2 
C             +   3/2 0.1(300kB)npt  SUM [ 1/(E-U_i) ]
c     monomer internal energies    nearest monomers harmonic attraction
c     
c     *** Note - in the current implementation the selection does not work!

      double precision sener,gamma
      integer nselec,npt,f
      double precision ptms(npt),time
      double precision dsall(3,*),r(3,*),d0(*),e0(*)
      integer pointr(*),ipick(*)
      logical fixend,debug
      
c     
c     common block for COORDinates and potential ENERGY derivatives
c     
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/SDEL.BLOCK'
      include 'COMMON/SGB.BLOCK'

C     Use this for calling matvec
      include 'COMMON/SCNDRV_SDEL.BLOCK'
      include 'COMMON/GETVEC.BLOCK'


c     pdq variables

C     p() carries our momenta of each sys.

      double precision dv1(3,MAXPT*LGRID)
      double precision p(LGRID),p_tmp(LGRID)
      double precision d0_tmp(LGRID),e0_tmp(LGRID), e1(LGRID)

C     Okay first optimization is to get rid of this
C     and just calculate ds() once..
      double precision dS(3,MAXPT*LGRID)
      
      integer i,j,k,l
      double precision tmp,alpha

      double precision clo,eij(lgrid)
      integer loop, Nloop

c     
c     Initialize our variables 

      alpha = 1.d-6
      sener = 0.d0

      do  i = 1,(pseg+2)*npt
        do l=1,3
          dsall(l,i) = 0.d0
        end do 
      end do

      do i = 1,pseg+1
         timea(i) = 0.0d0
      end do

      call Get_Dls(pseg,npt,r,d0,ipick)
C      call GetDave(dave,d0)
c     calculate contribution from nearest neighbours
      Nloop = pseg
      if (last) Nloop = pseg + 1
      tmp = 0.d0

C      write(6,*)"<dl>: ",dave
      do j = 1,Nloop
         tmp  = tmp + (d0(j) - dave)**2
C        write(6,*)"Dl: ",d0(j)
      end do
      sener = sener + gamma*tmp


C     include logarithm term:
      tmp  = 0.d0
      do j = 1,Nloop
         tmp =tmp -dlog(d0(j))
      enddo
      tmp = tmp * clo
      sener = sener +tmp

c     also save the exponents for derivatives calculations
c     
C      write(6,*)'dave',dave
      do 11 j = 1, pseg+1
         e1(j) = 2.d0*gamma*(1.d0 - dave/d0(j))
C        write(6,*)'d0(j)',d0(j)
C        write(6,*)'e1(j)',j,e1(j)
 11   continue
      
c     For the logarithmic term:
      do j=1, pseg+1
         eij(j) = -clo/(d0(j)**2)
      enddo
c     
c     calculate derivatives from nearest neighbour distance terms.
c     
      do j = 1,pseg+1
         k  = (j-1)*npt
         l  = k +npt
C     $DOIT IVDEP
         do i = 1,npt
            if (ipick(i).gt.0) then
              do loop=1,3
                tmp = e1(j)*(r(loop,k+i)-r(loop,l+i))
                dsall(loop,i+k)      = dsall(loop,i+k)   + tmp
                dsall(loop,i+l)      = dsall(loop,i+l)   - tmp
                tmp = eij(j)*(r(loop,k+i)-r(loop,l+i))
                dsall(loop,i+k)      = dsall(loop,i+k)   + tmp
                dsall(loop,i+l)      = dsall(loop,i+l)   - tmp
              end do
            end if
         end do
      end do
        
C       Calculate the energies and forces
                
        call Get_EnergyDerivatives(pseg,npt,r,dv1,e0)
c       write(6,*)"energy computed"

C     Add penalty for potential energy approaching the Total energy (ENERGYC)
C  this is due algorithm stability, because if potential energy exceeds
C total energy the momenta would be imaginary and algorithm would crash
        tmp = 0.1d0*kboltzmann*300*1.5d0*npt
        tmp = FORCE_ENERGY * tmp
        do j = 2,pseg+1
          sener = sener + tmp/(ENERGYC-e0(j))
C     store also multiplicative factor for derivative calculation
          e1(j) = tmp/(ENERGYC - e0(j))**2
        enddo

C  Add the derivatives with respect to the previous penalty
C        d.../dx(j) = dv(j)/(ENERGYC-e0(i))
        do j = 2,pseg+1
          k  = (j-1)*npt
          do i = 1,npt
            if (ipick(i).gt.0) then
              do loop=1,3
                dsall(loop,i+k) = dsall(loop,i+k) + e1(j)*dv1(loop,i+k)
              end do
            end if
         end do
        enddo

C     --------------------------------------------------------------
C     Now calculate our Action
               pdqU = 0.0d0
               pdqS = 0.0d0
               pdqA = 0.0d0
               pdqT = 0.0d0

C     First calculate our momenta for each molecule 
C     and save it 
          call Get_Momenta(pseg,ENERGYC,e0,p)
          call Communicate_Momenta(p)
          
C       Calculate 2* our action to pdqS
              do i = 1,Nloop
                  pdqS = pdqS + ( p(i)+p(i+1) ) * d0(i)
              enddo

C     Now calculate dS
              call Get_dS(pseg,npt,r,dv1,p,d0,pdqA,dS)                  
              call Communicate_Positions(npt,dS,.false.)
              call Communicate_Positions(npt,dv1,.false.)
            
              sener =  sener + pdqA
                
              call T_deriv(dsall,r,d0,p,dv1,dS,pseg,npt)

789        return
            end

            
C*********************************************************************
        subroutine Get_EnergyDerivatives(pseg,npt,rr,ddv,e0)
C*********************************************************************
        implicit none

        double precision rr(3,*),ddv(3,*),e0(*)
        integer pseg, npt
        
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SDEL.BLOCK'
        include 'COMMON/ENERGY.BLOCK'

        integer j,i,l,k
C     First calculate energies and forces
            do  5 j = 2,pseg+1

C     Lets convert this back to non-massweighted before our energy calls
C     and nbondm etc ..
               k = (j-1)*npt
               do  i = 1,npt
                 do l=1,3
                   coor(l,i) =  rr(l,k+i)*massfac(i)
                 end do
               end do
c
c     energy call for the internal energy of monomer j
                call eforce()
C               call wener(stdo)

c     e0 stores the internal energy of the monomers
c     e_total is obtained from ENER common block

                e0(j) = e_total
C              write(*,*) 'Energy: ',j, e0(j)

C     Since we are taking the derivative w/ respect to the
C     mass weighted coordinates we take a 1/sqrt(ptms)
                  do  i = 1,npt
                    do l =1,3
                      ddv(l,k+i) = dpot(l,i)*massfac(i)
                    end do
                  end do

c     *** End of internal energy calculations
 5          continue
        return
        end


C*********************************************************************
        subroutine Get_dS(pseg,npt,rr,ddv,pp,dd0,dS2,dS)
C*********************************************************************
C     Now calculate dS
        implicit none

        integer pseg,npt
        double precision rr(3,*),ddv(3,*),pp(*),dd0(*),dS2,dS(3,*)

        include 'COMMON/LENGTH.BLOCK'
        
        integer j,k,m,loop
        double precision dU_dqjm,t1,t2

        dS2=0.d0

        do j = 2,pseg+1
           k = (j-1) * npt
           do m = 1,npt
             do loop = 1,3
C     Here we calculate dS only
               dU_dqjm= - ddv(loop,k+m)
               t1 = rr(loop,k+m) - rr(loop,k+npt+m)
               t2 = rr(loop,k+m) - rr(loop,k-npt+m)

               dS(loop,k+m) = dU_dqjm /pp(j) * (dd0(j-1) + dd0(j))
     $                      + ( pp(j) + pp(j+1) ) * t1/dd0(j)
     $                      + ( pp(j) + pp(j-1) ) * t2/dd0(j-1)

               dS2 = dS2 + dS(loop,k+m)**2
c     write (6,*) 'debugging pdqA ',term0,p(i),t1,t2,pdqA
c     write (6,*) i,p(i),term0
             end do
           end do
        end do
        return
        end

C*********************************************************************
        subroutine Get_Momenta(pseg,maxE,ee,pp)
C*********************************************************************
        implicit none
        
        integer pseg,i
        double precision maxE,ee(*),pp(*)

          do i = 2,pseg+1
            pp(i) = 2 * (maxE-ee(i))
            if (pp(i).le.0) then
               write (6,*) 'Our momentum is < 0',i,ee(i)
               stop
            end if
            pp(i) = sqrt(pp(i))
          end do
          return
         end

C*********************************************************************
        subroutine Get_Dls(pseg,npt,rr,dls,ipick)
C*********************************************************************
c     calculate now distances between all i,i+1 pairs
c     and the corresponding energy terms

      implicit none

      integer pseg,npt,ipick(*)
      double precision rr(3,*),dls(*)

      integer i,j,k,loop
      double precision tmpx

      do 56 j = 1,pseg+1
         dls(j)=0.d0
         k = (j-1) * npt
         do 55 i = 1,npt
C           write(6,*)"ipick(i)",i,ipick(i)
            if (ipick(i).gt.0) then
              do loop =1,3
                 tmpx = rr(loop,i+k+npt) - rr(loop,i+k)
C                write(6,*)"rr",rr(loop,i+k)
                 dls(j) = dls(j) + tmpx**2
              end do
            end if
 55      continue

         dls(j) = dsqrt(dls(j))
c         write(6,*)"dls(j)",j,dls(j)
 56   continue

      return
      end
