        subroutine sds(sener,d0,ipick,npt,clo)
        implicit none
c
c  a subroutine to calculate the energy of a polymer of system copies
c  following the overdamped Onsager Machlup action.
c
c  S_l =        SUM i=1..n sqrt( Hamilt + (dU/dx_i)^2 ) *d0_i  


c
c common block for COORDinates and potential ENERGY derivatives
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/PATH.BLOCK'

c
c Local
c e0  -  vector of the individual energies of the monomers
c
        integer npt,i,j,k,l,ipick(*)
        
        double precision dx1,dx2
        double precision dStmp(lgrid), dStmp2(lgrid)
        double precision dStmp3(lgrid,3,maxpt),tmp(3,maxpt)
        double precision dU2, sener,d0(*),clo
c
c initialiaztion
c
        call sds_constraint(sener,d0,ipick,npt,clo)

C        return

C       write(6,*)"Constrain, gC: ",sener, gC
        do 5 i=2,pseg+1
          k = (i-1)*npt
          
          do 3 j=1,npt
            do 3 l=1,3
                coor(l,j) = r(l,j+k)
3         continue
        
          call eforce()
C          call wener(6)
          e0(i)=e_total

          do j=1,npt
            do l=1,3
              tmp(l,j)   = dpot(l,j)
            end do
          end do

          call Fast_2nd_deriv(tmp(1,1))
          do j=1,npt
            do l=1,3
              dStmp3(i,l,j)   = tmp(l,j)
            end do
          end do
          
          dU2 = 0.d0
          do j=1,npt
            do l=1,3
              dU2 = dU2 + dpot(l,j)**2
            end do
          end do

          dStmp(i) = dsqrt(Hamilt + dU2)
          sener = sener + dStmp(i)*d0(i) 

          if( d0(i)**2 .gt. 1.d-12) then
            dStmp2(i) = dStmp(i)/d0(i)
          else
            dStmp2(i) = 0.0
            write(6,*)"Warning, d0(i) too small",i,d0(i)
          endif
          
5       continue

C        return

          call Communicate_dStmp(dStmp)
          dStmp2(1) = dStmp(1)/d0(1)

        do 15 i=2,pseg+1
            k = (i-1)*npt
            do j=1,npt
              do l=1,3
                dx1 = r(l,k+j+npt) - r(l,k+j)
                dx2 = r(l,k+j)     - r(l,k+j-npt)
                dsall(l,k+j)=dsall(l,k+j)+dStmp2(i-1)*dx2-dStmp2(i)*dx1
              end do
            end do

            do j=1,npt
              do l=1,3
                dsall(l,k+j) = dsall(l,k+j)+dStmp3(i,l,j)/dStmp(i)*d0(i)
              end do
            end do

15       continue

        return
        end
