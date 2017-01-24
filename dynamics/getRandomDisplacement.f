        subroutine getRandomDisplacement(dr,dv)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/LD.BLOCK'

        double precision dr(3,maxpt),dv(3,maxpt)
        double precision rr(6*maxpt)
        integer i,j


        call gauss2(6*npt,rr)
        
        do i = 1,npt
          do j = 1,3
            dr(j,i) = sigmaR * rr((i-1)*6+2*j-1) / dsqrt(ptms(i))
            dv(j,i) = sigmaV / dsqrt(ptms(i)) * 
     &      (crv*rr((i-1)*6+2*j-1) + dsqrt(1.d0-crv**2)*rr((i-1)*6+2*j))
C            write(6,*)"dr,dv:",dr(j,i),dv(j,i)
          end do
        end do
        

        return

        end
