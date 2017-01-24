        subroutine Getvec_Approx(vector,dv,rr)
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'
        
        double precision vector(3,*),dv(3,*),rr(3,*)
        double precision aaa
        integer i,l

        aaa=1.d-6

        do  i = 1,npt
          do l=1,3
             coor(l,i) =  rr(l,i)*massfac(i)+aaa*vector(l,i)
          end do
        end do

        call eforce()
        
        do  i = 1,npt
          do l=1,3
             vector(l,i) = 1.d0/aaa * ( dpot(l,i) - dv(l,i)/massfac(i) )
          end do
        end do

        return
        end
