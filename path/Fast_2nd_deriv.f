        subroutine Fast_2nd_deriv(vector)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'

        double precision vector(3,maxpt)
        double precision grad1(3,maxpt),a
        integer i

        a = 1.d-9 

        do i=1,npt
c       generate temp coordinates
           coor(1,i) = coor(1,i) + a*vector(1,i)
           coor(2,i) = coor(2,i) + a*vector(2,i)
           coor(3,i) = coor(3,i) + a*vector(3,i)
        end do

C       store initial gradient
        do i=1,npt
           grad1(1,i) = dpot(1,i)
           grad1(2,i) = dpot(2,i)
           grad1(3,i) = dpot(3,i)
        end do 

        call eforce()
        
        do i=1,npt
c       recover original coordinates
           coor(1,i) = coor(1,i) - a*vector(1,i)
           coor(2,i) = coor(2,i) - a*vector(2,i)
           coor(3,i) = coor(3,i) - a*vector(3,i)
        end do

        do i=1,npt
c       compute H_U * vector
           vector(1,i) = (dpot(1,i)-grad1(1,i))/a
           vector(2,i) = (dpot(2,i)-grad1(2,i))/a
           vector(3,i) = (dpot(3,i)-grad1(3,i))/a
           if(abs(dpot(1,i)-grad1(1,i)).gt.0.0001)
     &       write(6,'(a,f8.3,1x,f8.4,1x,e10.3,1x,e10.3)')
     &       "Warning: ",dpot(1,i),grad1(1,i),a,vector(1,i)
        end do

        return
        end
