        subroutine deriv_test(nptstart,nptend)

        implicit none
        integer nptstart,nptend

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
c
c local
        integer i,j
        double precision dr,e_plus,e_minus
        double precision grad1(3,maxpt)

        dr = 0.001d0

        call eforce()
        do i = nptstart,nptend
         do j = 1,3
          grad1(j,i) = dpot(j,i)
         end do
        end do

            do i=nptstart,nptend
             do j = 1,3
c       generate temp coordiantes
               coor(j,i) = coor(j,i) + dr
               call eforce()
               e_plus = e_total
               coor(j,i) = coor(j,i) -2.d0*dr
               call eforce()
               e_minus = e_total
               grad1(j,i) = grad1(j,i) - (e_plus-e_minus)/(2.d0*dr)
               coor(j,i) = coor(j,i) + dr
            end do
           end do

            call eforce()

            do i=nptstart,nptend
             do j = 1,3
              write(*,*)' i j dpot DDpot ',i,j,dpot(j,i),grad1(j,i)
             end do
            end do
             
        return
        end
