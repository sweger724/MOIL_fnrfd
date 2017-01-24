C     Routine to mass weight our coordinate vector (r()) 

      subroutine mass_weight_coordinates(NN,rr)
      implicit none
      double precision rr(3,*)
      integer i,iloc,j,NN
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/MASSFAC.BLOCK'

      do i = 1,NN
         iloc = (i-1)*npt
         do j = 1,npt
            rr(1,iloc + j) = rr(1,iloc+j) /massfac(j)
            rr(2,iloc + j) = rr(2,iloc+j) /massfac(j)
            rr(3,iloc + j) = rr(3,iloc+j) /massfac(j)
         end do
      end do
      return 
      end
