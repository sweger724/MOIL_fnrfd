        subroutine calculate_kinetic_energy()
          implicit none

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/ENERGY.BLOCK'
          include 'COMMON/CONNECT.BLOCK'
          include 'COMMON/VELOC.BLOCK'

          integer i
          
           enkin = 0.d0

           do i = 1,npt
             enkin = enkin + ptms(i)*(velo(1,i)*velo(1,i)+
     1               velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
           end do

           enkin = 0.5d0*enkin

          return
        end
