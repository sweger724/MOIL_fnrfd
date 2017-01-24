       subroutine initPath()
       implicit none
c
c calculate path using SDEL algorithm
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/SDP.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'

        integer namel,i,l
        character*9 name
        double precision dU2

        namel = 9
        name = 'init_sdel'
        MASSWEIGHT = 0

c Check whether positive hamiltonian
        if (Hamilt.lt. 1.d-10) then
         call alert(name,namel,'Hamiltonian must be positive',16,1)
        else
          write(stdo,*)"Constant hamiltonian: ", Hamilt
        endif

        if (first) then
          do i =1, npt
            do l =1,3
              coor(l,i)=r_initial(l,i)
            end do
          end do

          call eforce()

          dU2 = 0.d0
          do i=1,npt
            do l=1,3
              dU2 = dU2 + dpot(l,i)**2
            end do
          end do
          dStmp_initial = dsqrt(Hamilt + dU2) 
        end if

        if (last) then
          do i =1, npt
            do l =1,3
              coor(l,i)=r_final(l,i)
            end do
          end do
          call eforce()

          dU2 = 0.d0
          do i=1,npt
            do l=1,3
              dU2 = dU2 + dpot(l,i)**2
            end do
          end do
          dStmp_final = dsqrt(Hamilt + dU2)
        end if 

      return
      end
