        subroutine lan_zen()
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'

        integer n
        real x(1)
        logical first
        data first/.true./

        if (first) then
         do 1 n=1,nmb
          probabi(n) = .true.
1       continue
        end if

        do 9 n=1,nmb
         if (emyes(n)) then
          call ranlux(x,1)
          if (x(1).lt.photon_per_step .and. e_morseb(n).lt.-20.) then
            emyes(n) = .false.
            repyes(n) = .true.
            write(*,*)' Absorbing a photon! '
            go to 9
           end if
         end if
         temp(1) = coor(1,imb1(n)) - coor(1,imb2(n))
         temp(2) = coor(2,imb1(n)) - coor(2,imb2(n))
         temp(3) = coor(3,imb1(n)) - coor(3,imb2(n))
        R = DSQRT(temp(1)*temp(1) + temp(2)*temp(2) + temp(3)*temp(3))
     
          
        if ( (R .gt. Rcros1) .and. (R .lt. Rcros2) ) then
           if ( probabi(n) ) then
            call probability(n)
           end if
           probabi(n) = .false.
          else
             probabi(n) = .true.
          end if
9      continue

        return
        end
