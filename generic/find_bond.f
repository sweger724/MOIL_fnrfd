        subroutine find_bond(iat1,iat2,ib)
c
c find bond index (ib) given two atomic indices iat1 & iat2
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'

        integer iat1,iat2,ib

c local
        integer ii1, jj1, i

        if (iat1.lt.iat2) then
                ii1 = iat1
                jj1 = iat2
        else
                ii1 = iat2
                jj1 = iat1
        end if

        do i=1,nb_all
                if (ii1.eq.ib1(i) .and. jj1.eq.ib2(i)) then
                        ib = i
                        return
                end if
        end do

        write(*,*)' iat1 iat2  = ',iat1,iat2
        call alert('find_bond',9,' Bond not found ',16,1)
        return
        end
