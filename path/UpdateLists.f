        subroutine UpdateLists(rr,npt,pseg)

C    Right now, we assume just one set of nb lists for all structures
C    at given processor... this is generated with respect to the
C    middle-one porcessor
C
          implicit none

          integer npt, pseg
          double precision rr(3,*)

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/COORD.BLOCK'
          include 'COMMON/ENERGY.BLOCK'
          include 'COMMON/MASSFAC.BLOCK'

          integer j,k,i,l

          j = pseg/2 +1
               k = (j-1)*npt
               do  i = 1,npt
                 do l=1,3
                   coor(l,i) =  rr(l,k+i)*massfac(i)
                 end do
               end do
c
                if (esymyes) call squeeze()
                call nbondm()
                if (esymyes) call syminit()

        return
        end
