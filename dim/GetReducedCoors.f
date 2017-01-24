      subroutine GetReducedCoors()

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SEARCH.BLOCK'
      include 'COMMON/COORD.BLOCK'

      integer i
      double precision CGtorsion

        do i = 1, Nreduced
          myReduced(i) = CGtorsion(torsions(i,1),torsions(i,2),
     &         torsions(i,3),torsions(i,4),coor)
        end do
        
      end


      subroutine GetReducedDerivs()

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SEARCH.BLOCK'
      include 'COMMON/COORD.BLOCK'

      integer i
      double precision D_one_dihedral

        do i = 1, Nreduced
          myReduced(i) = D_one_dihedral(torsions(i,1),torsions(i,2),
     &         torsions(i,3),torsions(i,4),coor,myDReduced(1,1,i))
        end do

      end
