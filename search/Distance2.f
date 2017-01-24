      subroutine Distance(c1,c2,rms)

      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'

      double precision c1(3,maxpt), c2(3,maxpt), rms

      integer i

        call rmsd_weight(npt,c2(1,1),c1(1,1),rms,.false.,mmm)
        
      return
      end
