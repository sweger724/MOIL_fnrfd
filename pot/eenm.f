      subroutine eenm()
        implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ELASTIC.BLOCK'

      integer i

C In order to speed upt the computation it is presumed that these global 
C  ELASTIC.BLOCK parameters were precomputed before calling eenm()
        call network(U1,gradU1,g1list1,g1list2,req1,maxU,.false.)
C        write(6,*)"UUU:",U1
        limU = max(U1,100*enm_beta)
        call network(U2,gradU2,g2list1,g2list2,req2,limU,.true.)
C        write(6,*)"UU2:",U2
        call G_deriv(U1,U2,enm_beta)

        e_enm = 0.5d0*( U1+U2 - dsqrt( (U1-U2)**2 + 4.d0 *enm_beta**2 ))

        do 10 i=1,npt
            dpot(1,i) = dpot(1,i) + 
     &                  dG_dU1 * gradU1(1,i) + dG_dU2 * gradU2(1,i)
            dpot(2,i) = dpot(2,i) + 
     &                  dG_dU1 * gradU1(2,i) + dG_dU2 * gradU2(2,i)
            dpot(3,i) = dpot(3,i) + 
     &                  dG_dU1 * gradU1(3,i) + dG_dU2 * gradU2(3,i)
10      continue 
        return
        end
