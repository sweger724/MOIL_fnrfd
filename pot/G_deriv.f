      subroutine G_deriv(Uu1,Uu2,beta)
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ELASTIC.BLOCK'
          double precision b,s,s3, Uu1, Uu2, beta

          b = beta**2
          s = 1.d0 / dsqrt( (Uu1-Uu2)**2 + 4.d0 * b)
          s3 = s**3

          dG_dU1 = 0.5d0 * (1.d0 - (Uu1 - Uu2)*s)
          dG_dU2 = 0.5d0 * (1.d0 + (Uu1 - Uu2)*s)
          d2G_d2U1 = - 2.d0 * b * s3
          d2G_d2U2 = d2G_d2U1
          d2G_dU1dU2 = - d2G_d2U1

      return
      end

