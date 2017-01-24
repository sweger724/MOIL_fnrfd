      subroutine norm_test(il,ii,E,s,udeb)

      implicit none

      integer ii,udeb,il
      double precision E,s,T

      T = (dble(ii))**(0.5)*
     &    ((s/(dble(ii))-(E/dble(ii))**2)/dble(ii))**(0.5)/
     &    (E/dble(ii))

      write (udeb,*) '- gaussian test -',il,ii,T

      return
      end
