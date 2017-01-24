      subroutine draw_normal(output,ndegf,temper)
        implicit none
       include 'COMMON/CONVERT.BLOCK'
c
c initialize velocities according to Boltzmann distribution, assuming 
c the masses are 1.
c
      integer i
      integer ndegf
      double precision factor,temper,output(*)
      real gaussian


      factor = dsqrt(temper*kboltzmann)
      call gauss2(ndegf,output)
      do i=1,ndegf
        output(i) = output(i)*factor
      end do 
c     
      return
      end
      

      
      subroutine gauss2(N,v)

        implicit none

        include 'COMMON/LENGTH.BLOCK'

        integer N, i,i1,i2,m
        double precision v(*), pi
        real x(6*maxpt)
         
         m = N
         if (mod(m,2).eq.1) m=m+1
         
         pi = 3.14159265359
         call RANLUX(x,m)
         
         do i =1, (m/2)
           i1 = 2*i-1
           i2 = 2*i
           v(i1) = dsqrt(-2.d0 * log(x(i1)))*cos(2.d0*pi*x(i2))
           v(i2) = dsqrt(-2.d0 * log(x(i2)))*sin(2.d0*pi*x(i1))
         end do

        return 
        end
