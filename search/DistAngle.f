        function DistAngle (x1,x2,N)
        
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SEARCH.BLOCK'

        double precision x1(max_reducedCoor), x2(max_reducedCoor)
        double precision d1, d2, d3, pi2, DistAngle

        double precision tmp, deltax
        integer i, N

         pi2 = 3.14159265 * 2.d0

         DistAngle = 0.d0

         do i = 1, N

           deltax = x1(i) - x2(i)
           d1 = abs (deltax)
           d2 = abs (deltax - pi2)
           d3 = abs (deltax + pi2)
           tmp = min(d1,d2,d3)
           tmp = tor_weight(i) * tmp
           DistAngle = DistAngle + tmp**2
         end do
         
        end
