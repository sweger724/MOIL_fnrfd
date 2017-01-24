      subroutine radiusG()
      
      implicit none
            
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c      local variables
      integer i 
      double precision Cx,Cy,Cz,M,R,dx,dy,dz
c
c  Calculate center of mass
c
        Cx = 0.d0
        Cy = 0.d0
        Cz = 0.d0
         M = 0.d0
       do i= 1,npt
           Cx = Cx + coor(1,i) * ptms(i)
           Cy = Cy + coor(2,i) * ptms(i)
           Cz = Cz + coor(3,i) * ptms(i)

           M = M + ptms(i)
       end do

        Cx = Cx / M
        Cy = Cy / M
        Cz = Cz / M
        
        R = 0.d0 

C  Calculate radius of gyration
       do i= 1,npt
           dx = coor(1,i) - Cx
           dy = coor(2,i) - Cy
           dz = coor(3,i) - Cz
           R = R + ptms(i) * (dx**2 + dy**2 + dz**2)
       end do

        R = sqrt(R/M)
      
      write(6,*)"Radius of gyration: ",R
      write(6,*)"Center: ", Cx, Cy, Cz
     
      end
