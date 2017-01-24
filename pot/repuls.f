        subroutine repuls(n)
        include 'COMMON/LENGTH.BLOCK'   
        include 'COMMON/ENERGY.BLOCK'  
        include 'COMMON/COORD.BLOCK'   
        include 'COMMON/CONNECT.BLOCK'   
        include 'COMMON/SPECL.BLOCK'   
        
        double precision tmp,tempe,tempx,tempy,tempz
        double precision temprx(maxmorsb),tempry(maxmorsb)
        double precision temprz(maxmorsb)
        double precision tempdx(maxpt),tempdy(maxpt),tempdz(maxpt)
        integer i,j,namel,level,n
        character*6 name
        data name/'repuls'/
        data namel/6/
        tote_repuls=0.d0
         e_repuls(n) = 0. d0
        if (Arep(n).le.0..or.beta1(n).le.0..or. Brep(n) .le.0.) then
          level = 1
          call alert(name,namel,'Arep or Brep or beta is 0.',17,level)
        end if
         i = imb1(n)
         j = imb2(n)
         tempx = coor(1,i) - coor(1,j)
         tempy = coor(2,i) - coor(2,j)
         tempz = coor(3,i) - coor(3,j)
         dist(n) = DSQRT(tempx*tempx + tempy*tempy + tempz*tempz)
         tmp = - beta1(n) * dist(n)
         tempe = DEXP(tmp)
         e_repuls(n) = Arep(n) * tempe - Brep(n)
         temprx(n) = tempx / dist(n)
         tempry(n) = tempy / dist(n)
         temprz(n) = tempz / dist(n)
         tempdx(i) = - beta1(n) * Arep(n) * tempe * temprx(n)
         tempdx(j) =   beta1(n) * Arep(n) * tempe * temprx(n)
         tempdy(i) = - beta1(n) * Arep(n) * tempe * tempry(n)
         tempdy(j) =   beta1(n) * Arep(n) * tempe * tempry(n)
         tempdz(i) = - beta1(n) * Arep(n) * tempe * temprz(n)
         tempdz(j) =   beta1(n) * Arep(n) * tempe * temprz(n)

         dpot(1,i) = dpot(1,i) + tempdx(i)
         dpot(1,j) = dpot(1,j) + tempdx(j)
         dpot(2,i) = dpot(2,i) + tempdy(i)
         dpot(2,j) = dpot(2,j) + tempdy(j)
         dpot(3,i) = dpot(3,i) + tempdz(i)
         dpot(3,j) = dpot(3,j) + tempdz(j)
        return
        end
        
