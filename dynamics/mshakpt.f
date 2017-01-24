        subroutine mshakpt(dt,dt2)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
        double precision dt,dt2
c
        integer i,j,iwaters,iter
        integer io,ih1,ih2
        double precision tmp
        double precision roh1(3),roh2(3),rhh(3),error(3),lambda(3)
        double precision roh12,roh22,rh1h2
        logical success
        

        do 8 iwaters=1,nwaters
c
c loop over TIP3/SPC/SPCE/TIP4 water molecules
c mdivision is assumed to be "on"
c
                io  = dpoipt(idxtip3(iwaters))-2
                ih1 = io + 1
                ih2 = io + 2

c
c get the 3 bond vectors.
c
                do 1 i=1,3
                        roh1(i) = coor(i,io)  - coor(i,ih1)
                        roh2(i) = coor(i,io)  - coor(i,ih2)
                        rhh(i)  = coor(i,ih1) - coor(i,ih2)
1               continue

c
c start iteration over current water
c calculate current distances
c
                do 7 iter=1,200

                roh12 = 0.d0
                roh22 = 0.d0
                rh1h2 = 0.d0
                
c
c Here velo is the current step
c
                do 2 i=1,3
                        tmp   = roh1(i) + velo(i,io) - velo(i,ih1)
                        roh12 = roh12 + tmp*tmp
                        tmp   = roh2(i) + velo(i,io) - velo(i,ih2)
                        roh22 = roh22 + tmp*tmp
                        tmp   = rhh(i) + velo(i,ih1) - velo(i,ih2)
                        rh1h2 = rh1h2 + tmp*tmp
2               continue

c
c calculate the vector of distance errors
c
                error(1) = roh12 - reqoh2
                error(2) = roh22 - reqoh2
                error(3) = rh1h2 - reqhh2
c
c check for convergence (current water)
c
                do 3 i=1,3
                if (dabs(error(i)).gt.0.4.or.iter.eq.19) then
      write(6,*)"current ERROR",iwaters,iter,error(1),error(2),error(3)
c reconstruct water molecule
                                coor(1,ih1) = 9999.d0
                                coor(1,ih2) = 9999.d0
                                call pos(success,ih1)
                                call pos(success,ih2)
                                do 35 j=1,3
                                        velo(j,ih1) = 0.d0
                                        velo(j,ih2) = 0.d0
                                        velo(j,io)  = 0.d0
35                              continue
                                go to 31
                else if (dabs(error(i)).gt.tolcons) then
c@
c                       write(*,*)' iwater error ', iwaters,error
                        go to 4
                end if
3               continue
31              continue
c
c if this point was reached the constraints are satisfied
c for the present water. Go to the next water
c
                go to 8
4               continue
c
c compute the Lagrange's multipliers
c
                do 5 i=1,3
                 lambda(i) = 0.d0
                  do 5 j=1,3
                   lambda(i) = lambda(i) - cmat(i,j)*error(j)
5               continue

c
c correct current step
c
                do 6 i=1,3
                 velo(i,io)  = velo(i,io)  +
     1  0.5d0*invmo*( lambda(1)*roh1(i)+lambda(2)*roh2(i))
                 velo(i,ih1) = velo(i,ih1) +
     2  0.5d0*invmh*(-lambda(1)*roh1(i)+lambda(3)*rhh(i))
                 velo(i,ih2) = velo(i,ih2) +
     3     0.5d0*invmh*(-lambda(2)*roh2(i)-lambda(3)*rhh(i))
6               continue

7       continue

c
c If this point was reached the constraint are NOT satisfied after
c maximum number of iterations. Call alert.
c
        call alert('mshakpt',7,'Mshak not converged!',20,1)

8       continue

        return
        end
