      subroutine mshakvl(dt)
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      double precision dt
c     
      integer i,j,iter,iwaters
      integer io,ih1,ih2

      double precision roh1(3),roh2(3),rhh(3),error(3),lambda(3)

c     
c     ********************************************
c     at this point, all coords are next coords, all
c     velocities are UNCONSTRAINED next velocities
c     (but do include constraints from coord update)
c     ********************************************
c     
c     
      do 8 iwaters=1,nwaters
c     
c     loop over tip3 water molecules to constrain their velocities.
c     mdivision is assumed on.
c     
	 io  = dpoipt(idxtip3(iwaters))-2
	 ih1 = io + 1
	 ih2 = io + 2
c     
c     get the 3 bond vectors
c     
	 do 1 i=1,3
            roh1(i) = coor(i,io)  - coor(i,ih1)
            roh2(i) = coor(i,io)  - coor(i,ih2)
            rhh(i)  = coor(i,ih1) - coor(i,ih2)
 1       continue

c     
c     compute the error vector
c     

	 error(1) = 0.d0
	 error(2) = 0.d0
	 error(3) = 0.d0

	 do 2 i=1,3
            error(1) = error(1) + (velo(i,io)  - velo(i,ih1))*roh1(i)
            error(2) = error(2) + (velo(i,io)  - velo(i,ih2))*roh2(i)
            error(3) = error(3) + (velo(i,ih1) - velo(i,ih2))*rhh(i)
 2       continue

         do 7 iter=1,20

c     
c     compute Lagrange's multipliers
c     
            do 3 i=1,3
               lambda(i) = 0.d0
               do 3 j=1,3
                  lambda(i) = lambda(i) - cmat(i,j)*error(j)
 3             continue

c     
c     correct velocities
c     
               do 4 i=1,3
                  velo(i,io)  = velo(i,io)  +
     $                 invmo*(lambda(1)*roh1(i)+lambda(2)*roh2(i))
                  velo(i,ih1) = velo(i,ih1) +
     $                 invmh*(-lambda(1)*roh1(i)+lambda(3)*rhh(i))
                  velo(i,ih2) = velo(i,ih2) +
     $                 invmh*(-lambda(2)*roh2(i)-lambda(3)*rhh(i))
 4             continue

c     
c     check convergence
c     
               error(1) = 0.d0
               error(2) = 0.d0
               error(3) = 0.d0

               do 5 i=1,3
                  error(1) = error(1)+(velo(i,io)-velo(i,ih1))*roh1(i)
                  error(2) = error(2)+(velo(i,io)-velo(i,ih2))*roh2(i)
                  error(3) = error(3)+(velo(i,ih1)-velo(i,ih2))*rhh(i)
 5             continue

               do 6 i=1,3
                  if (error(i).gt.tolcons) then
                     go to 7
                  end if
 6             continue
               go to 8
 7          continue
            call alert('mshakvl',7,' Not converged ',24,1)
 8       continue
         return
         end
