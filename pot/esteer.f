      subroutine steer()
        implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/STEER.BLOCK'

c calculate center of mass constraint for selected particles.


      integer i,j,level,namel
      character*6 name
      logical first
      double precision invsteer
      double precision xcent,ycent,dx,dy
      data first/.true./
      save first,invsteer

      namel = 6
      name = 'esteer'
      if (nsteer.eq.0) then
       level =3
       call alert(name,namel,'Nothing to steer',16,level)
       return
      end if

      if (first) then
        invsteer = dble(1.d0/nsteer)
	first = .false.
        xcent = 0.d0
        ycent =0.d0
        do i=1,nsteer
         j = pick_steer_compress(i)
         xcent = xcent + coor(1,j)
         ycent = ycent + coor(2,j)
        end do
        xcent = xcent*invsteer
        ycent = ycent*invsteer
      end if
      
      do 100 i=1,nsteer
         j = pick_steer_compress(i)
         dx = coor(1,j)-xcent
         dy = coor(2,j)-ycent
         dpot(1,j) = dpot(1,j) + kxy*dx
         dpot(2,j) = dpot(2,j) + kxy*dy
         dpot(3,j) = dpot(3,j) + kz
         e_steer=e_steer + 0.5d0*(kxy*(dx*dx+dy*dy))+kz*coor(3,j)
100   continue
     

      return
      end
