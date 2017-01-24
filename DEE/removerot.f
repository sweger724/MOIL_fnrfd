      subroutine removerot(indpos,rotfirst,rotlast)
c
c Subroutine that actually removes the rotamers from the list. It works
c on the rotamers of a given enhanced position. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c
c      
      character*9 name
      integer namel,level
c
      integer indpos,j,rotlastreal
      integer rot1,rot1real,rotfirst,rotlast
      integer nmonenhstart
c
      name = 'removerot'
      namel = 9
      level=1
c
c
c     
      rot1=rotfirst
      rot1real=rotaux(rot1)
      rotlastreal=rotaux(rotlast)
c
      do while (rot1real.le.rotlastreal)
c     
c........check if it was removed in this iteration and then 
c........remove it from the auxiliary matrices. 
         if (.not.kept(rot1real)) then
c
            do j=indpos,nposenh
               poirotenhaux(j)=poirotenhaux(j)-1
            end do
c     
            do j=rot1+1,nmonenhleft
               rotaux(j-1)=rotaux(j)
            end do
c
            nmonenhleft=nmonenhleft-1
c     
         else 
c
            rot1=rot1+1
c     
         end if
c
         if (rot1.le.nmonenhleft) then
c
            rot1real=rotaux(rot1)
c     
         else
c
            rot1real=rotlastreal+1
c
         endif
c
      end do
c
c
      return
      end
c     
c
