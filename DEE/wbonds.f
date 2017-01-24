      subroutine wbonds(uwcon)
c
c Write the connectivity list in terms of a bond list. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c      
      character*6 name
      integer namel,level
c
      integer uwcon,i
c
      name = 'wbonds'
      namel = 6
c
      write(uwcon,*) npt
      write(uwcon,*) nb
      do i=1,nb
         write(uwcon,*) ib1(i),ib2(i)
      end do
c
      return
      end
