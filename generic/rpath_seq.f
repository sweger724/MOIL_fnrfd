      subroutine rpath_seq(upath,nstru)
      implicit none
c
c subroutine to read PATH format coordinate files
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'

      integer upath,nstru

c local
      character*9 name
      integer namel,level,i,j,k
      double precision e0
      
      name  = 'rpath_seq'
      namel = 9

                 do i=1,npt
                 chain(i) = '*'
                 enddo


      read(upath,err=999,end=999)e0,((coor(j,i),i=1,npt),j=1,3)
      if (nstru.gt.0) then
         write(stdo,100)nstru,e0
 100     format(1x,' *** READING PATH FILE ',i5,' ENERGY = ',e15.5)
      end if
      return
c
 999  continue
      level = 1
      call alert(name,namel,'Error while reading Path file',29,level)
c
      return
      end
