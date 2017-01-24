      subroutine AddCenter(j,Csize,expandNode)

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/SEARCH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

      double precision rms
      integer Csize,i,k,l,kk,j, expandNode

      ! could be made more effectively by considering
      ! global set of relevant centers within 2*cutoff from the
      ! expanded node
      k = (j-1)*npt
      do i = 1, npt
        coor(1,i) = structures(1,k+i)
        coor(2,i) = structures(2,k+i)
        coor(3,i) = structures(3,k+i)
      end do
     
      do i = 1,Csize
        call Distance(coor(1,1),centers(1,(i-1)*npt+1),rms)
        if (rms. lt. cutoff) return
      end do

      write(6,*)"Adding center."

      Csize = Csize + 1
      kk = (Csize-1)*npt
      do k = 1,npt
        do l = 1, 3
          centers(l,kk+k) = coor(l,k)
        end do
      end do
      ord(Csize) = Csize
      nodeType(Csize) = nodeType(expandNode)
      return

      end
