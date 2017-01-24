      integer FUNCTION overlappair_idx(i,j) 
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
c args
      integer i, j
c locals
      integer tempi, tempj, temp, jstart, jend, jj
c begin
      tempi = i
      tempj = j
      if ( tempj .eq. tempi ) then
          overlappair_idx = 0
	  return
      end if

      if ( tempj .le. tempi ) then
          temp = tempi
	  tempi = tempj
	  tempj = temp
      end if

      jstart = overlappair_ptrs(tempi)
      jend = overlappair_ptrs(tempi+1) - 1
      if ( jend .lt. jstart ) then
          overlappair_idx = 0
	  return
      end if 

      do jj = jstart, jend 
          if ( overlappair(jj) .eq. tempj ) then
	     overlappair_idx = jj
	     return
          end if
      end do
	overlappair_idx = 0

      return
      END 
