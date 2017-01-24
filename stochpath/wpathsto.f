      subroutine wpathsto(ucrd)
c
c.v1.0 (last changed 24/2/97)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      integer ucrd
c
      integer i,j,k
c

c.....rewinds the output file. This is done to have a cleaner 
c.....output, at the risk of loosing the work if the system 
c.....crashes while writting. NEED TO FLUSH AT THE END.
c
c 
      write(stdo,*) ' Rewind and writing unit', ucrd
c
      rewind ucrd
c
c.....write only the structures actualy used. 
c.....(bc1) Two points are kept fixed at each extrema but only 
c...........one is used (and printed).
c.....(bc2) Two points are kept fixed at each extrema and both 
c...........used and printed. 
c
      k=0
      if (first) then
         if (bc2) then
            write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1           (r(k+2+j),j=0,npt3-3,3),
     2           (r(k+3+j),j=0,npt3-3,3)
         end if
         k = k + npt3
         write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1        (r(k+2+j),j=0,npt3-3,3),
     2        (r(k+3+j),j=0,npt3-3,3)
      end if
c     
      k = npt3
      do 10 i=1,pseg
         k = k + npt3
         write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1        (r(k+2+j),j=0,npt3-3,3),
     2        (r(k+3+j),j=0,npt3-3,3)
 10   continue
c     
      if (last) then
         k = k + npt3
         write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1        (r(k+2+j),j=0,npt3-3,3),
     2        (r(k+3+j),j=0,npt3-3,3)
         if (bc2) then
            k = k + npt3
            write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1           (r(k+2+j),j=0,npt3-3,3),
     2           (r(k+3+j),j=0,npt3-3,3)
         end if
      end if
c     
c     
c FLUSH is a command to force writing to the disk even with a cost
c       performance. This makes sure that data will not be lost
c       however, the flush routine is machine dependent
c
      write(stdo,*)'closing unit', ucrd
      call close_open_bin(ucrd)
cout     call flush1(stdo)
c     
      return
      end
      
