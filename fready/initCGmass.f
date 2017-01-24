      subroutine initCGmass()
      
      implicit none
            
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c      local variables
      integer i,getCGid
c
c  Iniailize side chain masses
c

       CMmass(1) = 76.d0  + 6.d0
       CMmass(2) = 90.d0  + 11.d0
       CMmass(3) = 56.d0  + 2.d0
       CMmass(4) = 54.d0  + 4.d0
       CMmass(5) = 0.d0  
       CMmass(6) = 12.d0  + 3.d0
       CMmass(7) = 36.d0  + 6.d0
       CMmass(8) = 100.d0 + 7.d0
       CMmass(9) = 122.d0 + 8.d0
       CMmass(10) = 44.d0 + 3.d0
       CMmass(11) = 48.d0 + 9.d0
       CMmass(12) = 48.d0 + 9.d0
       CMmass(13) = 68.d0 + 7.d0 
       CMmass(14) = 36.d0 + 7.d0
       CMmass(15) = 62.d0 + 11.d0
       CMmass(16) = 68.d0 + 4.d0
       CMmass(17) = 66.d0 + 6.d0
       CMmass(18) = 28.d0 + 3.d0
       CMmass(19) = 40.d0 + 5.d0
       CMmass(20) = 84.d0 + 7.d0

       do i= 1,npt
         CGid(i)=getCGid(moname(poimon(i)))
         if (ptnm(i).eq."CA  ") then
           ptms(i) = 47.d0
         else
           ptms(i) = CMmass(CGid(i))
         endif
         invms(i) = 1.d0/ptms(i)
       enddo
        
       return

      end
