      subroutine str_velinit(strtemp)
c     
c     irand is provided from DYNA common BLOCK
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c     
c     initialize velocities according to Boltzmann distribution
c     using a gaussian distribution with a variance of 12/12=1
c     
c     ::::: -:- Modified from dynamics project by ZVA
c     
      double precision strtemp
      integer i,k,l
      double precision factor,tfactor
      real gaussian
c     
c      write(stdo,'(10x,a,f8.2)') 
c     >     "Start init fast velocities ! T_in=",strtemp
c     
      tfactor = dsqrt(strtemp*1.9878d-3)
c
c
c     
      do 10  k=1,inofrz
         i = nofreez(k)
         factor=tfactor/sqrt(ptms(i))
         velo(1,i) = gaussian(irand)*factor
         velo(2,i) = gaussian(irand)*factor
         velo(3,i) = gaussian(irand)*factor
 10   continue
c     
      call veleqtemp(strtemp)
c     
c      write(stdo,'(10x,a,f8.2)') 
c     >     "  End init fast velocities ! T_in=",strtemp
      return
      end
c     
      subroutine veleqtemp(strtemp)
c     
c     irand is provided from DYNA common BLOCK
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c     
c     ::::: -:- Modified from dynamics project by ZVA
c     
      double precision strtemp
      integer i,k,l
      double precision wkinene0,wkinene1,wscale
c     
      wkinene0 = 1.5d0 * inofrz * 1.9878d-3 * strtemp 
c     
c     calculate current temperature(s), and scale velocities
c     to obtain desired one(s)
c     
      wkinene1=0.0d0
c     
      do 413 k=1,inofrz
         l = nofreez(k)
         wkinene1 = wkinene1 +
     >        ptms(l)*(velo(1,l)**2+velo(2,l)**2+velo(3,l)**2)
 413  continue
c     
      wkinene1 = wkinene1 * 0.5d0
c     
c     write(stdo,'(10x,2(a,e10.3))') 
c     >     "Scale Kinetic energy=", wkinene1," to =>", wkinene0
c     
      wscale=sqrt(wkinene0/wkinene1)
      do 414 k=1,inofrz
         l = nofreez(k)
         velo(1,l) = velo(1,l)*wscale
         velo(2,l) = velo(2,l)*wscale
         velo(3,l) = velo(3,l)*wscale
 414  continue
c     
      return
      end
                
                
