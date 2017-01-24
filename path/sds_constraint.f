      subroutine sds_constraint(sener,d0,ipick,npt,clo)

      implicit none

c     a subroutine to calculate the action associated with PATH constarins.
c     ThePATH constraint is given by
c     
c     S =     gamma SUM  (d(i,i+1) - <d>)^2 
C             - clo SUM  log(d(i,i+1))
c     
c     *** Note - in the current implementation the selection does not work!

      integer nselec,npt,f
      double precision sener, d0(*)
      integer ipick(*)
      
c     
c     common block for COORDinates and potential ENERGY derivatives
c     
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/SDEL.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'


      double precision e1(LGRID)


      integer i,j,k,l
      double precision tmp

      double precision clo,eij(lgrid)
      integer loop, Nloop

c     
c     Initialize our variables 
      sener = 0.d0

      do  i = 1,(pseg+2)*npt
        do l=1,3
          dsall(l,i) = 0.d0
        end do 
      end do

      call Get_Dls(pseg,npt,r,d0,ipick)
C      call GetDave(dave,d0)
c     calculate contribution from nearest neighbours
      Nloop = pseg
      if (last) Nloop = pseg + 1
      tmp = 0.d0

      do j = 1,Nloop
         tmp  = tmp + (d0(j) - dave)**2
      end do
      sener = sener + gamma*tmp


C     include logarithm term:
      tmp  = 0.d0
      do j = 1,Nloop
         tmp =tmp -dlog(d0(j))
      enddo
      tmp = tmp * clo
      sener = sener +tmp

c     also save the exponents for derivatives calculations
      do j = 1, pseg+1
         e1(j) = 2.d0*gamma*(1.d0 - dave/d0(j))
      end do
      
c     For the logarithmic term:
      do j=1, pseg+1
         eij(j) = -clo/(d0(j)**2)
      enddo
c     
c     calculate derivatives from nearest neighbour distance terms.
      do j = 1,pseg+1
         k  = (j-1)*npt
         l  = k +npt
         do i = 1,npt
            if (ipick(i).gt.0) then
              do loop=1,3
                tmp = e1(j)*(r(loop,k+i)-r(loop,l+i))
                dsall(loop,i+k)      = dsall(loop,i+k)   + tmp
                dsall(loop,i+l)      = dsall(loop,i+l)   - tmp
                tmp = eij(j)*(r(loop,k+i)-r(loop,l+i))
                dsall(loop,i+k)      = dsall(loop,i+k)   + tmp
                dsall(loop,i+l)      = dsall(loop,i+l)   - tmp
              end do
            end if
         end do
      end do
                
      return
      end
