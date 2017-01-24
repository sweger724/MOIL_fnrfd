        subroutine Get_Dls(pseg,npt,rr,dls,ipick)
c     calculate now distances between all i,i+1 pairs
c     and the corresponding energy terms

      implicit none

        include "COMMON/LENGTH.BLOCK"
        include "COMMON/SDEL.BLOCK"

      integer pseg,npt,ipick(*)
      double precision rr(3,*),dls(*)

      integer i,j,k,loop
      double precision tmpx

C       allmass=0.d0
C       do i=1, npt
C          allmass = allmass + 1.d0/massfac(i)**2
C       enddo

      do 56 j = 1,pseg+1
         dls(j)=0.d0
         k = (j-1) * npt
         do 55 i = 1,npt
C            write(6,*)"ipick(i)",i,ipick(i),npt
            if (ipick(i).gt.0) then
              do loop =1,3
                 tmpx = rr(loop,i+k+npt) - rr(loop,i+k)
C                 write(6,*)"rr",rr(loop,i+k)
                 dls(j) = dls(j) + tmpx**2
              end do
            end if
 55      continue

         dls(j) = dsqrt(dls(j))
C         write(6,*)"dls(j) ",j,dls(j)/dsqrt(dble(npt))
 56   continue

      return
      end
