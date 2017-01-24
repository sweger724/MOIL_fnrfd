C           
C  GetVelocities is part of old code, we might want to use it later
C      it is just data-keeping subroutine, computes the velovities
C      of a given frame


C   variable igrid should get the name change
C          
           subroutine GetVelocities(igrid,npt,rr,pp,ptms)

           implicit none

           integer igrid
           double precision rr(3,*),pp(3,*),npt,ptms(npt)

           include 'COMMON/LENGTH.BLOCK'
           
           double precision normvector, drvector(3,MAXPT)
           integer i,f,l,k,km1,kp1
           
C     Lets calculate each individual momentum from the position vectors
C     First lets get the norm of the position vector (dr ) and
C     position vec.

        do i =2,igrid+1
          km1 = (i-2) * npt
          k   = (i-1) * npt
          kp1 =  i    * npt 
                 
          normdrvector = 0.0d0
          do f = 1,npt
            do l=1,3
              drvector(l,f) = rr(l,kp1 + f) - rr(l,km1 + f)
              normdrvector = normdrvector+drvector(l,f)**2
            end do
          end do
          normdrvector = sqrt(normdrvector)

          do f = 1,npt
            do l=1,3
              velocities(l,k+f) =
     $        pp(i)*drvector(l,f)/( normdrvector * sqrt(ptms(f)) )
            end do
          end do
          
        end do
