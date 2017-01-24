        subroutine ovrlpck(coor,coor2,dpot,velo,jpick,iorie,rms)
        implicit none
c
c overlap coor2 with respect to coor
c such that their mass weighted rms is a minimum.
c for dynamics rotate the forces and velocities with the same rotation matrix
c In addition, the center of mass of both coordinate sets
c is set to zero
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
C        include 'COMMON/PARALLEL.BLOCK'

        integer jpick(*)
        double precision coor(3,*),coor2(3,*)
        double precision dpot(3,*),velo(3,*),rms
        integer iorie
c local
        integer i,namel
        character*7 name
        double precision tmp(3), tmp_mass(maxpt)
        
        namel = 7
        name  = 'ovrlpck'

        if (iorie .ne. npt) then 
          do i = 1, npt
            tmp_mass(i) = ptms(i)
            ptms(i) = 0.d0
          end do

          do i = 1, iorie
            ptms(jpick(i)) = tmp_mass(jpick(i))
          end do
        end if

        call rmsd_weight(npt,coor,coor2,rms,.false.,ptms)

        if (iorie .ne. npt) then
          do i = 1,npt
            ptms(i) = tmp_mass(i)
          end do
        end if

c rotate also the force vector dpot and velocities

        do i=1,npt

          tmp(1) =  rotat(1,1)*dpot(1,i) + rotat(1,2)*dpot(2,i)
     1          + rotat(1,3)*dpot(3,i)
          tmp(2) =  rotat(2,1)*dpot(1,i) + rotat(2,2)*dpot(2,i)
     1          + rotat(2,3)*dpot(3,i)
          tmp(3) =  rotat(3,1)*dpot(1,i) + rotat(3,2)*dpot(2,i)
     1          + rotat(3,3)*dpot(3,i)

         dpot(1,i) = tmp(1)
         dpot(2,i) = tmp(2)
         dpot(3,i) = tmp(3)

          tmp(1) =  rotat(1,1)*velo(1,i) + rotat(1,2)*velo(2,i)
     1          + rotat(1,3)*velo(3,i)
          tmp(2) =  rotat(2,1)*velo(1,i) + rotat(2,2)*velo(2,i)
     1          + rotat(2,3)*velo(3,i)
          tmp(3) =  rotat(3,1)*velo(1,i) + rotat(3,2)*velo(2,i)
     1          + rotat(3,3)*velo(3,i)

         velo(1,i) = tmp(1)
         velo(2,i) = tmp(2)
         velo(3,i) = tmp(3)
         
        end do

        return
        end
