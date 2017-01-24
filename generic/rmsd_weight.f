        subroutine rmsd_weight(Npart,coor,coor2,rms,centerCoor,pick)
        implicit none
c
c overlap coor2 with respect to coor (i.e. coor2 is moving!)
c such that their mass weighted rms is a minimum.
c In addition, the center of mass of both coordinate sets
c is set to zero
c Following Kabsch Acta. Cryst. A32,922(1976)
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'

        double precision coor(3,*),coor2(3,*),pick(*)
        logical centerCoor
        double precision rms
        integer Npart

c local
        double precision kabs(3,3),kabs2(3,3),a(3),e(3),b(3,3)
        double precision trash(3,3)
        double precision det,norm
        double precision tmp(3), wmonmass
        double precision small, allmass
        integer level,ERROR,ismall
        integer i,j,k,namel
        character*11 name
        
        namel = 11
        name  = 'rmsd_weight'
        small = 1.d20

c calculate center of mass for x,y,z and x1,y1,z1
c
c@
        ERROR = 0
        do 1 i=1,3
         a(i) = 0.d0
         e(i) = 0.d0
         do 1 j=1,3
          kabs(i,j)  = 0.d0
          kabs2(i,j) = 0.d0
          b(i,j)     = 0.d0
1       continue
          
        Acm(1)     = 0.d0
        Acm(2)     = 0.d0
        Acm(3)     = 0.d0
        Bcm(1)    = 0.d0
        Bcm(2)    = 0.d0
        Bcm(3)    = 0.d0
        allmass   = 0.d0
        do 2 i=1,Npart
         Acm(1)    = pick(i)*coor(1,i)  + Acm(1)
         Acm(2)    = pick(i)*coor(2,i)  + Acm(2)
         Acm(3)    = pick(i)*coor(3,i)  + Acm(3)
         Bcm(1)    = pick(i)*coor2(1,i) + Bcm(1)
         Bcm(2)    = pick(i)*coor2(2,i) + Bcm(2)
         Bcm(3)    = pick(i)*coor2(3,i) + Bcm(3)
         allmass = allmass + pick(i)
2       continue

        Acm(1)  = Acm(1)  / allmass
        Acm(2)  = Acm(2)  / allmass
        Acm(3)  = Acm(3)  / allmass
        Bcm(1) = Bcm(1) / allmass
        Bcm(2) = Bcm(2) / allmass
        Bcm(3) = Bcm(3) / allmass


c2

c set the center of mass positions to zero
        do 3 i=1,Npart
         coor(1,i)  = coor(1,i)  - Acm(1)
         coor(2,i)  = coor(2,i)  - Acm(2)
         coor(3,i)  = coor(3,i)  - Acm(3)
         coor2(1,i) = coor2(1,i) - Bcm(1)
         coor2(2,i) = coor2(2,i) - Bcm(2)
         coor2(3,i) = coor2(3,i) - Bcm(3)
3       continue

c calculate the kabsch matrix
        
        do 6 i=1,Npart
         do 6 k=1,3
          do 6 j=1,3
           kabs(j,k) = kabs(j,k) + pick(i)*coor(j,i)*coor2(k,i)
6       continue
c check that the detrminant is not zero

        det = kabs(1,1)*kabs(2,2)*kabs(3,3) - kabs(1,1)*kabs(2,3)
     1          *kabs(3,2) - kabs(1,2)*kabs(2,1)*kabs(3,3) +
     2          kabs(1,2)*kabs(2,3)*kabs(3,1) +
     3          kabs(1,3)*kabs(2,1)*kabs(3,2) -
     4          kabs(1,3)*kabs(2,2)*kabs(3,1)
        if (dabs(det).lt.1.d-10) then
         level = 1
         write(stdo,*)' det = ',det
         call alert(name,namel,' Zero determinant1 ',18,level)
        end if

c construct a positive definite matrix by multiplying kabs
c by its transpose

        do 7 i=1,3
         do 7 j=1,3
          do 7 k=1,3
           kabs2(i,j) = kabs2(i,j) + kabs(k,i)*kabs(k,j)
7       continue



c diagonalize kabs2 (matrix RtR in the above reference)
c on output kabs2 is destroyed and the eigenvectors are provided
c in kabs2 - kabs2(i=1,3 ; j) is the j-th eigenvector.
c "a" includes the eigenvalues, "e" is used 
c as a temporary vector
      
        call house(kabs2,3,3,a,e,ERROR)

        if (error.gt.0) then
         write(stdo,100)error
100      format(1x,' Error in diagonalization, error para ',i5)
         level = 1
         call alert(name,namel,' House failed ',14,level)
        end if

c sort the vectors such that a(1) > a(2)> a(3)
        do 8 i=1,3
         if (a(i).lt.small) then
          ismall = i
          small  = a(i)
         end if
8       continue


c generate the so called b vectors
         do 13 j=1,3
         norm = 1.d0/dsqrt(a(j))
         do 12 i=1,3
          do 12 k=1,3
           b(i,j) = b(i,j) + kabs(i,k)*kabs2(k,j)*norm
12       continue
13      continue
c
c   avijit fix (04/12/2000)
c
      if (a(3).lt.1.0d-4) then
         a(3) = 0
C         write (stdo,*) 'zeroing out 3rd eigenvalue'
         kabs2(1,3) = kabs2(2,1) * kabs2(3,2) - 
     $        kabs2(2,2) * kabs2(3,1)
         kabs2(2,3) = kabs2(3,1) * kabs2(1,2) - 
     $        kabs2(3,2) * kabs2(1,1)
         kabs2(3,3) = kabs2(1,1) * kabs2(2,2) -
     $        kabs2(1,2) * kabs2(2,1)
         
         b(1,3) = b(2,1) * b(3,2) - 
     $        b(2,2) * b(3,1)
         b(2,3) = b(3,1) * b(1,2) - 
     $        b(3,2) * b(1,1)
         b(3,3) = b(1,1) * b(2,2) -
     $        b(1,2) * b(2,1)
               
      end if

c calculate the rotation matrix rotat
14      continue

        do 15 i=1,3
         do 15 j=1,3
          rotat(i,j) = 0.d0
15      continue

        do 16 i=1,3
         do 16 j=1,3
          do 16 k=1,3
           rotat(i,j) = rotat(i,j) + b(i,k)*kabs2(j,k)
16      continue

c
c check that the rotation matrix is unitary
c
       do i=1,3
        do j=1,3
         trash(i,j) = 0
         do k=1,3
          trash(i,j) = trash(i,j) + rotat(i,k)*rotat(j,k)
         end do
        end do
       end do

      if (debug) then
         write(stdo,*)' rotation matrix '
         do 17 i=1,3
          write(stdo,*)(rotat(i,j),j=1,3)
17       continue
      end if


c check if the determinant is negative or zero

        det = rotat(1,1)*rotat(2,2)*rotat(3,3) - rotat(1,1)*rotat(2,3)
     1          *rotat(3,2) - rotat(1,2)*rotat(2,1)*rotat(3,3) +
     2          rotat(1,2)*rotat(2,3)*rotat(3,1) +
     3          rotat(1,3)*rotat(2,1)*rotat(3,2) -
     4          rotat(1,3)*rotat(2,2)*rotat(3,1)
        if (dabs(det).lt.1.d-10) then
         level = 1
         call alert(name,namel,' Zero determinant ',18,level)
        end if

c If the determinant is negative, invert all elements to 
c obtain rotation transformation (excluding inversion)

C       det = dabs(det)
        if (det.lt.0) then
         level = 0
         call alert(name,namel,' Negative determinant ',22,level)
         write(stdo,*) ' Negative det detected, reconstructing '
         do 18 i=1,3
          b(i,ismall) = -b(i,ismall)
18       continue
         go to 14
        end if

        if (dabs(det-1.d0).gt.3.d-1) then
         level = 1
         write(stdo,*) ' determinant of rotation matrix = ',det
         call alert(name,namel,'Nonunitary rotation ?!',22,level)
        end if

        det=(rotat(1,1)-1)**2 + (rotat(2,2)-1)**2 + (rotat(3,3)-1)**2
        det=det+rotat(2,1)**2+rotat(1,2)**2+rotat(2,3)**2
     &     +    rotat(3,2)**2+rotat(3,1)**2+rotat(1,3)**2
c        if (det.gt. 0.7) then 
c          write(6,*)"DET:",det
c          write(stdo,*)' rotation matrix '
c          do i=1,3
c            write(stdo,*)(rotat(i,j),j=1,3)
c          end do
c        end if

c rotate x1,y1,z1 
        do 20 i=1,Npart
          tmp(1) =  rotat(1,1)*coor2(1,i) + rotat(1,2)*coor2(2,i)
     1          + rotat(1,3)*coor2(3,i)
          tmp(2) =  rotat(2,1)*coor2(1,i) + rotat(2,2)*coor2(2,i)
     1          + rotat(2,3)*coor2(3,i)
          tmp(3) =  rotat(3,1)*coor2(1,i) + rotat(3,2)*coor2(2,i)
     1          + rotat(3,3)*coor2(3,i)
         coor2(1,i) = tmp(1)
         coor2(2,i) = tmp(2)
         coor2(3,i) = tmp(3)
20      continue        
C       return

c calculate current rms

        rms = 0.d0
        if (Npart .ne. npt) then 
        
          do i=1,Npart
            tmp(1) = coor2(1,i) - coor(1,i)
            tmp(2) = coor2(2,i) - coor(2,i)
            tmp(3) = coor2(3,i) - coor(3,i)
            rmsarray(i) = pick(i)*( tmp(1)**2 + tmp(2)**2 + tmp(3)**2 )
            rms = rms + rmsarray(i)
          end do

        else

          do j=1,totmon
            rmsarray(j) = 0.d0
            wmonmass = 0.d0
            do i=poipt(j-1)+1,poipt(j)
              tmp(1) = coor2(1,i) - coor(1,i)
              tmp(2) = coor2(2,i) - coor(2,i)
              tmp(3) = coor2(3,i) - coor(3,i)
              rmsarray(j) = rmsarray(j) + pick(i)*(tmp(1)**2 + tmp(2)**2
     $                    + tmp(3)**2 )
              wmonmass = wmonmass + pick(i)
            end do 
            rms=rms+rmsarray(j)
            if(wmonmass.gt.0.0d0) then
              rmsarray(j) = dsqrt(rmsarray(j)/wmonmass)
            endif
          end do

        end if

 300  continue

c       write(stdo,*) 'Npart:',Npart
C       return
        rms = dsqrt(rms/allmass)
c       write(stdo,*) 'rms:',rms
        
C  Return coor back to initial position
        if (centerCoor) then
          do i=1,Npart
            coor(1,i)  = coor(1,i) +  Acm(1)
            coor(2,i)  = coor(2,i) +  Acm(2)
            coor(3,i)  = coor(3,i) +  Acm(3)    
          end do
        end if

c        write(stdo,*) 'Center 1: ',Acm(1) ,Acm(2), Acm(3)
c        write(stdo,*) 'Center 2: ',Bcm(1) ,Bcm(2), Bcm(3) 
c        write(stdo,*) ' RMSD with homolog is = ',rms

        return
        end

