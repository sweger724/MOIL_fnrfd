	subroutine Find_Rotation(npt,coor1,coor2,rotat)
        implicit none
c
C  Calculates rotation between coor1, coor2 of npt particles
C  it does not change positions of coors, it just returns the 
C  rotation matrix. It assumes that both coors are centred to [0,0,0]
C
C
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

	double precision coor1(3,*),coor2(3,*),rotat(3,3)
	integer npt

c local
	double precision kabs(3,3),kabs2(3,3),a(3),e(3),b(3,3)
	double precision allmass,det,norm
	double precision rms
	integer level,error
	integer i,j,k,l,namel
	character*16 name
	
	namel = 16
	name  = 'Find_Rotation'

c
	do 1 i=1,3
	 a(i) = 0.d0
	 e(i) = 0.d0
	 do 1 j=1,3
	  kabs(i,j)  = 0.d0
	  kabs2(i,j) = 0.d0
	  b(i,j)     = 0.d0
1	continue
	  

c calculate the kabsch matrix
	
	do 6 l=1,npt
	   i = l
	 do 6 k=1,3
	  do 6 j = 1,3
           kabs(j,k) = kabs(j,k) + coor1(j,i)*coor2(k,i)
6	continue

	if (.false.) then
         write(stdo,*)' kabs matrix '
         do i=1,3
          write(stdo,*)(kabs(i,j),j=1,3)
         end do
	 write(stdo,*)' coor1 matrix '
	 do i=1,3
          write(stdo,*)(coor1(i,j),j=1,3)
         end do
	 write(stdo,*)' coor2 matrix '
	 do i=1,3
          write(stdo,*)(coor2(i,j),j=1,3)
         end do
        end if

c check that the detrminant is not zero

	det = kabs(1,1)*kabs(2,2)*kabs(3,3) - kabs(1,1)*kabs(2,3)
     1		*kabs(3,2) - kabs(1,2)*kabs(2,1)*kabs(3,3) +
     2		kabs(1,2)*kabs(2,3)*kabs(3,1) +
     3		kabs(1,3)*kabs(2,1)*kabs(3,2) -
     4		kabs(1,3)*kabs(2,2)*kabs(3,1)
C	if (dabs(det).lt.1.d-10) then
C	 level = 1
C	 write(stdo,*)' det = ',det
C	 call alert(name,namel,' Zero determinant ',18,level)
C	end if

c construct a positive definite matrix by multiplying kabs
c by its transpose

	do 7 i=1,3
	 do 7 j=1,3
	  do 7 k=1,3
	   kabs2(i,j) = kabs2(i,j) + kabs(k,i)*kabs(k,j)
7	continue

c diagonalize kabs2 (matrix RtR in the above reference)
c on output kabs2 is destroyed and the eigenvectors are provided
c in kabs2 - kabs2(i=1,3 ; j) is the j-th eigenvector.
c "a" includes the eigenvalues, "e" is used 
c as a temporary vector
      
	call house(kabs2,3,3,a,e,ERROR)

	if (error.gt.0) then
	 write(stdo,100)error
100	 format(1x,' Error in diagonalization, error para ',i5)
	 level = 1
	 call alert(name,namel,' House failed ',14,level)
	end if

	do 11 j=1,3

c generate the so called b vectors
	 norm = 1.d0/dsqrt(a(j))
	 do 10 i=1,3
	  do 10 k=1,3
	   b(i,j) = b(i,j) + kabs(i,k)*kabs2(k,j)*norm
10	 continue
11	continue
c
c   avijit fix (04/12/2000)
c
      if (a(3).lt.1.0d-4) then
         a(3) = 0
C         write (6,*) 'zeroing out 3rd eigenvalue'
         kabs2(1,3) = kabs2(2,1) * kabs2(3,2) - 
     $        kabs2(2,2) * kabs2(3,1)
         kabs2(2,3) = kabs2(3,1) * kabs2(1,2) - 
     $        kabs2(3,2) * kabs2(1,1)
         kabs2(3,3) = kabs2(1,1) * kabs2(2,2) -
     $        kabs2(1,2) * kabs2(2,1)
         
         b(1,3) = b(2,1) * b(3,2) - b(2,2) * b(3,1)
         b(2,3) = b(3,1) * b(1,2) - b(3,2) * b(1,1)
         b(3,3) = b(1,1) * b(2,2) - b(1,2) * b(2,1)
               
      end if

c calculate the rotation matrix rotat
114      continue

	do 115 i=1,3
	 do 115 j=1,3
	   rotat(i,j) = 0.d0
115      continue


	do 12 i=1,3
	 do 12 j=1,3
	  do 12 k=1,3
	   rotat(i,j) = rotat(i,j) + b(i,k)*kabs2(j,k)
12	continue

	 if (debug) then
	 write(stdo,*)' rotation matrix '
	 do 125 i=1,3
	  write(stdo,*)(rotat(i,j),j=1,3)
125	 continue
	end if
c check if the determinant is zero or negative

	det = rotat(1,1)*rotat(2,2)*rotat(3,3) - rotat(1,1)*rotat(2,3)
     1		*rotat(3,2) - rotat(1,2)*rotat(2,1)*rotat(3,3) +
     2		rotat(1,2)*rotat(2,3)*rotat(3,1) +
     3		rotat(1,3)*rotat(2,1)*rotat(3,2) -
     4		rotat(1,3)*rotat(2,2)*rotat(3,1)
	if (dabs(det).lt.1.d-10) then
	 level = 1
	 call alert(name,namel,' Zero determinant ',18,level)
	end if

	return
	end

