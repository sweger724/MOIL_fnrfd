	subroutine comc(coor,grdp,divms,natom,pointr,
     1  nselec,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,dmass1,dmass2,
     2  coor1,sigma,orth,debug)
c
c calculate gradient of rigid body constraints
c and then orthnormalize the constraints vectors.
c
c *** input
c coor     - the coordinates of the reference system for which
c            center of mass translation and rotation are compared
c grdp     - the path slope at the position x0,y0,z0
c divms    - double precision vector of 1.d0 over masses
c orth     - boolean, if true, 7 constraints are used, not 6.
c debug    - logical. true=print a lot of debugging information
c
c *** output
c constraints value at initial configuration, should remain
c      		constants.
c grdcmx - gradient of the translation of center of mass (X)
c grdcmy                   and Y
c grdcmz                   and Z
c
c grdlx  - gradient of rotation (X)
c grdly                         (Y)
c grdlz                         (Z)
c
c NOTE THAT THE VECTORS ABOVE UNDERWENT ORTHONORMALIZTION WITH
c RESPCT TO EACH OTHER.
c
	implicit none

	include 'COMMON/UNITS.BLOCK'

	double precision coor(3,*)
	double precision grdp(3,*)
	double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
	double precision grdcmx(3,*),grdcmy(3,*),grdcmz(3,*)
	double precision sigma(*)
	double precision divms(*),dmass1(*),dmass2(*),coor1(3,*)
	integer natom,nselec,nc,ncc
	integer pointr(*)
	logical orth, debug
c
c	local
c
	double precision norm,test(7)
	integer i,j,k


	if (orth) then
	   nc = 7
	   ncc = 7
	else
	   nc = 6
	   ncc = 3
	end if
	
c       
c The zero order values of the constraints (which are not zero)
c are stored in sigma(1-7)
c sigma(1-3) - center of mass positions for x y z
c sigma(4-6) - Orientation value for x y & z ( zero before
c           orthonormalization)
c sigma(7) - (if orth) scalar product of path slope and initial
c            coordinates

	do 11 i=1,nc
	   sigma(i)=0.d0
 11	continue
c
c calculate double precision mass and picking the
c selected atoms to a separate vector
c

c       NB: the "deselection" in coor1 is no longer needed,
c       since it is done in getmlst.f. -AW. 6-Dec-2007.
	do 1 j=1,nselec
	   i=pointr(j)
	   dmass1(j) =1.d0/divms(i)
	   dmass2(j) =divms(i)
	   do 87 k = 1,3
c	      coor1(k,j)=coor(k,i)
	      coor1(k,j)=coor(k,j)
 87	   continue
1	continue


c calculate the gradient of the translation vector.
c Note the strange normalization: If a1(i) and a2(i) are
c vectors of constraint gradients and m(i) is the mass,
c the scalar product is defined by
c SUM a1(i)*a2(i)/m(i) 
c
	norm=0.d0
	do 2 j=1,nselec
	   grdcmx(1,j)=dmass1(j)
	   grdcmx(2,j)=0.d0
	   grdcmx(3,j)=0.d0
	   
	   grdcmy(1,j)=0.d0
	   grdcmy(2,j)=dmass1(j)
	   grdcmy(3,j)=0.d0
	   
	   grdcmz(1,j)=0.d0
	   grdcmz(2,j)=0.d0
	   grdcmz(3,j)=dmass1(j)

	   grdlx(1,j)=0.d0
	   grdlx(2,j)= dmass1(j)*coor1(3,j)
	   grdlx(3,j)=-dmass1(j)*coor1(2,j)

	   grdly(1,j)=-grdlx(2,j)
	   grdly(2,j)=0.d0
	   grdly(3,j)=dmass1(j)*coor1(1,j)

	   grdlz(1,j)=-grdlx(3,j)
	   grdlz(2,j)=-grdly(3,j)
	   grdlz(3,j)=0.d0
 2	continue


	do 25 j=1,nselec
c          actually it is: norm=norm+grdcm[xyz](i)*grdcm[xyz](i)/dmass(i)
	   norm = norm + dmass1(j)

c          opportunity for constraints' values calculation
	   do 65 k=1,3
	      sigma(k) = sigma(k) + dmass1(j)* coor1(k,j)
 65	   continue

c          What about sigma(4-6)? -TF, 8-Mar-2005
	   if (orth) then
	      do 665 k=1,3
		 sigma(7) = sigma(7) + coor1(k,j)*grdp(k,j)
 665	      continue
	   end if
 25	continue



c
c normalize grdcmx grdcmy grdcmz and the constraint functions
c
	norm = 1.d0 / dsqrt(norm)

	do 54 k = 1,3
	   sigma(k) = sigma(k) * norm
 54	continue

	do 3 j=1,nselec
	   do 4 k=1,3
	      grdcmx(k,j) = grdcmx(k,j) * norm
	      grdcmy(k,j) = grdcmy(k,j) * norm
	      grdcmz(k,j) = grdcmz(k,j) * norm
 4	   continue
 3	continue
	
	if (debug) then
	   write(stdo,*) ' grdp '
	   write(stdo,1000)((grdp(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdcmx '
	   write(stdo,1000)((grdcmx(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdcmy '
	   write(stdo,1000)((grdcmy(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdcmz '
	   write(stdo,1000)((grdcmz(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdlx '
	   write(stdo,1000)((grdlx(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdly '
	   write(stdo,1000)((grdly(k,i),k=1,3),i=1,nselec)
	   write(stdo,*) ' grdlz '
	   write(stdo,1000)((grdlz(k,i),k=1,3),i=1,nselec)
c 1000	   format(1x,3(f9.6))
 1000	   format(1x,3(f15.6))

 1001	   format(20X,3(f15.6))
	   
	   do 31 i=1,nc
	      test(i)=0.d0
 31	   continue
	   do 32 j=1,nselec
	      do 33 k=1,3
		 test(1)=test(1)+grdcmx(k,j)*grdcmx(k,j)*dmass2(j)
		 test(2)=test(2)+grdcmy(k,j)*grdcmy(k,j)*dmass2(j)
		 test(3)=test(3)+grdcmz(k,j)*grdcmz(k,j)*dmass2(j)
		 test(4)=test(4)+grdcmx(k,j)* grdlx(k,j)*dmass2(j)
		 test(5)=test(5)+grdcmx(k,j)* grdly(k,j)*dmass2(j)
		 test(6)=test(6)+grdcmx(k,j)* grdlz(k,j)*dmass2(j)
		 if (orth) then
		    test(7)=test(7)+grdcmx(k,j)*grdp(k,j)*dmass2(j)
		 end if
 33	      continue
 32	   continue
	   write(stdo,*)' before orthog. norms of grdcm[xyz]'
	   write(stdo,*)' scalar prod. of grdcmx & grdl[x-z]'
	   if (orth) write(stdo,*) 'grdp'
	   write(stdo,*)(test(i),i=1,nc)
	end if
c       
c       orthogonalize the vectors with mass weighting. Obviously after
c       orthonormaliztion, the meaning of rotation in a given direction
c       may be lost.

c ** grdlx
	call orthg(grdcmx,grdlx,nselec,dmass2,sigma,1,4)
	call orthg(grdcmy,grdlx,nselec,dmass2,sigma,2,4)
	call orthg(grdcmz,grdlx,nselec,dmass2,sigma,3,4)
c ** grdly
	call orthg(grdcmx,grdly,nselec,dmass2,sigma,1,5)
	call orthg(grdcmy,grdly,nselec,dmass2,sigma,2,5)
	call orthg(grdcmz,grdly,nselec,dmass2,sigma,3,5)
	call orthg(grdlx ,grdly,nselec,dmass2,sigma,4,5)
c ** grdlz
	call orthg(grdcmx,grdlz,nselec,dmass2,sigma,1,6)
	call orthg(grdcmy,grdlz,nselec,dmass2,sigma,2,6)
	call orthg(grdcmz,grdlz,nselec,dmass2,sigma,3,6)
	call orthg(grdlx ,grdlz,nselec,dmass2,sigma,4,6)
	call orthg(grdly ,grdlz,nselec,dmass2,sigma,5,6)
c ** grdp: note order of args 1, 2, 6, and 7. -TF 8-Mar-2005
	if (orth) then
	   call orthg(grdcmx,grdp,nselec,dmass2,sigma,1,7)
	   call orthg(grdcmy,grdp,nselec,dmass2,sigma,2,7)
	   call orthg(grdcmz,grdp,nselec,dmass2,sigma,3,7)
	   call orthg(grdlx ,grdp,nselec,dmass2,sigma,4,7)
	   call orthg(grdly ,grdp,nselec,dmass2,sigma,5,7)
	   call orthg(grdlz ,grdp,nselec,dmass2,sigma,6,7)
	end if



c
c check that everything is orthonormal
c
	if (debug) then
	   do 61 i=1,nc
	      test(i)=0.d0
 61	   continue
	   do 7 j=1,nselec
	      do 77 k=1,3
		 test(1)=test(1)+grdcmx(k,j)*grdcmx(k,j)*dmass2(j)
		 test(2)=test(2)+grdcmy(k,j)*grdcmy(k,j)*dmass2(j)
		 test(3)=test(3)+grdcmz(k,j)*grdcmz(k,j)*dmass2(j)
		 test(4)=test(4)+ grdlx(k,j)* grdlx(k,j)*dmass2(j)
		 test(5)=test(5)+ grdly(k,j)* grdly(k,j)*dmass2(j)
		 test(6)=test(6)+ grdlz(k,j)* grdlz(k,j)*dmass2(j)
		 if (orth) then
		    test(7)=test(7)+grdp(k,j)*grdp(k,j)*dmass2(j)
		 end if
 77	      continue
 7	   continue
	   write(*,*)'normalization of grdcm[x,y,z] grdl[x,y,z]'
	   if (orth) write(*,*) 'and grdp'
	   write(*,*) (test(i),i=1,nc)
	   
	   do 8 i=1,ncc
	      test(i)=0.d0
 8	   continue
	   do 9 j=1,nselec
	      do 90 k=1,3
		 test(1)=test(1)+grdlx(k,j)*grdcmx(k,j)*dmass2(j)
		 test(2)=test(2)+grdlx(k,j)*grdcmy(k,j)*dmass2(j)
		 test(3)=test(3)+grdlx(k,j)*grdcmz(k,j)*dmass2(j)
		 if (orth) then
		    test(4)=test(4)+grdlx(k,j)*  grdp(k,j)*dmass2(j)
		    test(5)=test(5)+ grdp(k,j)*grdcmx(k,j)*dmass2(j)
		    test(6)=test(6)+ grdp(k,j)*grdcmy(k,j)*dmass2(j)
		    test(7)=test(7)+ grdp(k,j)*grdcmz(k,j)*dmass2(j)
		 end if
 90	      continue
 9	   continue
	   write(*,*)'scalar products (should be zero) of:'
	   write(*,*)'grdlx-grdcmx grdlx-grdcmy grdlx-grdcmz'
	   if (orth) then
	      write(*,*)'grdlx-grdp grdp-grdcmx'
	      write(*,*)'grdp-grdcmy grdp-grdcmz'
	   end if
	   write(*,*)(test(i),i=1,ncc)
	end if
c       --- end of debug statement
	
	   return
	   end
