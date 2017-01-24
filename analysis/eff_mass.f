      program effms
c     
c     calculate effective mass along a reaction coordinate stored
c     in a PATH file
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/CCRD.BLOCK'
c     of     - integer function for Opening a File, returned value is the
c     assigned unit number
c     geti   - integer function to get integer value from a command line
c     getd   - double precision function to get DP value from a line
c     nstru  - number of structures in dynamics file
c     namel  - length of program name
c     level  - level of error found, level=0 only warning will be issued,
c     level=1 program stop
c     urcrd,ucon - units of dynamics coord and connectivity files
c     name - name of program (character) = contact
c     find - find a charcter in line (logical function)
c     fopen - check if file is open (logical function
c     pickpt - true if pick instruction was found in input
c     
      integer of,geti,nstru
      integer namel,i,j,k,l,level
      integer urcrd,ucon
      double precision getd

      integer error,size
      parameter (size=9*maxpt2d*maxpt2d)
      double precision invmass(size)
      double precision inv2(size)
      double precision testmp(1)
      double precision store(size)
      double precision vec1(3*maxpt2d), vec2(3*maxpt2d)
      double precision step,rc_mass,tmp,norm
      character*5 name
      logical find,fopen
      logical pickpt
      data ucon,urcrd/2*99/

      lpstr = 1
      norew = .false.

      stdi=5
      stdo=6
      totmon=0
      npt=0
      name='effms'
      namel=5
c     open junk file for rline
c     
      jnkf=25
      open(unit=jnkf,status='scratch')
c     default parameters
      nstru= 1
      step = 0.001d0
      pickpt=.false.

 1    continue
      call rline(name,namel,stdi)
      if (find('norw')) norew=.true.
      if (find('file')) then
         if (find ('conn')) then  
            ucon=of()
*     get connectivity 
            call rconn(ucon)
         end if
	 if (maxpt2d.lt.npt*npt) then
		write(*,*)' maxpt2d too small maxpt2d ',maxpt2d
		write(*,*) ' npt*npt = ',npt*npt
		level = 1
         	call alert(name,namel,'maxpt2d too small',17,level)
		stop
	 end if
		
*     get coordinate file
         if (find ('rcrd')) then
            if (npt.eq.0) then
               level = 1
               call alert(name,namel,'Must read con file first',
     1              24,level)
            end if
            urcrd=of()
         end if
      else 
         nstru=geti('#str',nstru)
         if (find ('action')) goto  3
      end if
      go to 1
 3    continue
c     initialized the nofreez vector
c     
      if (.not. fopen(ucon)) then
         level=1
         call alert(name,namel,'ucon not opened',15,level)
      else if (.not. fopen(urcrd)) then
         level=1
         call alert(name,namel,'urcrd not opened',16,level)
      end if
      
*     read paths structures
      rewind urcrd
      write(stdo,*) ' structure # effective mass '
      do 15 l=1,nstru
         do 2 i=1,size
            invmass(i)=0.d0
 2       continue
         rewind urcrd
         if (l.eq.1) then
            call rpath(urcrd,1)
            do 4 i=1,npt
               coor2(1,i) = coor(1,i)
               coor2(2,i) = coor(2,i)
               coor2(3,i) = coor(3,i)
 4          continue
            call rpath(urcrd,2)
            do 5 i=1,npt
               coor(1,i) = coor(1,i) - coor2(1,i)
               coor(2,i) = coor(2,i) - coor2(2,i)
               coor(3,i) = coor(3,i) - coor2(3,i)
 5          continue		
         else if (l.eq.nstru) then
            call rpath(urcrd,l-1)
            do 6 i=1,npt
               coor2(1,i) = coor(1,i)
               coor2(2,i) = coor(2,i)
               coor2(3,i) = coor(3,i)
 6          continue
            call rpath(urcrd,l)
            do 7 i=1,npt
               coor(1,i) = coor(1,i) - coor2(1,i)
               coor(2,i) = coor(2,i) - coor2(2,i)
               coor(3,i) = coor(3,i) - coor2(3,i)
 7          continue 
         else
            call rpath(urcrd,l-1)
            do 8 i=1,npt
               coor2(1,i) = coor(1,i)
               coor2(2,i) = coor(2,i)
               coor2(3,i) = coor(3,i)
 8          continue
            call rpath(urcrd,l+1)
            do 9 i=1,npt
               coor(1,i) = coor(1,i) - coor2(1,i)
               coor(2,i) = coor(2,i) - coor2(2,i)
               coor(3,i) = coor(3,i) - coor2(3,i)
 9          continue 
         end if
         norm = 0.d0
         do 11 i=1,npt
            norm = norm + coor(1,i)*coor(1,i) + coor(2,i)*coor(2,i)
     1           + coor(3,i)*coor(3,i)
 11      continue
         norm = 1.d0/dsqrt(norm)
         do 12 i=1,npt
            coor(1,i) = coor(1,i)*norm
            coor(2,i) = coor(2,i)*norm
            coor(3,i) = coor(3,i)*norm
            tmp = coor(1,i)*coor(1,i) + coor(2,i)*coor(2,i)
     1           + coor(3,i)*coor(3,i)
            if (dabs(tmp-1.d0).lt.1.d-12) then
               write(*,*)' DIAGONAL MASS MATRIX - RESULTS TRIVIAL'
c     stop
            end if
 12      continue
c     copy the path unit vector (currently stored in coor) to 
c     the first column and first row of the matrix - mass
         call vdcopy(coor,invmass,3*npt)
c     write(*,*)' after vdcopy '
c     call wmat(invmass,3*npt)
c     generate the rest of the matrix using a set of vectors
c     linearly independent of the rest of the vectors
         call makmat(invmass,3*npt)
c     write(*,*)' after makmat '
c     call wmat(invmass,3*npt)
c     orthonormalize the rest of the vectors
c     
         call orthmat(invmass,3*npt)
c     write(*,*)' after orthmat '
c     call wmat(invmass,3*npt)
c     create the off diagonal mass matrix 
c     
         call offdm(invmass,store,inv2,ptms,npt)
c     write(*,*)' after offdm '
c     write(*,*)' ptms = ',(ptms(i),i=1,npt)
c     call wmat(store,3*npt)
c     diagonalize the off-diagonal mass matrix as a way to generate
c     the inverse 
c     
         call house(store,3*npt,3*npt,vec1,vec2,error)
c     write(*,*)' error = ',error
c     write(*,*)' eigenvalues = ',(vec1(j),j=1,3*npt)
c     write(*,*)' after house #1 '
c     call wmat(store,3*npt)
c     use the new eigenvector and eigenvalues to generate a matrix which is
c     an inverse to the offdm 
c     
         call invmat(store,invmass,testmp,vec1,3*npt)
c     write(*,*)' after invmat '
c     call wmat(invmass,3*npt)
c     take the current matrix from 2 to 3*npt and diagonalize it
c     as a starting point for finding the inverse of the submatrix.
c     

c     (1) Copy invmass(2->3*npt) to store(1->3*npt-1)
         call shftmat(store,invmass,3*npt-1,1)
c     call shftmat(testmp,invmass,3*npt-1,1)
c     write(*,*)' after shftmat '
c     call wmat(store,3*npt-1)

c     (2) Diagonalize store
         call house(store,3*npt-1,3*npt-1,vec1,vec2,error)

c     (3) Compute the inverse of the submatrix of the inverse 
         call invmat(store,inv2,testmp,vec1,3*npt-1)

c     !!!! compute the effective mass !!!!
c     
         call effmass(rc_mass,inv2,invmass,3*npt)
         write(stdo,*)'** ',l,rc_mass
 15   continue
      
      stop
      end
