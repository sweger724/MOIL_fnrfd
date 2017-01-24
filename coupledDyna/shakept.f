	subroutine shakept(maxit)
c
c correct bond distances to statisfy shake constraints
c epsilon - allowed average error in distance
c maxit   - maximum number of iteration
c
c	double precision epsilon
	integer maxit

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

C  Uses SHAKE algorithm to update coordinates and maintain
C  distance contstraints.  This implementation handles atoms with
C  zero mass gracefully and correctly.

C  INPUT --  nshak:                  number of shake constraints
C            ishak1, ishak2          pairs of atoms to apply shake constraints
C            dissq                   squares of shake distances
C            ptms                    atomic masses
C            x     , y     , z       coordinates after unconstrained move
C            xref ,  yref ,  zref    coordinates before unconstrained move
C            epsilon                 convergence criterion
C            maxit                 maximum number of iterations

C  OUTPUT -- velo  updated step


      double precision rx, ry, rz, r2
      double precision Mg, fac1, fac2
      double precision mass1(maxshak),mass2(maxshak)
      double precision ovrlx
      integer iat1, iat2, k, iterate, level
      integer istart
      logical first,again
      parameter (ovrlx=1.3d0)
      data first/.true./

      save first
      save mass1,mass2


      again = .false.
      istart = 1

      if (first) then
       first = .false.
       do 1 k=1,nshak
	iat1 = ishak1(k)
	iat2 = ishak2(k)
	mass1(k) = -ptms(iat2)/(ptms(iat1)+ptms(iat2))
	mass2(k) = 1.d0 + mass1(k)
1      continue
      end if

C  check each constraint at most maxit times
      do 1000 iterate=1,maxit
C          adjust each constraint in order
          do 100 k=istart,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)
C                compute difference vectors
              rx = cooref(1,k) + velo(1,iat1) -  velo(1,iat2)
              ry = cooref(2,k) + velo(2,iat1) -  velo(2,iat2)
              rz = cooref(3,k) + velo(3,iat1) -  velo(3,iat2)

C                  get update factors
              Mg = ((rx*rx + ry*ry + rz*rz) - dissq(k)) /
     1        (2.d0*(rx*cooref(1,k) + ry*cooref(2,k) + rz*cooref(3,k)))

	      fac1 = ovrlx*Mg*mass1(k)
	      fac2 = ovrlx*Mg*mass2(k)

C                  update coordinates along original difference vector
c@		if (k.eq.istart) then
c@			write(*,*)' ************* corrections '
c@			write(*,*)' k iat1 iat2 ',k,iat1,iat2
c@			write(*,*) ' Mg fac1 fac2 ',Mg,fac1,fac2
c@			write(*,*)' corrections = '
c@			write(*,*)' 1 = ',fac1*cooref(1,k)
c@			write(*,*)' 2 = ',fac2*cooref(1,k)
c@			write(*,*)' ************* corrections end'
c@		end if
	      velo(1,iat1) = velo(1,iat1) + fac1*cooref(1,k)
	      velo(2,iat1) = velo(2,iat1) + fac1*cooref(2,k)
	      velo(3,iat1) = velo(3,iat1) + fac1*cooref(3,k)

	      velo(1,iat2) = velo(1,iat2) + fac2*cooref(1,k)
	      velo(2,iat2) = velo(2,iat2) + fac2*cooref(2,k)
	      velo(3,iat2) = velo(3,iat2) + fac2*cooref(3,k)

100       continue
c chen
c	if (prll_on_off) then
c	    call put_my_V()
c	    call update_ex_shared()
c	    call get_others_V()
c	end if

C          check for convergence of all distances and return if ok
          do 200 k=istart,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)
              rx = cooref(1,k) + velo(1,iat1)-velo(1,iat2)
              ry = cooref(2,k) + velo(2,iat1)-velo(2,iat2)
              rz = cooref(3,k) + velo(3,iat1)-velo(3,iat2)
	      r2   = rx*rx + ry*ry + rz*rz - dissq(k)
              if (dabs(r2) .gt. shak_diff(k)) then
c chen
c		if (prll_on_off) call shak_unhappy()
c@		write(*,*)' iat1 iat2 diff',iat1,iat2,diff
c@		write(*,*)' rx ry rz ',rx,ry,rz
c@		write(*,*)' cooref = ',cooref(1,k),cooref(2,k),cooref(3,k)
c@		write(*,*)' iterate ',iterate
c@		write(*,*)' velo iat1 ',velo(1,iat1),velo(2,iat1),velo(3,iat1)
c@		write(*,*)' velo iat2 ',velo(1,iat2),velo(2,iat2),velo(3,iat2)
		istart = k
	        again = .true.
		go to 1000
	      end if
200       continue
	  if (again) then
	     again = .false.
	     istart= 1
c chen
c		if (prll_on_off) call shak_unhappy()
	     go to 1000
	  end if
c chen
c	if (prll_on_off) then
c	    call shak_happy()
c	    if (everybody_happy(1) .eq. 1) then
c		    call update_v()
c debug
c	write(stdo,*)'shakpt converged after ',iterate, ' iterations'


c       return
c       else
c		    go to 1000
c	    end if
c	end if
	     return
 1000	  continue

C
C If this point is reached, SHAKE did not converge! Do something about errors!
C
      level = 1
      call alert('shakept',7,' Unconverged after maxit steps',30,level)
      return
      end

c---------------------------------------------------------------------
c	subroutine     get_others_V()
c	The  subroutine copies the velocities (calculated in the previous iteration
c	by other PEs) of the shared particals.
c
c 	include 'COMMON/LENGTH.BLOCK'
c 	include 'COMMON/COORD.BLOCK'
c 	include 'COMMON/VELOC.BLOCK'
c 	include 'COMMON/CONNECT.BLOCK'
c 	include 'COMMON/UNITS.BLOCK'
c 	include 'COMMON/SHAKE.BLOCK'
c 	include 'COMMON/PARALLEL.BLOCK'
c 	include 'COMMON/DEBUG.BLOCK'


c 	integer i, pt,  pe, j
c         double precision  z
c 	do  1 i = 1, n_shared_total
c 	    pt = shared_list(i, 1)
c 	    pe = shared_list(i, 2)
c 	    j =  shared_list(i, 3)
c 	    velo(1,pt) = ex_shared(1,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j)
c  	    velo(2,pt) = ex_shared(2,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j)
c  	    velo(3,pt) = ex_shared(3,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j)
c1	continue
c	end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c An alternative version that averages the velocities instead of exchanging them
c	subroutine     get_others_V()
c	The  subroutine copies the velocities (calculated in the previous iteration
c	by other PEs) of the shared particals.
c
c	include 'COMMON/LENGTH.BLOCK'
c	include 'COMMON/COORD.BLOCK'
c	include 'COMMON/VELOC.BLOCK'
c	include 'COMMON/CONNECT.BLOCK'
c	include 'COMMON/UNITS.BLOCK'
c	include 'COMMON/SHAKE.BLOCK'
c	include 'COMMON/PARALLEL.BLOCK'
c	include 'COMMON/DEBUG.BLOCK'
c
c
c	integer i, pt,  pe, j
c        double precision  z,temp
c	do  1 i = 1, n_shared_total
c	    pt = shared_list(i, 1)
c	    pe = shared_list(i, 2)
c	    j =  shared_list(i, 3)
c	    velo(1,pt) = (velo(1,pt) + ex_shared(1,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j))/2.d0
c 	    velo(2,pt) = (velo(2,pt) + ex_shared(2,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j))/2.d0
c 	    velo(3,pt) = (velo(3,pt) + ex_shared(3,pe * mx_shk_shr_tot + my_pe *
c     *		                     mx_shk_shr + j))/2.d0
c1	continue
c	end
c---------------------------------------------------------------------
c	subroutine     put_my_V()
c	The  subroutine copies the velocities (calculated in this iteration)
c	of the shared particals.

c	include 'COMMON/LENGTH.BLOCK'
c	include 'COMMON/COORD.BLOCK'
c	include 'COMMON/VELOC.BLOCK'
c	include 'COMMON/CONNECT.BLOCK'
c	include 'COMMON/UNITS.BLOCK'
c	include 'COMMON/SHAKE.BLOCK'
c	include 'COMMON/PARALLEL.BLOCK'
c	include 'COMMON/DEBUG.BLOCK'


c	integer i, pt,  pe, j

c	do  1 i = 1, n_shared_total
c	    pt = shared_list(i, 1)
c	    pe = shared_list(i, 2)
c	    j =  shared_list(i, 3)
c	    ex_shared(1,my_pe * mx_shk_shr_tot + pe * mx_shk_shr + j)
c     *	             = velo(1,pt)
c	    ex_shared(2,my_pe * mx_shk_shr_tot + pe * mx_shk_shr + j)
c     *	             = velo(2,pt)
c	    ex_shared(3,my_pe * mx_shk_shr_tot + pe * mx_shk_shr + j)
c     *	             = velo(3,pt)
c1	continue
c	end
c---------------------------------------------------------------------
c	subroutine     update_v()
c	The  subroutine merges the velocities calculated by all the PEs

c	include 'COMMON/LENGTH.BLOCK'
c	include 'COMMON/COORD.BLOCK'
c	include 'COMMON/VELOC.BLOCK'
c	include 'COMMON/CONNECT.BLOCK'
c	include 'COMMON/UNITS.BLOCK'
c	include 'COMMON/SHAKE.BLOCK'
c	include 'COMMON/PARALLEL.BLOCK'
c	include 'COMMON/DEBUG.BLOCK'


c	integer i, pt,  pe, j

c	do 250 i = 1, n_shaked_prtcls
c	    j = shaked_prtcls(i)
c	    shak_velo(1,j) = velo(1,j)
c	    shak_velo(2,j) = velo(2,j)
c	    shak_velo(3,j) = velo(3,j)
c250	continue
c	call add_v()
c	do 260 i = 1, n_shared_total
c	    j =  shared_list(i, 1)
c	    velo(1,j) = velo(1,j) / 2.d0
c	    velo(2,j) = velo(2,j) / 2.d0
c	    velo(3,j) = velo(3,j) / 2.d0
c260	continue
c	end




