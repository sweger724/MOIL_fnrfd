	subroutine shakevl(epsilon,maxit)
c
c correct bond distances to statisfy shake constraints
c epsilon - allowed average error in distance
c maxit   - maximum number of iteration
c      
	double precision epsilon
	integer maxit

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'

C  Uses SHAKE algorithm to update coordinates and maintain 
C  distance contstraints.  This implementation handles atoms with
C  zero mass gracefully and correctly. 

C  INPUT --  nshak:                  number of shake constraints
C            ishak1, ishak2          pairs of atoms to apply shake constraints
C            ptms                    atomic masses
C            vx     , vy     , vz    velocities after unconstrained move
C            xref ,  yref ,  zref    coordinates after unconstrained move
C            epsilon                 convergence criterion
C            maxit                 maximum number of iterations

C  OUTPUT -- vx, vy, vz  updated velocities


      double precision vrx, vry, vrz
      double precision Mg, fac1, fac2, diff
c Note that the dimension of fac[x,y,z] must be the larger
c number between maxpt/maxshak
c
      double precision facx(maxshak),facy(maxshak),facz(maxshak)
      double precision mass1(maxshak),mass2(maxshak),fac(maxshak)
      double precision ovrlx
      integer iat1, iat2, k, iterate, level, maxiter, istart
      logical first,again
      parameter (ovrlx=1.3d0)
      data first /.true./

      save first
      save mass1,mass2,fac

	if (first) then
	 do 1 k=1,nshak
          iat1 = ishak1(k)
          iat2 = ishak2(k)
	  mass1(k) = -ptms(iat2)/(ptms(iat1)+ptms(iat2))
	  mass2(k) = 1.d0 + mass1(k)
	  fac(k)   = 1.d0/dissq(k)
1	 continue
	 first = .false.
	end if
	maxiter = maxit
	istart  = 1
	again   = .false.
      
        do 2 k=1,nshak
                facx(k) = cooref(1,k)*fac(k)
                facy(k) = cooref(2,k)*fac(k)
                facz(k) = cooref(3,k)*fac(k)
2       continue

C  check each constraint at most loopmax times

      do 1000 iterate=1,maxiter

C          adjust each constraint in order
          do 100 k=istart,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)
              vrx = velo(1,iat1) - velo(1,iat2)
              vry = velo(2,iat1) - velo(2,iat2)
              vrz = velo(3,iat1) - velo(3,iat2)
              
C                  get update factors
              Mg =  (vrx*facx(k) + vry*facy(k) + vrz*facz(k))

	      fac1 = ovrlx*mass1(k)*Mg
	      fac2 = ovrlx*mass2(k)*Mg

C                  update coordinates along original difference vector
               velo(1,iat1) = velo(1,iat1) +  fac1 * cooref(1,k)
               velo(2,iat1) = velo(2,iat1) +  fac1 * cooref(2,k)
               velo(3,iat1) = velo(3,iat1) +  fac1 * cooref(3,k)
               velo(1,iat2) = velo(1,iat2) +  fac2 * cooref(1,k)
               velo(2,iat2) = velo(2,iat2) +  fac2 * cooref(2,k)
               velo(3,iat2) = velo(3,iat2) +  fac2 * cooref(3,k)
100       continue
c chen
c	if (prll_on_off) then
c	    call put_my_V()
c	    call update_ex_shared()
c	    call get_others_V()
c	end if

C          check for convergence of all distances and return if ok
          do 200 k=1,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)
              vrx = velo(1,iat1) - velo(1,iat2)
              vry = velo(2,iat1) - velo(2,iat2)
              vrz = velo(3,iat1) - velo(3,iat2)
              diff = (vrx*facx(k)+vry*facy(k)+vrz*facz(k))
	      if (dabs(diff).gt.epsilon) then
c chen
c		if (prll_on_off) call shak_unhappy()
		istart = k
		again = .true.
		go to 1000
	      end if
200       continue
	  if (again) then
		again  = .false.
		istart = 1
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
c	write(stdo,*)'shakpvel converged after ',iterate, ' iterations'
c		    return
c	      else
c		    go to 1000
c	    end if
c	end if
          return
1000	  continue
	  level = 1
      call alert('shakevl',7,' Velocity shake not converged',29,level)
	  return
	  end
