	subroutine shakecons(epsilon,maxit,iac1,iac2)
c
c correct bond distances to statisfy shake constraints
c epsilon - allowed average error in distance
c maxit   - maximum number of iteration
c
	double precision epsilon
	integer maxit,iac1,iac2

	include 'COMMON/LENGTH.BLOCK'
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
C            dissq                   squares of shake distances
C            ptms                    atomic masses
C            x     , y     , z       coordinates after unconstrained move
C            xref ,  yref ,  zref    coordinates before unconstrained move
C            epsilon                 convergence criterion
C            maxit                 maximum number of iterations

C  OUTPUT -- x, y, z  updated coordinates


      double precision rx, ry, rz
      double precision Mg, fac1, fac2, diff
      integer iat1, iat2, k, iterate, level
      logical ok

C  check each constraint at most loopmax times
      do 1000 iterate=1,maxit

C          adjust each constraint in order
          do 100 k=1,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)

	      if((iat1.eq.iac1).and.(iat2.eq.iac2)) then
C                compute difference vectors
              rx = coor(1,iat1) - coor(1,iat2)
              ry = coor(2,iat1) - coor(2,iat2)
              rz = coor(3,iat1) - coor(3,iat2)

C                  get update factors
              Mg = (dissq(k) - (rx*rx + ry*ry + rz*rz))
     1          /(2.*(rx*cooref(1,k) + ry*cooref(2,k) + rz*cooref(3,k)))

C                 If both atoms have zero mass, this if check is necessary.

              if (ptms(iat2) .eq. ptms(iat1)) then
                  fac1 = 0.5 * Mg
                  fac2 = 0.5 * Mg
              else
                  fac1 = Mg*ptms(iat2)/(ptms(iat2)+ ptms(iat1))
                  fac2 = Mg*ptms(iat1)/(ptms(iat2)+ ptms(iat1))
              endif

C                  update coordinates along original difference vector
              coor(1,iat1) = coor(1,iat1) + fac1 * cooref(1,k)
              coor(2,iat1) = coor(2,iat1) + fac1 * cooref(2,k)
              coor(3,iat1) = coor(3,iat1) + fac1 * cooref(3,k)
              coor(1,iat2) = coor(1,iat2) - fac2 * cooref(1,k)
              coor(2,iat2) = coor(2,iat2) - fac2 * cooref(2,k)
              coor(3,iat2) = coor(3,iat2) - fac2 * cooref(3,k)
	      end if
100       continue

C          check for convergence of all distances and return if ok
          ok = .true.
          do 200 k=1,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)

	      if((iat1.eq.iac1).and.(iat2.eq.iac2)) then

              rx = coor(1,iat1)-coor(1,iat2)
              ry = coor(2,iat1)-coor(2,iat2)
              rz = coor(3,iat1)-coor(3,iat2)
	      diff = (1.d0 - (rx*rx+ry*ry+rz*rz)/dissq(k))
              ok = ok .and. dabs(diff).lt.epsilon
	      end if
200       continue
	  if (ok) return
1000  continue
165       format(2I8,E15.5)

C
C If this point is reached, SHAKE did not converge! Do something about errors!
C
      level = 1
      call alert('shakept',7,' Unconverged after maxit steps',30,level)
      return
      end
