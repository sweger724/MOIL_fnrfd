        subroutine shakept(maxit)
        implicit none
c
c correct bond distances to statisfy shake constraints
c epsilon - allowed average error in distance
c maxit   - maximum number of iteration
c
c       double precision epsilon
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

      !write(6,*)"NNN",maxit,nshak

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
              if (rx .gt. 1d10 .or. ry .gt. 1d10 .or. rz .gt. 1d10) then
          write(6,*) " *** Shake diverged in ",iterate, "iteration! ***"
                stop
              end if
              
C                  get update factors
              Mg = ((rx*rx + ry*ry + rz*rz) - dissq(k)) /
     1        (2.d0*(rx*cooref(1,k) + ry*cooref(2,k) + rz*cooref(3,k)))

              fac1 = ovrlx*Mg*mass1(k)
              fac2 = ovrlx*Mg*mass2(k)
              !write(6,*)"FFF",k,dissq(k),
!     1         ((rx*rx + ry*ry + rz*rz) - dissq(k))

C                  update coordinates along original difference vector
c@              if (k.eq.istart) then
c@                      write(*,*)' ************* corrections '
c@                      write(*,*)' k iat1 iat2 ',k,iat1,iat2
c@                      write(*,*) ' Mg fac1 fac2 ',Mg,fac1,fac2
c@                      write(*,*)' corrections = '
c@                      write(*,*)' 1 = ',fac1*cooref(1,k)
c@                      write(*,*)' 2 = ',fac2*cooref(1,k)
c@                      write(*,*)' ************* corrections end'
c@              end if
              !write(6,*)"VVV",iat1,velo(1,iat1)
              velo(1,iat1) = velo(1,iat1) + fac1*cooref(1,k)
              velo(2,iat1) = velo(2,iat1) + fac1*cooref(2,k)
              velo(3,iat1) = velo(3,iat1) + fac1*cooref(3,k)

              velo(1,iat2) = velo(1,iat2) + fac2*cooref(1,k)
              velo(2,iat2) = velo(2,iat2) + fac2*cooref(2,k)
              velo(3,iat2) = velo(3,iat2) + fac2*cooref(3,k)
              !write(6,*)"UUU",iat1,velo(1,iat1)

100       continue

C          check for convergence of all distances and return if ok
          do 200 k=istart,nshak
              iat1 = ishak1(k)
              iat2 = ishak2(k)
              rx = cooref(1,k) + velo(1,iat1)-velo(1,iat2)
              ry = cooref(2,k) + velo(2,iat1)-velo(2,iat2)
              rz = cooref(3,k) + velo(3,iat1)-velo(3,iat2)
              r2   = rx*rx + ry*ry + rz*rz - dissq(k)
              !write(6,*)"RRR:",r2,velo(1,iat1),iat1
              if (dabs(r2) .gt. shak_diff(k)) then
                istart = k
                again = .true.
                go to 1000
              end if
200       continue
          if (again) then
             again = .false.
             istart= 1
             go to 1000
          end if
             return
 1000     continue

C
C If this point is reached, SHAKE did not converge! Do something about errors!
C
      level = 1
      write(6,*)"Shake failed on",istart,ishak1(istart),ishak2(istart)
      call alert('shakept',7,' Unconverged after maxit steps',30,level)
      return
      end
