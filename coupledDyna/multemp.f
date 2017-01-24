	subroutine multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1		curtemp,ntemp)
c
c calculate values of multiple temperatures: groups are defined by tpo
c vx vy vz - velocities
c natom    - number of particles (dimension of the vectors)
c inofrz   - number of particles that are not freezed
c nofreez  - vector pointer to the position of the particles
c
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/LENGTH.BLOCK'

	double precision velo(3,*)
	double precision ptms(*),curtemp(*)
	integer inofrz,ntemp,nofreez(*),tpo(*),tgroup(*)
c
c local
c
	integer i,l,n
	double precision factor
	data factor/1.98768d-3/

	do 1 n=1,ntemp
	 curtemp(n) = 0.d0
c@@
c	write(*,*)' n tgroup(n) ',n, tgroup(n)
C@@
1	continue
	do 2 i=1,inofrz
	  l = nofreez(i)
	  n = tpo(l)
	  curtemp(n) = curtemp(n) + ptms(l) * (velo(1,l)*velo(1,l)
     1		+ velo(2,l)*velo(2,l) + velo(3,l)*velo(3,l))
2	continue

	do 3 n=1,ntemp
	   curtemp(n)=curtemp(n)/(factor*tgroup(n))
c@@
c	write(*,*) ' curtemp(n) = ', curtemp(n)
c@@
3	continue
	return
	end

	subroutine multemp_prll(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1		curtemp,ntemp)
c
c calculate values of multiple temperatures: groups are defined by tpo
c vx vy vz - velocities
c natom    - number of particles (dimension of the vectors)
c inofrz   - number of particles that are not freezed
c nofreez  - vector pointer to the position of the particles
c
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/SHAKE.BLOCK'

	double precision velo(3,*)
	double precision ptms(*),curtemp(*)
	integer inofrz,ntemp,nofreez(*),tpo(*),tgroup(*)
c
c local
c
	integer l,n
	double precision factor
	data factor/1.98768d-3/

	do 1 n=1,ntemp
	 curtemp(n) = 0.d0
c@@
c	write(*,*)' n tgroup(n) ',n, tgroup(n)
C@@
1	continue
	do 2 l=pt_start,pt_end
c	  l = nofreez(i)
	  n = tpo(l)
	  curtemp(n) = curtemp(n) + ptms(l) * (velo(1,l)*velo(1,l)
     1		+ velo(2,l)*velo(2,l) + velo(3,l)*velo(3,l))
2	continue

	do 3 n=1,ntemp
	   call reduce_1(curtemp(n))
	   curtemp(n)=curtemp(n)/(factor*tgroup(n))
c@@
c	write(*,*) ' curtemp(n) = ', curtemp(n)
c@@
3	continue
	return
	end

