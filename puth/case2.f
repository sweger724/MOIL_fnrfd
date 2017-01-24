	subroutine case2(ri,rj,rk,rl,dij,djk,dkl,cijk,cjkl,cpijkl)
	implicit none
c
c building the coordinate of atom(i) based on wcon entries and the coordinates
c of nearby atoms j,k, and l.
c dnm  - the equilbrium distance between atom n and m
c cjkl - cosine of theta(kjl)
c cpisjkl - cosine of the torsion angle ijkl
c
	double precision ri(3),rj(3),rk(3),rl(3)
	double precision dij,djk,dkl,cijk,cjkl,cpijkl

c Local
c la,lb,lc : linear coefficients expression ri in terms of other vectors
c vjk & vjl : vector of differences
c vprod : vector product
c
c rij = la*rjk + lb*rkl + lc*(rjk x rkl)

	double precision sijk,sjkl
	double precision la,lb,lc
	double precision vjk(3),vkl(3)
	double precision vprod(3)
	integer i
	character*5 name
	integer namel
	integer level

	name = 'case2'
	namel = 5

	if ( rj(1).gt.9998. .or.
     1    rk(1).gt.9998. .or. rl(1).gt.9998) then
		level = 3
		call alert(name,namel,' > 1 crd missing',16,level)
		write(*,*)'  rj rk rl  ',rj,rk,rl
		return
	end if

	sijk = sqrt(1-cijk*cijk)
	sjkl = sqrt(1-cjkl*cjkl)
	la = dij*sijk*cpijkl/(dkl*sjkl)
	lb = (-dij*cijk+la*dkl*cjkl)/djk
	lc = dij*dij-la*la*dkl*dkl+2*la*lb*djk*dkl
     1   *cjkl-lb*lb*djk*djk
	lc = lc/(djk*djk*dkl*dkl*(1-cjkl*cjkl))
        if (dabs(lc).lt.1.d-8) lc = 0.d0
c@        write(*,*)' la lb lc ',la,lb,lc
	lc = sqrt(lc)

c
c compute vector differences
c
	do i=1,3
	 vjk(i)=rj(i)-rk(i)
	 vkl(i)=rk(i)-rl(i)
	end do
c
c compute vector product rkl x rjk
c

	vprod(1)=vkl(2)*vjk(3)-vkl(3)*vjk(2)
	vprod(2)=vkl(3)*vjk(1)-vkl(1)*vjk(3)
	vprod(3)=vkl(1)*vjk(2)-vkl(2)*vjk(1)
c
c compute rij (placed temporarily in ri)
c
	do i=1,3
         ri(i) = la*vkl(i)+lb*vjk(i)+lc*vprod(i)
	end do
c
c compu ri
c
	do i=1,3
	 ri(i) = ri(i) + rj(i)
	end do

	return
	end
