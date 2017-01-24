	subroutine case1(ri,rj,rk,rl,dij,djk,djl,ckjl,cijk,cijl)
	implicit none
c
c building the coordinate of atom(i) based on wcon entries and the coordinates
c of nearby atoms j,k, and l.
c dnm  - the equilbrium distance between atom n and m
c ckjl - cosine of theta(kjl)
c
	double precision ri(3),rj(3),rk(3),rl(3)
	double precision dij,djk,djl,ckjl,cijk,cijl

c Local
c la,lb,lc : linear coefficients expression ri in terms of other vectors
c vjk & vjl : vector of differences
c vprod : vector product
c
c rij = la*rjk + lb*rjl + lc*(rjk x rjl)

	double precision la,lb,lc
	double precision vjk(3),vjl(3)
	double precision vprod(3)
	integer i,j,k,l
	integer level

	character*5 name
	integer namel

	name = 'case1'
	namel = 5

	if ( rj(1).gt.9998. .or.
     1    rk(1).gt.9998. .or. rl(1).gt.9998) then
		level = 3
		write(*,*)'  rj rk rl  ',rj,rk,rl
		call alert(name,namel,' > 1 crd missed',15,level)
		return
	end if

	lb = dij*(cijk*ckjl-cijl)/(djl*(1.-ckjl*ckjl))
	la = (-dij*cijk-lb*djl*ckjl)/djk
	lc = (dij*dij-la*la*djk*djk-
     1      lb*lb*djl*djl-2*la*lb*djk*djl*ckjl)/
     2      (djk*djk*djl*djl*(1-ckjl*ckjl))
        if (lc.lt.0) then
         lc= 0.d0
        else
	 lc = sqrt(lc)
        end if

c
c compute vector differences
c
	do i=1,3
	 vjk(i)=rj(i)-rk(i)
	 vjl(i)=rj(i)-rl(i)
	end do
c
c compute vector product rjk x rjl
c
	vprod(1)=vjk(2)*vjl(3)-vjk(3)*vjl(2)
	vprod(2)=-vjk(1)*vjl(3)+vjk(3)*vjl(1)
	vprod(3)=vjk(1)*vjl(2)-vjk(2)*vjl(1)

c
c compute rij (placed temporarily in ri)
c
	do i=1,3
         ri(i) = la*vjk(i)+lb*vjl(i)+lc*vprod(i)
	end do
c
c compu ri
c
	do i=1,3
	 ri(i) = ri(i) + rj(i)
	end do

	return
	end
