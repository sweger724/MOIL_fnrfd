	subroutine special

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/CONSPECL1.BLOCK'
	include 'COMMON/CONSPECL2.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/ENERGY.BLOCK'

        double precision rr(maxmorsb)
        double precision fe(maxmorsb),dr(maxmorsb)
        double precision f2(maxmorsb),temprx(maxmorsb)
	double precision tempry(maxmorsb)
	double precision temprz(maxmorsb)
	double precision tmp
	integer i,j,n
       
c@
        double precision dx3(maxpt),dy3(maxpt),dz3(maxpt)
ccccccccccccccccccccccccccccccccccccccccccc
        poitype1(0)= 0
c poitype1(n-1)+1 is the first particle of n-th special unit which 
c has morse bond
c poitype1(n )is the last partical of n-th special unit which has 
c morse bond
	do 3 i=1,poitype1(nmb)
	 dx1(i)=0.d0 
	 dy1(i)=0.d0
	 dz1(i)=0.d0
	 sdx(i)=0.d0
	 sdy(i)=0.d0
	 sdz(i)=0.d0
3       continue
ccccccccccccccccccccccccccccccccccccccccccc
       do 777 n = 1,nmb
	if (emyes(n)) then
         temprx(n) = coor(1,imb1(n)) - coor(1,imb2(n))
         tempry(n) = coor(2,imb1(n)) - coor(2,imb2(n))
         temprz(n) = coor(3,imb1(n)) - coor(3,imb2(n))
	 rr(n) = temprx(n)*temprx(n) + tempry(n)*tempry(n) 
     1         + temprz(n)*temprz(n)
	 dist(n) = DSQRT(rr(n))
	 dr(n) = dist(n) - rcut(n)
	 fe(n) = -lamda(n)*dr(n)
	 f(n) = 1.d0 + DEXP(fe(n))
	 f(n) = 1.d0/f(n)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if (debug) then
          print*,'n dist(n)=',n,dist(n) 
          print*,'f(n)=',f(n)
        end if

	 f2(n) = f(n)*f(n)
	 tmpf(n) = f2(n)*lamda(n)*DEXP(fe(n))
	 tmp = 1.d0/dist(n)
	 tempdfx(n) = tmpf(n)*temprx(n)*tmp
	 tempdfy(n) = tmpf(n)*tempry(n)*tmp
	 tempdfz(n) = tmpf(n)*temprz(n)*tmp
c
       if((ebyes .and. snb1(n).gt.0).or.(ebyes .and. snb2(n) .gt.0))
     1	call ebond_specl(n)
	e_bond = e_bond - se_bond1 + f(n)*se_bond1 +(1-f(n))*se_bond2
	e_total= e_total- se_bond1 + f(n)*se_bond1 +(1-f(n))*se_bond2
	if ((ethyes .and. snang1(n).gt.0) .or. (ethyes .and. snang2(n)
     1	.gt.0))
     1	call etheta_specl(n)
	e_theta = e_theta - se_theta1 + f(n)*se_theta1 + (1-f(n))*
     1	se_theta2
	e_total = e_total - se_theta1 + f(n)*se_theta1 + (1-f(n))*
     1	se_theta2
c ---------------------------------------------
c        if (ctrue) then
c	 call cdie_spcl(n)
c the coding here is not great. The problem is that the special component
c is computed first by the ordinary energy and then re-calculated in 
c lower_excited part but with a different defintion of the cutoff. Common
c sense suggests to calculate it ONLY in the lower_upper part, and perhaps to
c further modify the cutoff scheme there.
c At present it is left "as is".
c
c        write(*,*)' el1 el2 vdw1 dvw2 ',se_el1,se_el2,se_vdw1,se_vdw2
c        write(*,*)' Before e_el = ',e_el
c	e_el    =  e_el - se_el1 + f(n)*se_el1 + (1-f(n))*se_el2 
c        write(*,*) 'After e_el ' ,e_el
c        write(*,*)' Before e_vdw = ',e_vdw
c	e_vdw   =  e_vdw - se_vdw1 + f(n)*se_vdw1 + (1-f(n))*se_vdw2 
c        write(*,*)' After e_vdw = ',e_vdw
c	e_total = e_total - se_vdw1 -se_el1 + f(n)*(se_vdw1+se_el1) 
c     1            + (1-f(n))*(se_el2 + se_vdw2)
c	end if
      
c ------------------------------------------         
	  
        do 20 j=poitype1(n-1)+1,poitype1(n)
	 i=newpoit(j)
	 dpot(1,i) = dpot(1,i) - dx1(j) + sdx(j)
	 dpot(2,i) = dpot(2,i) - dy1(j) + sdy(j)
	 dpot(3,i) = dpot(3,i) - dz1(j) + sdz(j) 
20      continue
ccccccccccccccccccccccccccccccccccccccccccccccccc
       end if
777   continue
cccccccccccccccccccccccccccccccccccccccccccccccccc
c ************************************************************
	return
	end
