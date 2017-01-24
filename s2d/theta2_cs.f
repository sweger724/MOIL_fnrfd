      subroutine theta2()
c---This is theta2.cos written by  M. Kara-Ivanov
c--- On the basis of Ron's program!
c-- it will help to get rid of singularities sintheta approaching 0
c---in the case if the angles are close to 180 deg !
c---Idea is to use E=k(cos - coseq)*(cos-coseq): i.e. cos instead of angles!
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SCNDRV.BLOCK'
      include 'COMMON/UNITS.BLOCK'

c calculate second derivative of angle energy. The second
c derivatives are stored in a condensed (but not MAXIMALLY
c CONDENSED) way. Diagonal elements (in atom index) are added
c to the diag(6*maxpt) vector of the second derivative matrix.
c The off-diagonal elements are stored as derivatives of individual
c angles. I.e. for each angle we associate an off-diagonal block
c for the 3 pairs i,j i,k j,k possible in the angle defined by the
c three particles i-j-k. For each such pair there are nine elements
c for all x/y/z combinations. Length of off-diagonal vector is therefore
c 27*maxangl. It is compact since none of the elements is always zero. It is
c not fully compressed since particles may appear more than once, since the
c storage is done according to the angle index and not according to the
c particle index. In the present formulation of the energy it is
c quite convenient.
c

	integer ith,i,j,k,ii,jj,kk,iith
	double precision qxij,qyij,qzij,qxkj,qykj,qzkj
	double precision rijrkj,rij2,rkj2,rij,rkj,r,r2
	double precision c1,c2,delijx,delijy,delijz
	double precision delkjx,delkjy,delkjz,cost
	double precision scalar1,scalar2,aix,aiy,aiz
	double precision ajx,ajy,ajz,akx,aky,akz,a11
	double precision gamijx,gamijy,gamijz,gamkjx,gamkjy,gamkjz

	character*6 name
	integer namel, level
	data name/'theta2'/
	data namel,level/6,1/




c initialize loop over angles in ichunk chunks
      do 100 ith=1,nangl

c maintain the pointer to the "real" vector

      		i=iangl1(ith)
      		j=iangl2(ith)
      		k=iangl3(ith)

      		qxij   = coor(1,i)-coor(1,j)
      		qyij   = coor(2,i)-coor(2,j)
      		qzij   = coor(3,i)-coor(3,j)
      		qxkj   = coor(1,k)-coor(1,j)
      		qykj   = coor(2,k)-coor(2,j)
      		qzkj   = coor(3,k)-coor(3,j)

		rijrkj = (qxij*qxkj + qyij*qykj + qzij*qzkj)
      		rij2   = 1.d0/(qxij*qxij+qyij*qyij+qzij*qzij)
      		rkj2   = 1.d0/(qxkj*qxkj+qykj*qykj+qzkj*qzkj)
      		rij    = dsqrt(rij2)
      		rkj    = dsqrt(rkj2)
		r      = (rij*rkj)
		r2     = r*r
		c1     = rijrkj*rij2
		c2     = rijrkj*rkj2
		gamijx = qxij*rij2
		gamijy = qyij*rij2
		gamijz = qzij*rij2
		gamkjx = qxkj*rkj2
		gamkjy = qykj*rkj2
		gamkjz = qzkj*rkj2

		delijx = r*gamijx
		delijy = r*gamijy
		delijz = r*gamijz
		delkjx = r*gamkjx
		delkjy = r*gamkjy
		delkjz = r*gamkjz

		cost   = rijrkj*r
		if (dabs(cost).ge.1.d0) then
		  level = 0
		  call alert(name,namel,' Linear angle ',14,level)
		end if

c---Now: e=k(cst-coseq)*(cst-coseq) This is angle energy!
c---Derivative: de/dx=2k(cst-coseq)*dcst/dx
c---Second Derivative: d2e/dxi*dxj=scalar1*dcst/dxi*dcst/dxj  +
c--  +  scalar2*d2cst/dxidxj
		scalar1 = 2.d0*kangl(ith)
		scalar2 = 2.d0*kangl(ith)*(cost-dcos(angleq(ith)))

		ii        = 6*(i-1)+1
		jj        = 6*(j-1)+1
		kk        = 6*(k-1)+1
		iith      = 27*(ith-1)+1

		aix       = r*(qxkj - c1*qxij)
		aiy       = r*(qykj - c1*qyij)
		aiz       = r*(qzkj - c1*qzij)

		ajx       = r*(-(qxij+qxkj) + c1*qxij + c2*qxkj)
		ajy       = r*(-(qyij+qykj) + c1*qyij + c2*qykj)
		ajz       = r*(-(qzij+qzkj) + c1*qzij + c2*qzkj)

		akx       = r*(qxij - c2*qxkj)
		aky       = r*(qyij - c2*qykj)
		akz       = r*(qzij - c2*qzkj)


C diagonal elements: first comes i then j and then k
C
C here comes the six elements of i,i

		a11       = scalar1*aix*aix - scalar2*(2.d0*delijx*qxkj 
     1		 - 3.d0*cost*gamijx*gamijx)
c xx
		diag(ii)  = diag(ii) + a11 - scalar2*cost*rij2
c yx
		diag(ii+1)= diag(ii+1) + scalar1*aiy*aix - scalar2*(
     1		 delijy*qxkj + delijx*qykj - 3.d0*cost*gamijx*gamijy)
c zx
		diag(ii+2)= diag(ii+2) + scalar1*aiz*aix - scalar2*(
     1		 delijz*qxkj + delijx*qzkj - 3.d0*cost*gamijx*gamijz)
		a11       = scalar1*aiy*aiy - scalar2*(2.d0*delijy*qykj
     1		 - 3.d0*cost*gamijy*gamijy)
c yy
		diag(ii+3)= diag(ii+3) + a11 - scalar2*cost*rij2
c zy
		diag(ii+4)= diag(ii+4) + scalar1*aiz*aiy - scalar2*(
     1		 delijz*qykj + delijy*qzkj - 3.d0*cost*gamijz*gamijy)
		a11       = scalar1*aiz*aiz - scalar2*(2.d0*delijz*qzkj 
     1		 - 3.d0*cost*gamijz*gamijz)
c zz
		diag(ii+5)= diag(ii+5) + a11 - scalar2*cost*rij2
C
C here come the six elements of j,j

		a11       = scalar1*ajx*ajx +
     1		 scalar2*(-2.d0*(delkjx+delijx)*(qxkj+qxij)+2.d0*r
     2		 + cost*(gamijx+gamkjx)**2 + 2.d0*cost*(gamijx*gamijx
     3		 +gamkjx*gamkjx))
c xx
		diag(jj)  = diag(jj) + a11 - scalar2*cost*(rij2+rkj2)
c yx
		diag(jj+1)= diag(jj+1) + scalar1*ajy*ajx + scalar2*(
     1		 -(delijy+delkjy)*(qxkj+qxij)-(delijx+delkjx)*
     2		 (qykj+qyij)+cost*(gamkjy+gamijy)*(gamkjx+gamijx)
     3		 +2.d0*cost*(gamijy*gamijx+gamkjy*gamkjx))
c zx
		diag(jj+2)= diag(jj+2) + scalar1*ajz*ajx + scalar2*(
     1		 -(delijz+delkjz)*(qxkj+qxij)-(delkjx+delijx)*
     2		 (qzkj+qzij)+cost*(gamkjz+gamijz)*(gamkjx+gamijx)
     3		 +2.d0*cost*(gamijz*gamijx+gamkjz*gamkjx))
c yy
		a11       = scalar1*ajy*ajy +
     1		 scalar2*(-2.d0*(delkjy+delijy)*(qykj+qyij)+2.d0*r
     2		 + cost*(gamijy+gamkjy)**2 + 2.d0*cost*(gamijy*gamijy
     3		 +gamkjy*gamkjy))
		diag(jj+3)= diag(jj+3) + a11 - scalar2*cost*(rij2+rkj2)
c yz
		diag(jj+4)= diag(jj+4) + scalar1*ajy*ajz + scalar2*(
     1		 -(delijy+delkjy)*(qzkj+qzij)-(delijz+delkjz)*
     2		 (qykj+qyij)+cost*(gamkjy+gamijy)*(gamkjz+gamijz)
     3		 +2.d0*cost*(gamijy*gamijz+gamkjy*gamkjz))
c zz
		a11       = scalar1*ajz*ajz +
     1		 scalar2*(-2.d0*(delkjz+delijz)*(qzkj+qzij)+2.d0*r
     2		 + cost*(gamijz+gamkjz)**2 + 2.d0*cost*(gamijz*gamijz
     3		 +gamkjz*gamkjz))
		diag(jj+5)= diag(jj+5) + a11 - scalar2*cost*(rij2+rkj2)

C Here come the k,k elements

c xx
		diag(kk)  = diag(kk) + scalar1*akx*akx + scalar2*(
     1		 -2.d0*delkjx*qxij+3.d0*cost*gamkjx*gamkjx-cost*rkj2)
c yx
		diag(kk+1)= diag(kk+1) + scalar1*akx*aky + scalar2*(
     1		 -delkjy*qxij-delkjx*qyij+3.d0*cost*gamkjx*gamkjy)
c zx
		diag(kk+2)= diag(kk+2) + scalar1*akx*akz + scalar2*(
     1		 -delkjz*qxij-delkjx*qzij+3.d0*cost*gamkjx*gamkjz)
c yy
		diag(kk+3)= diag(kk+3) + scalar1*aky*aky + scalar2*(
     1		 -2.d0*delkjy*qyij+3.d0*cost*gamkjy*gamkjy-cost*rkj2)
c yz
		diag(kk+4)= diag(kk+4) + scalar1*akz*aky + scalar2*(
     1		 -delkjy*qzij-delkjz*qyij+3.d0*cost*gamkjz*gamkjy)
c zz
		diag(kk+5)= diag(kk+5) + scalar1*akz*akz + scalar2*(
     1		 -2.d0*delkjz*qzij+3.d0*cost*gamkjz*gamkjz-cost*rkj2)

C Off diagonal elements: there are nine of them for each pair
C i,j i,k j,k overall 3*9=27 elements for each angle. It is
C possible that the same atom will have contribution from other angles
C However contributions here are calculated from
C individual angles. Only in final squeezing to atomic based matrix
C simple reference to particle indices can be made.
C

C pair i,j

c xi,xj (==xj,xi)
		d2theta(iith) = scalar1*aix*ajx + scalar2*(
     1		 delijx*(qxkj+qxij)-r+(delijx+delkjx)*qxkj -
     2		 cost*(gamijx+gamkjx)*gamijx-2.d0*cost*gamijx*gamijx
     3		 +cost*rij2)
c xi,yj
		d2theta(iith+1) = scalar1*aix*ajy + scalar2*(
     1		 delijx*(qykj+qyij)+(delijy+delkjy)*qxkj -cost*
     2		 (gamijy+gamkjy)*gamijx-2.d0*cost*gamijx*gamijy)
c xi,zj
		d2theta(iith+2) = scalar1*aix*ajz + scalar2*(
     1		 delijx*(qzkj+qzij)+(delijz+delkjz)*qxkj -cost*
     2		 (gamijz+gamkjz)*gamijx-2.d0*cost*gamijx*gamijz)
c yi,xj
		d2theta(iith+3) = scalar1*aiy*ajx + scalar2*(
     1		 delijy*(qxkj+qxij)+(delijx+delkjx)*qykj -cost*
     2		 (gamijx+gamkjx)*gamijy-2.d0*cost*gamijy*gamijx)
c yi,yj
		d2theta(iith+4) = scalar1*aiy*ajy + scalar2*(
     1		 delijy*(qykj+qyij)-r+(delijy+delkjy)*qykj -
     2		 cost*(gamijy+gamkjy)*gamijy-2.d0*cost*gamijy*gamijy
     3		 +cost*rij2)
c yi,zj
		d2theta(iith+5) = scalar1*aiy*ajz + scalar2*(
     1		 delijy*(qzkj+qzij)+(delijz+delkjz)*qykj -cost*
     2		 (gamijz+gamkjz)*gamijy-2.d0*cost*gamijy*gamijz)
c zi,xj
		d2theta(iith+6) = scalar1*aiz*ajx + scalar2*(
     1		 delijz*(qxkj+qxij)+(delijx+delkjx)*qzkj -cost*
     2		 (gamijx+gamkjx)*gamijz-2.d0*cost*gamijz*gamijx)
c zi,yj
		d2theta(iith+7) = scalar1*aiz*ajy + scalar2*(
     1		 delijz*(qykj+qyij)+(delijy+delkjy)*qzkj -cost*
     2		 (gamijy+gamkjy)*gamijz-2.d0*cost*gamijz*gamijy)
c zi,zj
		d2theta(iith+8) = scalar1*aiz*ajz + scalar2*(
     1		 delijz*(qzkj+qzij)-r+(delijz+delkjz)*qzkj -
     2		 cost*(gamijz+gamkjz)*gamijz-2.d0*cost*gamijz*gamijz
     3		 +cost*rij2)

C pair i,k

c xi,xk
		d2theta(iith+9) = scalar1*aix*akx + scalar2*(
     1		 -delijx*qxij + r - delkjx*qxkj + cost*gamijx*gamkjx)
c xi,yk
		d2theta(iith+10)= scalar1*aix*aky + scalar2*(
     1		 -delijx*qyij - delkjy*qxkj + cost*gamijx*gamkjy)
c xi,zk
		d2theta(iith+11)= scalar1*aix*akz + scalar2*(
     1		 -delijx*qzij - delkjz*qxkj + cost*gamijx*gamkjz)
c yi,xk
		d2theta(iith+12) = scalar1*aiy*akx + scalar2*(
     1		 -delijy*qxij - delkjx*qykj + cost*gamijy*gamkjx)
c yi,yk
		d2theta(iith+13)= scalar1*aiy*aky + scalar2*(
     1		 -delijy*qyij + r - delkjy*qykj + cost*gamijy*gamkjy)
c yi,zk
		d2theta(iith+14)= scalar1*aiy*akz + scalar2*(
     1		 -delijy*qzij - delkjz*qykj + cost*gamijy*gamkjz)
c zi,xk
		d2theta(iith+15) = scalar1*aiz*akx + scalar2*(
     1		 -delijz*qxij - delkjx*qzkj + cost*gamijz*gamkjx)
c zi,yk
		d2theta(iith+16)= scalar1*aiz*aky + scalar2*(
     1		 -delijz*qyij - delkjy*qzkj + cost*gamijz*gamkjy)
c zi,zk
		d2theta(iith+17)= scalar1*aiz*akz + scalar2*(
     1		 -delijz*qzij +r - delkjz*qzkj + cost*gamijz*gamkjz)

C Pair j,k

c xj,xk 
		d2theta(iith+18) = scalar1*ajx*akx + scalar2*(
     1		 delkjx*(qxkj+qxij)-r+(delijx+delkjx)*qxij -
     2		 cost*(gamijx+gamkjx)*gamkjx-2.d0*cost*gamkjx*gamkjx
     3		 +cost*rkj2)
c xj,yk
		d2theta(iith+19) = scalar1*ajx*aky + scalar2*(
     1		 delkjy*(qxkj+qxij)+(delijx+delkjx)*qyij -cost*
     2		 (gamijx+gamkjx)*gamkjy-2.d0*cost*gamkjx*gamkjy)
c xj,zk
		d2theta(iith+20) = scalar1*ajx*akz + scalar2*(
     1		 delkjz*(qxkj+qxij)+(delijx+delkjx)*qzij -cost*
     2		 (gamijx+gamkjx)*gamkjz-2.d0*cost*gamkjx*gamkjz)

c yj,xk 
		d2theta(iith+21) = scalar1*ajy*akx + scalar2*(
     1		 delkjx*(qykj+qyij)+(delijy+delkjy)*qxij -cost*
     2		 (gamijy+gamkjy)*gamkjx-2.d0*cost*gamkjy*gamkjx)
c yj,yk
		d2theta(iith+22) = scalar1*ajy*aky + scalar2*(
     1		 delkjy*(qykj+qyij)-r+(delijy+delkjy)*qyij -cost*
     2		 (gamijy+gamkjy)*gamkjy-2.d0*cost*gamkjy*gamkjy
     3		 +cost*rkj2)
c yj,zk
		d2theta(iith+23) = scalar1*ajy*akz + scalar2*(
     1		 delkjz*(qykj+qyij)+(delijy+delkjy)*qzij -cost*
     2		 (gamijy+gamkjy)*gamkjz-2.d0*cost*gamkjy*gamkjz)

c zj,xk 
		d2theta(iith+24) = scalar1*ajz*akx + scalar2*(
     1		 delkjx*(qzkj+qzij)+(delijz+delkjz)*qxij -cost*
     2		 (gamijz+gamkjz)*gamkjx-2.d0*cost*gamkjz*gamkjx)
c jj,yk
		d2theta(iith+25) = scalar1*ajz*aky + scalar2*(
     1		 delkjy*(qzkj+qzij)+(delijz+delkjz)*qyij -cost*
     2		 (gamijz+gamkjz)*gamkjy-2.d0*cost*gamkjz*gamkjy)
c zj,zk
		d2theta(iith+26) = scalar1*ajz*akz + scalar2*(
     1		 delkjz*(qzkj+qzij)-r+(delijz+delkjz)*qzij -cost*
     2		 (gamijz+gamkjz)*gamkjz-2.d0*cost*gamkjz*gamkjz
     3		 +cost*rkj2)

100	continue

      return	

      end
