	double precision function projf(dpot,g,nselec,pointr)

	implicit none
	integer i, j, k, nselec, pointr(*)
	double precision dpot(3,*), g(3,*)

c       projf is the scalar product between the force (-dpot)
c	and the normalized reaction path gradient (g, or "unitgrad")

	projf = 0.d0
	do 1 j=1,nselec
	   i = pointr(j)
	   do 2 k = 1,3
	      projf = projf + dpot(k,i) * g(k,j)
 2	continue
 1	continue

	projf = -projf

	return 
	end 
