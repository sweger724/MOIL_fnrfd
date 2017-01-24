	subroutine nbondm_spcl()
      	include 'COMMON/LENGTH.BLOCK'
      	include 'COMMON/CONSPECL1.BLOCK'
      	include 'COMMON/CONSPECL2.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	integer nlist1,nlist2

C-------local variables--------------------------
	integer ki,i,j
	double precision rxyz,r2
	integer pointm1(maxsmon+1)
	integer pointm2(maxsmon+1)
	integer nblistm1(sihugemo)
	integer nblistm2(sihugemo)

	if (.not. (evdyes.or.eelyes)) return

c  list of nonbond interaction when ligand does not bond
c ki = index pointing to first atom in monomer
	ki = 1
c nlist = index for monomer non-bonding list
	nlist1 = 0
c loop over first monomer
	if (debug) then
	 write(stdo,*)' in nbondm npt =',npt
	 write(stdo,1000)(poipt(i),i=1,totmon)
1000	 format(10(1x,i6))
	end if
        do 2 i=1,stmon1
         mono_cent(1,i) = 0.d0
         mono_cent(2,i) = 0.d0
         mono_cent(3,i) = 0.d0
         do 1 j=ki,spoipt1(i)
          mono_cent(1,i) = mono_cent(1,i) + coor(1,newpoit(j))
          mono_cent(2,i) = mono_cent(2,i) + coor(2,newpoit(j))
          mono_cent(3,i) = mono_cent(3,i) + coor(3,newpoit(j))
1        continue
         r2 = 1.d0/dfloat(spoipt1(i)-ki+1)
         mono_cent(1,i)  = mono_cent(1,i)*r2
         mono_cent(2,i)  = mono_cent(2,i)*r2
         mono_cent(3,i)  = mono_cent(3,i)*r2
         ki = spoipt1(i) + 1
2       continue

	
		
	do 100 i=1,stmon1 - 1 
c set pointers, include self in list
		nlist1 = nlist1+1
		pointm1(i)=nlist1
		nblistm1(nlist1) = i


c loop over second monomer
		do 200 j=i+1,stmon1
			rxyz = mono_cent(1,i) - mono_cent(1,j)
			r2 = rxyz*rxyz
			if (r2.gt.cutebig2) go to 200
			rxyz = mono_cent(2,i) - mono_cent(2,j)
			r2 = r2 + rxyz*rxyz
			if (r2.gt.cutebig2) go to 200
			rxyz = mono_cent(3,i) - mono_cent(3,i)
			r2 = r2 + rxyz*rxyz

c check distances^2. If smaller than cutebig2, include pair
c i,j in list, update pointers, go for next j.
			if(r2.le.cutebig2) then
				nlist1=nlist1+1
				nblistm1(nlist1)=j
			end if

200		continue
		ki=spoipt1(i)+1

100	continue

c set list for last monomer (only self)
	nlist1 = nlist1 + 1
	pointm1(stmon1)=nlist1
C the next sentence is just to insure the proper existence of
C pointm1(totmon)+1
	pointm1(stmon1+1)=nlist1 +1
	nblistm1(nlist1)=stmon1

c---------------------------------------------------------O

c  list of nonbond interaction when ligand bonded to MB

	ki = 1
        do 22 i=1,stmon1
         mono_cent(1,i) = 0.d0
         mono_cent(2,i) = 0.d0
         mono_cent(3,i) = 0.d0
         do 11 j=ki,spoipt2(i)
          mono_cent(1,i) = mono_cent(1,i) + coor(1,newpoit(j))
          mono_cent(2,i) = mono_cent(2,i) + coor(2,newpoit(j))
          mono_cent(3,i) = mono_cent(3,i) + coor(3,newpoit(j))
11        continue
         r2 = 1.d0/dfloat(spoipt2(i)-ki+1)
         mono_cent(1,i)  = mono_cent(1,i)*r2
         mono_cent(2,i)  = mono_cent(2,i)*r2
         mono_cent(3,i)  = mono_cent(3,i)*r2
         ki = spoipt2(i) + 1
22       continue

c nlist = index for monomer non-bonding list
	nlist2 = 0
c loop over first monomer
	if (debug) then
	 write(stdo,*)' in nbondm npt =',npt
	 write(stdo,1000)(poipt(i),i=1,totmon)
	end if
	do 1100 i=1,stmon2 - 1 
c set pointers, include self in list
		nlist2 = nlist2+1
		pointm2(i)=nlist2
		nblistm2(nlist2) = i


c loop over second monomer
		do 1200 j=i+1,stmon2

                        rxyz = mono_cent(1,i) - mono_cent(1,j)
                        r2 = rxyz*rxyz
                        if (r2.gt.cutebig2) go to 1200
                        rxyz = mono_cent(2,i) - mono_cent(2,j)
                        r2 = r2 + rxyz*rxyz
                        if (r2.gt.cutebig2) go to 1200
                        rxyz = mono_cent(3,i) - mono_cent(3,i)
                        r2 = r2 + rxyz*rxyz


c check distances^2. If smaller than cutebig2, include pair
c i,j in list, update pointers, go for next j.
			if(r2.le.cutebig2) then
				nlist2=nlist2+1
				nblistm2(nlist2)=j
			end if

1200		continue
		ki=spoipt2(i)+1

1100	continue

c set list for last monomer (only self)
	nlist2 = nlist2 + 1
	pointm2(stmon2)=nlist2
C the next sentence is just to insure the proper existence of
C pointm2(totmon)+1
	pointm2(stmon2+1)=nlist2 +1
	nblistm2(nlist2)=stmon2
	call nbond_spcl(pointm1,nblistm1,pointm2,nblistm2)

	return
	end
