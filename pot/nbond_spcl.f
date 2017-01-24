	subroutine nbond_spcl(pointm1,nblistm1,pointm2
     1  ,nblistm2)
      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/CONSPECL1.BLOCK'
      	include 'COMMON/CONSPECL2.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	integer pointm1(*),nblistm1(*)
	integer pointm2(*),nblistm2(*)
	double precision rxyz,r2
C Check for EVERY monomer.

	logical flag2
	integer i,j,k,l,m,ki,kj,jbeg,jend
	integer nlist1,nlist2
	integer jx,ll,mm
	integer jbegx,jendx


c list of nonbond interaction when ligand did not bond to MB	
	nlist1 = 0
	ki = 1
	do 100 i=1,stmon1
          jbeg = pointm1(i)
	  jend = pointm1(i+1) - 1

	  do 200 l=ki,spoipt1(i)
	     kj = 1

	     do 300 j=jbeg,jend
	        k=nblistm1(j)
		kj = max(kj,l)

	        do 400 m=kj,spoipt1(k)


c a particle cannot be a neigbour to itself
	if (l.eq.m) goto 400 
c skip test if bothe particles are frozen
        ll = newpoit(l)
	mm = newpoit(m)
	if (zerofrz(ll).eq.0 .and. zerofrz(mm).eq.0) go to 400
       		    flag2=((lesid(ll).eq.lesid(mm)).and.
     *               (spoimon1(l).ne.spoimon1(m))).and.
     *		     (lesid(ll).ne.0.and.lesid(mm).ne.0)
		     if (flag2) goto 400
c------------------------------------------------------------------
                    if (repyes0) then
                     if((ll.eq.imb1(1)) .or.(ll.eq.imb1(2)).or.
     2                 (ll.eq.imb1(3)) .or.(ll.eq.imb1(4))) then
                      if((mm.eq.imb2(1)) .or. (mm.eq.imb2(2)).or.
     1                (mm.eq.imb2(3)).or.(mm.eq.imb2(4))) go to 400
                     else if((mm.eq.imb1(1)) .or.(mm.eq.imb1(2)).or.
     2                 (mm.eq.imb1(3)) .or.(mm.eq.imb1(4))) then
                      if((ll.eq.imb2(1)) .or. (ll.eq.imb2(2)).or.
     1                (ll.eq.imb2(3)).or.(ll.eq.imb2(4))) go to 400
                     end if
                    end if
c--------------------------------------------------------------
		      jbegx=sexc11(l-1)+1
	              jendx=sexc11(l)
	              if (jbegx.le.jendx) then
		       do 1500 jx=jbegx,jendx
			if(sexc12(jx).eq.m) then
				go to 400
			end if
1500		       continue
		      end if
		        rxyz = coor(1,ll)-coor(1,mm)
			r2 = rxyz*rxyz
			if (r2.ge.cutebig2) go to 400
		        rxyz = coor(2,ll)-coor(2,mm)
			r2 = r2 + rxyz*rxyz
			if (r2.ge.cutebig2) go to 400
			rxyz = coor(3,ll)-coor(3,mm)
			r2 = r2 + rxyz*rxyz
c Since the "special" part usually includes short lists
c we do not separate the lists to electrostatic and van der Waasl
c only, here all neighbours are joined into a single list

                if(r2.lt.cutebig2)then
                   nlist1 = nlist1 + 1
                   slist11(nlist1) = m
                end if
                  
400	        continue

	        kj = spoipt1(k) + 1
300 	     continue
             spoint11(l)=nlist1

200	continue

	     ki = spoipt1(i) + 1
100	continue

c-----------------------------------------------------------

c list of nonbond interaction when ligand bonded to MB	

	nlist2 = 0
	ki = 1
	do 600 i=1,stmon2
          jbeg = pointm2(i)
	  jend = pointm2(i+1) - 1

	  do 800 l=ki,spoipt2(i)
	     kj = 1

	     do 900 j=jbeg,jend
	        k=nblistm2(j)
		kj = max(kj,l)

	        do 5000 m=kj,spoipt2(k)


c a particle cannot be a neigbour to itself
	if (l.eq.m) goto 5000 
c skip test if bothe particles are frozen
        ll = newpoit(l)
	mm = newpoit(m)
	if (zerofrz(ll).eq.0 .and. zerofrz(mm).eq.0) go to 5000
       		    flag2=((lesid(ll).eq.lesid(mm)).and.
     *               (spoimon2(l).ne.spoimon2(m))).and.
     *		     (lesid(ll).ne.0.and.lesid(mm).ne.0)
		     if (flag2) goto 5000
		      jbegx=sexc21(l-1)+1
	              jendx=sexc21(l)
	              if (jbegx.le.jendx) then
		       do 2500 jx=jbegx,jendx
			if(sexc22(jx).eq.m) then
				go to 5000
			end if
2500		       continue
		      end if
                        rxyz = coor(1,ll)-coor(1,mm)
                        r2 = rxyz*rxyz
                        if (r2.ge.cutebig2) go to 5000
                        rxyz = coor(2,ll)-coor(2,mm)
                        r2 = r2 + rxyz*rxyz
                        if (r2.ge.cutebig2) go to 5000
                        rxyz = coor(3,ll)-coor(3,mm)
                        r2 = r2 + rxyz*rxyz
c Since the "special" part usually includes short lists
c we do not separate the lists to electrostatic and van der Waasl
c only, here all neighbours are joined into a single list

                if(r2.lt.cutebig2)then
                   nlist2 = nlist2 + 1
                   slist21(nlist2) = m
                end if
                  
5000	        continue

	        kj = spoipt2(k) + 1
900 	     continue
             spoint21(l)=nlist2

800	continue

	     ki = spoipt2(i) + 1
600	continue
	return
	end
