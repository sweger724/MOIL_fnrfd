	subroutine addbond(ubond)
c
c subroutine for adding bonds which are hard to get
c otherwise. e.g. linking two molecules. Absolute addresses
c of atoms is required in one option, another option is specifing
c the monomer number and unique particle name (to that monomer)
c be added this restriction will be omitted.
c Input line to ubond, should have the format:
c bond atm1=[intg] atm2=[intg]
c  OR
c bond chem mono index [unique particle] mono index [unique particle]
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/MONOMERS.BLOCK'
	include 'COMMON/PROPERT.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/ENERGY.BLOCK'

	double precision getd
	integer ubond
	integer geti
	logical find,unbind
	integer intgr,search
c local
	character*7 name
	character*4 m1,m2,pt1,pt2
	integer ist1,ist2,ien1,ien2,i
	integer namel,id1,id2,level,im1,im2
	integer idm1, idm2
	integer i1,i2
	logical empty
	data name/'addbond'/
	data namel/7/

	unbind = .false.
	emyes0 = .false.
	nmb    = 0

	if (ubond.ne.5) rewind ubond
1	continue
	call rline(name,namel,ubond)
	if (find('*EOD')) then
	  if(nmb.gt.0) then
	     do 200 i=1,nmb
	        write(stdo,101)i,rmeq(i),imb1(i),imb2(i)
200	     continue 
101	     format(1x,'morse bond #no=',i4,'equilibrium distance'
     1             ,f5.2,' prt1 - prt2 ',2i7)
	  end if
	  return
	end if
c ==============================================================
        if (find('unbd')) unbind=.true.
	if (find('mors')) then
	  emyes0=.true.
	  idm1=geti('atm1',-1)
	  idm2=geti('atm2',-1)
	  if (idm1.lt.0 .or. idm2.lt.0) then
	   level = 1
	   call alert(name,namel,'Illegal morse bond input',18,level)
	  end if
	  nmb=nmb+1
	  rmeq(nmb) = getd('requ',-10.d0)
	  D(nmb)    = getd('Dmor',D(nmb))
	  alpha(nmb)= getd('alph',alpha(nmb))
	 if (idm1.lt. idm2 ) then
	  imb1(nmb) = idm1
	  imb2(nmb) = idm2
	 else
	  imb1(nmb) = idm2
	  imb2(nmb) = idm1
	 end if
	 if (unbind) go to 999
	 if (ptid(idm1).lt.ptid(idm2)) then
	  i1 = ptid(idm1)
	  i2 = ptid(idm2)
	 else
	  i1 = ptid(idm2)
	  i2 = ptid(idm1)
	 end if
	 idm1 = search(i1,i2,0,0,bondp(1,1),bondp(1,2),
     1    0,0,maxubnd,2)
	 rmeq(nmb)   = basereq(idm1)
999      continue
	  if (rmeq(nmb).lt.0.d0) then
	   level = 1
	   call alert(name,namel,'requ not found',13,level)
	  end if
	end if
cjarek
c	if (nmb.gt.0) then
c	 write(stdo,101)nmb,rmeq(nmb)
c101	 format(1x,'number of morse bonds= ',i5,' equilibrium distance'
c     1	  ,f5.2)
c	end if
c ===================================================================
c
	if (find('bond')) then
	 id1 = geti('atm1',-1)
	 id2 = geti('atm2',-1)
	 if (id1.lt.0 .or. id2.lt.0) then

	  if (find('chem')) then
	   call get4c(m1,empty)
	   if (empty) then
	    level = 1
	    call alert(name,namel,'Missing atom name',17,level)
	   end if
	   im1 = intgr()
	   call get4c(pt1,empty)
	   if (empty) then
	    level = 1
	    call alert(name,namel,'Missing atom name',17,level)
	   end if
	   call get4c(m2,empty)
	   if (empty) then
	    level = 1
	    call alert(name,namel,'Missing atom name',17,level)
	   end if
	   im2 = intgr()
	   call get4c(pt2,empty)
	   if (empty) then
	    level = 1
	    call alert(name,namel,'Missing atom name',17,level)
	   end if 
	   if (m1.ne.moname(im1) .or. m2.ne.moname(im2)) then
	    level = 1
	    write(stdo,100)m1,moname(im1),m2,moname(im2)
100	    format(1x,'monomer name does not match, read ',a4,
     1       ' expected ',a4,' read ',a4,' expected ',a4)
	    call alert(name,namel,'Mono do not match',17,level)
	   end if

	   ien1 = poipt(im1)
	   if (im1.eq.1) then
	    ist1 = 1
	   else
	    ist1 = poipt(im1-1)+1
	   end if

	   ien2 = poipt(im2)
	   if (im2.eq.1) then
	    ist2 = 1
	   else
	    ist2 = poipt(im2-1)+1
	   end if

	   do 2 i=ist1,ien1
	    if(pt1.eq.ptnm(i)) then
	     id1 = i
	     go to 3
	    end if
2	   continue
	   level = 1
	   call alert(name,namel,'Particle not found',18,level)
3	   continue
	   
	   do 4 i=ist2,ien2
	    if(pt2.eq.ptnm(i)) then
	     id2 =i
	     go to 5
	    end if
4	   continue
	   level = 1
	   call alert(name,namel,'Particle not found',18,level)
5	   continue
	  else  
	   level = 1
	   call alert(name,namel,'Illegal bond input',18,level)
	  end if

	 end if
	 nb = nb + 1
	 if (id1.lt.id2) then
	  ib1(nb) = id1
	  ib2(nb) = id2
	 else
	  ib1(nb) = id2
	  ib2(nb) = id1
	 end if
	end if
	 if (ptid(ib1(nb)).lt.ptid(ib2(nb))) then
	  i1 = ptid(ib1(nb))
	  i2 = ptid(ib2(nb))
	 else
	  i1 = ptid(ib2(nb))
	  i2 = ptid(ib1(nb))
	 end if
	id1 = search(i1,i2,0,0,bondp(1,1)
     1   ,bondp(1,2),0,0,maxubnd,2)
	kbond(nb) = basekb(id1)
	req(nb)   = basereq(id1)
	go to 1
	end
