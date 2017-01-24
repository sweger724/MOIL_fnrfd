	subroutine nbsym_metal(pointm,nblistm,lsym,ia,jb,kc)
c
c generating non-bonded list for symmetry related atoms
c note that the neighbour list is simply added on the top
c of the existing nonbonded list. Therefore the generation
c of the primary atom list MUST be done before this one.
c Also, the pointer to the symmetry related particles
c Key vectors
c
      include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	integer pointm(0:*),nblistm(*)
	integer lsym,ia,jb,kc
C Check for EVERY monomer.
C Later change to only i+1,2,3 in no les case, or all in les case

	logical oneig1,oneig2,oneig3
	integer jbeg,jend,lstart,lend
	integer i,j,k,l,m,ki,kj
	double precision rx,ry,rz,r2
	double precision da,dc
	double precision xl,yl,zl

	integer namel
        character*5 name
        data namel,name/5,'nbsym'/

        save namel,name

	
	da = ia*a
	dc = kc*c

c yael - parallel
       if (prll_on_off.and.(my_pe.ne.0)) then
          ki = dpoipt(monsym(lsym,my_pe))+1
       else
          ki = 1
       endif


       if (prll_on_off) then
          lstart = monsym(lsym,my_pe)+1
          lend = monsym(lsym,my_pe+1)
       else
          lstart = 1
          lend = totdmon
       endif


	jbeg = 1
	do 1 i=lstart,lend
	  jend = pointm(i)


	    do 3 l=ki,dpoipt(i)

c metal image charges are for charged particle only
c skip loop if particle is not charged
c only list3 is updated
c
	     if (dabs(ptchg(l)).lt.1.d-6) go to 3
	     oneig1 = .false.
	     oneig2 = .false.
	     oneig3 = .false.

		xl = coor(1,l) + da
		if (coor(2,l).lt.0) then
			yl = -b - coor(2,l)
		else
			yl = b - coor(2,l)
		end if
		zl = coor(3,l) + dc

	     do 2 j=jbeg,jend
	      k=nblistm(j)
	      if (k.eq.1) then
		 kj = 1
	      else
		 kj = dpoipt(k-1)+1
	      end if
	        do 4 m=kj,dpoipt(k)
		   rx = xl - coor(1,m)
		   ry = yl - coor(2,m)
		   rz = zl - coor(3,m)
		   r2 = rx*rx + ry*ry + rz*rz
c
c Only third list is used for metal images
c
c
                if (r2.lt.cutebig2 .and.
     1              (flagchr(l).and.flagchr(m))) then
		   	 oneig3 = .true.
                   nlist3 = nlist3 + 1
                   list3(nlist3) = m
                end if

4	        continue

	   if (nlist1.gt.ichgvdw) then
                write(*,*)' nlist1 ichgvdw ',nlist1,ichgvdw
                call alert(name,namel,'List1 too short',15,1)
         else if (nlist2.gt.ivdw) then
                write(*,*)' nlist2 ivdw ',nlist2,ivdw
                call alert(name,namel,'List2 too short',15,1)
         else if (nlist3.gt.ichg) then
                write(*,*)' nlist3 ichg ',nlist3,ichg
                call alert(name,namel,'List3 too short',15,1)
         end if

c
c psym is the pointer to lista. Thus psym(i) is the position in lista (lista(psym(i)))
c that is equal to the last atom interacting with the symmetry generated atoms - i
c The lista vector therefore includes the lista of real atom neigbors and neighbors
c to the symmetry generated particles as well.
c iblock(j) is the index of the last atom that is included in
c the current symmetry operation j
c symreal(k) is the corresponding real particle index of symmetry generated
c particle with index k
c
2	  continue
	     if (oneig1) then
	      indxsym1         = indxsym1 + 1
	      symreal1(indxsym1) = l
	     end if
	     if (oneig2) then
	      indxsym2         = indxsym2 + 1
	      symreal2(indxsym2) = l
	     end if
	     if (oneig3) then
	      indxsym3         = indxsym3 + 1
	      symreal3(indxsym3) = l
	     end if
	     if (indxsym1.ne.0) psym1(indxsym1)    = nlist1
	     if (indxsym2.ne.0) psym2(indxsym2)    = nlist2
	     if (indxsym3.ne.0) psym3(indxsym3)    = nlist3
3 	     continue
	  ki = dpoipt(i) + 1
	  jbeg = pointm(i) + 1
1	continue

	iblock1(lsym)     = indxsym1
	iblock2(lsym)     = indxsym2
	iblock3(lsym)     = indxsym3
	return
	end
