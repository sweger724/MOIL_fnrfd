      subroutine nbsym(pointm,nblistm,lsym,ia,jb,kc)
c     
c     generating non-bonded list for symmetry related atoms
c     note that the neighbour list is simply added on the top
c     of the existing nonbonded list. Therefore the generation
c     of the primary atom list MUST be done before this one.
c     Also, the pointer to the symmetry related particles
c     Key vectors
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
C     Check for EVERY monomer.
C     Later change to only i+1,2,3 in no les case, or all in les case

      logical oneig1,oneig2,oneig3,flag2
      integer jbeg,jend,lstart,lend
      integer i,j,k,l,m,ki,kj
      double precision rx,ry,rz,r2
      double precision da,db,dc
      double precision xl,yl,zl

      integer namel
      character*5 name
      data namel,name/5,'nbsym'/

      save namel,name

      
      da = ia*a
      db = jb*b
      dc = kc*c

c     yael - parallel
      if (prll_on_off.and.(my_pe.ne.0)) then
         ki = dpoipt(monsym(lsym,my_pe))+1
      else
         ki = 1
      endif

      if (lsym.eq.1) then
         if (prll_on_off) then
            nlist1 = point1(dpoipt(monp(my_pe+1)))
            nlist2 = point2(dpoipt(monp(my_pe+1)))
            nlist3 = point3(dpoipt(monp(my_pe+1)))
         else
            nlist1 = point1(npt-1)
            nlist2 = point2(npt-1)
            nlist3 = point3(npt-1)
         endif
      end if

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

            oneig1 = .false.
            oneig2 = .false.
            oneig3 = .false.

            xl = coor(1,l) + da
            yl = coor(2,l) + db
            zl = coor(3,l) + dc

            do 2 j=jbeg,jend
               k=nblistm(j)
               if (k.eq.1) then
                  kj = 1
               else
                  kj = dpoipt(k-1)+1
               end if
               do 4 m=kj,dpoipt(k)
                     flag2=((lesid(l).eq.lesid(m)).and.
     *               (cplbl(l).ne.cplbl(m))).and.
     *               (lesid(l).ne.0.and.lesid(m).ne.0)
                     if (flag2) goto 4

                  rx = xl - coor(1,m)
                  ry = yl - coor(2,m)
                  rz = zl - coor(3,m)
                  r2 = rx*rx + ry*ry + rz*rz
c     
c     check for hydrogens. Hydrogens have zero van-der Waals rad if arith=.false.
c     (default) and should be placed in the third (electrostatic only)
c     list
c     
                  flag2 = (.not.arith)
     $                 .and.(epsgm12(l).lt.1.d-3.or.epsgm12(m).lt.1.d-3)
     $                 .and.(r2.lt.cutebig2)
                  if (flag2 .and. flagchr(l) .and. flagchr(m)) then
                     oneig3 = .true.
                     nlist3 = nlist3 + 1
                     list3(nlist3) = m
                     go to 4
                  else if (flag2) then
                     go to 4
                  end if


c     make now a decision to which of the three
c     lists this pair belongs: list1 includes both
c     van der Waals and electrostatic and the max.
c     distance is cutvbig2. list2 uses the same
c     cutoff (cutvbig2) but intends to uncharged particle only.
c     list3 is for for distances larger than cutvdw2 and smaller
c     than cutebig2 and includes only electrostatic.
c     
                  if(r2.lt.cutvbig2)then
                     if (flagchr(l).and.flagchr(m)) then
                        oneig1 = .true.
                        nlist1 = nlist1 + 1
                        list1(nlist1) = m
                     else
                        oneig2 = .true.
                        nlist2 = nlist2 + 1
                        list2(nlist2) = m
                     end if
                  else if (r2.lt.cutebig2) then
                     if (flagchr(l).and.flagchr(m)) then
                        oneig3 = .true.
                        nlist3 = nlist3 + 1
                        list3(nlist3) = m
                     end if
                  end if

 4             continue

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
c     psym is the pointer to lista. Thus psym(i) is the position in lista (lista(psym(i)))
c     that is equal to the last atom interacting with the symmetry generated atoms - i
c     The lista vector therefore includes the lista of real atom neigbors and neighbors
c     to the symmetry generated particles as well.
c     iblock(j) is the index of the last atom that is included in
c     the current symmetry operation j
c     symreal(k) is the corresponding real particle index of symmetry generated
c     particle with index k
c     
 2          continue
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
 3       continue
         ki = dpoipt(i) + 1
         jbeg = pointm(i) + 1
 1    continue

      iblock1(lsym)     = indxsym1
      iblock2(lsym)     = indxsym2
      iblock3(lsym)     = indxsym3
      return
      end



