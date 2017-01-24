      subroutine rbackconn(uback) 
c
c Simplification  of rconn thar reads only  the relevant information 
c for DEE from the connectivity file of the original backbone molecule. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      
      integer uback
      character*9 name
      integer namel,level

      integer i,j,k,k1,iold,l,junk(30)
      integer nangl1,ntors1,nimp1
c
      character*3 version,vrs_cmp
      parameter (version='7.1')

      name = 'rbackconn'
      namel = 9
c
      rewind uback
c
      read(uback,1000,err=99,end=99)vrs_cmp
 1000 format(a3)
      if (vrs_cmp.ne.version) then
         level = 1
         call alert(name,namel,' con file too old',17,level)
      end if
      read(uback,*)
      read(uback,*)
c.....reads: totmon npt nb nangl ntors nimp totex totspe nbulk
      read(uback,1001)totmon,npt,nb,nangl1,ntors1,nimp1,totex,
     1        totspe,nbulk
 1001 format(1x,8i6,6x,i6)
c
      level = 1
      if (totmon.gt.maxmono) then
	 write(stdo,*)' totmon = ',totmon,'  maxmono = ',maxmono
	 call alert(name,namel,' number of mono exceeded ',25,level)
      else if (npt.gt.maxpt) then
	 write(stdo,*)' npt  = ',npt,'  maxpt = ',maxpt
	 call alert(name,namel,' number of prtc exceeded ',25,level)
      else if (nb.gt.maxbond) then
	 write(stdo,*)' nb = ',nb,'  maxbond = ',maxbond
	 call alert(name,namel,' number of bond exceeded ',25,level)
      else if(totex.gt.maxex) then
	 write(stdo,*)' totex = ',totex,'  maxex =',maxex
	 call alert(name,namel,' number of excl exceeded ',25,level)
      else if(totspe.gt.maxspec) then
	 write(stdo,*)' totspe = ',totspe,' maxspec = ',maxspec
	 call alert(name,namel,' number of spec exceeded ',25,level)
      else if(nbulk.gt.MAXBULK) then
	 write(stdo,*)' nbulk = ',nbulk,' MAXBULK = ',MAXBULK
	 call alert(name,namel,' number of BULK exceeded ',25,level)
      end if
c
      do i=1,nbulk
         read(uback,1002)BULK(i)
      end do
 1002 format(a)
c
      read(uback,*)
      read(uback,1003)(pbulk(i),i=1,nbulk)
 1003 format(10(1x,i5))
c
      read(uback,*)     
      read(uback,1004)(moname(i),i=1,totmon)
 1004 format(10(1x,a4))
      read(uback,*)
      read(uback,1005)(poipt(i),i=1,totmon)
 1005 format(10(1x,i5))
c
c.....Properties of particles list : 
c.....pt ptid lesid  ptnm   ptms   ptchg   epsgm6 epsgm12
c
      read(uback,*)
      read(uback,*)
      do 20 i=1, npt
         read(uback,1006) j,poimon(i),ptid(i),ptnm(i),ptms(i)
     1           ,ptchg(i),epsgm6(i),epsgm12(i)
 20   continue
 1006 format(1x,i5,1x,i5,1x,i3,5x,a4,1x,f7.2,2x,f9.5,1x,e12.5,e12.5)
c
c.....Bonds list: 
c.....ib1 ib2 kbond req 
      read(uback,*)
      read(uback,*)
      do 30 j= 1, nb
         read(uback,1007)ib1(j),ib2(j)
 30   continue
 1007 format(2(1x,i5),2(1x,f10.4))
c
c.....skip angles, torsions and improper torsions
      do i=1,nangl1+ntors1+nimp1+8
         read(uback,*)
      end do
c
c
c.....Exclusion list 1-2 1-3 1-4
c.....atom number, number of exclusions and list
      k = 0
      i = 0
      do while (k.lt.totex)
         read(uback,117)j,k1
 117     format(2(1x,i5))
         i = i + 1
         if (i.lt.j) then
            iold=i
            do i=iold,j-1
               exc1(i) = k
            end do
         else if (i.gt.j) then
            level = 1
            call alert(name,namel,'Fishy exclusion list!',21,level)
         end if
c........now i=j 
         exc1(j) = k + k1
         read(uback,1008)(exc2(l),l=k+1,exc1(j))
         k = k + k1
      end do
      do j=i+1,npt
         exc1(j) = exc1(j-1)
      end do
 1008 format(1x,10i5)
c
c.....Special list 1-4
      read(uback,*)
      read(uback,*)
      do i=1, totspe
         read(uback,1009)spec1(i),spec2(i)
      end do
 1009 format(2(1x,i5),3(1x,f15.5))
c
c.....Charge flags. T = charged atom , F = uncharged')
      read(uback,*)
      read(uback,1017)(flagchr(i),i=1,npt)
 1017 format(1x,15l2)
c
      return
c
c
 99   continue
      call alert(name,namel,'problems with reading',19,1)
c
      return
      end





