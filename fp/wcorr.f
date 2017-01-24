      program wcorr

c     this program calculates the correlation function of water-water
C     distances in a given input file. Picked atoms are assumed to
C     belong to the peptide, and they are assumed to be listed at the
C     beginning of each structure, so that the first atom of the first
C     water is atom "npick + 1".

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'


      character*10 name
      integer n, i, j, k,pointr(maxpt),ipick(maxpt),level
      integer namel,npick,ucon,urcrd,of,geti,nstr,nwat
      logical find, fopen
      double precision getd, e, coor2(3,maxpt),r2

      name   = 'wcorr'
      namel  = 5
      stdi   = 5
      stdo   = 6
      stderr = 0
      nstr   = 1

c     junk file (necessary)
      jnkf = 25
      open (unit=jnkf,status='scratch')

c     begin input loop
 1    continue
      call rline(name,namel,stdi)
      if (find('file')) then
         if (find('conn')) then
            ucon = of()
            call rconn(ucon)
         endif

         if (find('rcrd')) urcrd = of()
      end if

      if (find('selc')) then
         call pick(ipick,npick)
         npick=0
         do 3 i=1,npt
            npick=npick+ipick(i)
 3       continue
         j=0
         do 4 i=1,npt
            if (ipick(i).gt.0) then
               j=j+1
               pointr(j)=i
            end if
 4       continue
      end if

      nstr   = geti('nstr',nstr)

      if (find('acti')) go to 2
      go to 1
 2    continue
c     end input loop---------------------------------------

      if (.not.fopen(urcrd)) then
         level = 1
         call alert(name,namel,'urcrd not opened',15,level)
      end if

      nwat = npt - npick

c     use coor2 as initial structure
      read(urcrd) e, ((coor2(j,i),i=1,npt),j=1,3)

c     correlations in water-water distance
      do 10 n = 2,nstr
         read(urcrd) e, ((coor(j,i),i=1,npt),j=1,3)
         r2 = 0.d0
c        measure distance based on oxygens
         do 11 i = (npick+1),npt,3
            do 12 j = (i+3),npt,3
               do 13 k = 1,3
                  r2 = r2 + (coor2(k,i) - coor(k,j))**2
 13            continue
 12         continue
 11      continue
         r2 = r2 / nwat
         write(*,*) 'n, r2:', n, r2
 10   continue

      end
