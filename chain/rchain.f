        subroutine rchain(urcrd,r,igrid,styl,pointr,npick)

        implicit none
c
c a subrotuine to read initial chain coordinates. Three styles
c are supported: 
c (i)   DYNAmics (binary single precision format for compatability
c       with QUANTA)
c (ii)  PATH (double precision format recommend for computations.
c (iii) INIT (initialize a chain using linear interpolation, urcrd
c               contains the names of normal coordinate files for
c               reactants and products).
c (iv)  INTR (interpolation. Missing structures are interpolated
c             (linearly) between existing structure. Here urcrd IS NOT a
c             coordinate file but rather a set of directions
c             written in a free format (the crd file must be PATH)
c
c             [number of old structures to be read]
c             [indices of the old structures in the new set- 1 4 5.. ]
c             [name of old coordinate file]
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CCRD.BLOCK'

        integer urcrd,igrid
        logical find
        integer of
        integer npick
        integer pointr(*)
        double precision r(3,*)
        character*4 styl

c local
        integer i,j,k,l,i1,i2,npt3,npt2,namel,level
        integer u1,u2
        integer nstru,istru(lgrid)
        integer nofreez(maxpt)
        integer rbin
        double precision div,rms, dummy(3,maxpt),dummy2(3,maxpt)

        character*6 name

        rbin = 1

        norew = .false.
        lpstr = 1

        name  = 'rchain'
        namel = 6

        npt3 = 3*npt
        npt2 = 2*npt

        if (styl.eq.'DYNA') then
         i1 = npt
         do 1 i=1,igrid
          rewind urcrd
          j = (i-1)*npt
          call rdyncrd(urcrd,j,i1,nofreez,rbin)
          do 1 k=1,npt
           l = j + k
           r(1,l)  = coor(1,k)
           r(2,l)  = coor(2,k)
           r(3,l)  = coor(3,k)
1         continue
         do 15 i=1,igrid-1
          j = (i-1)*npt + 1
          k = i*npt + 1
          call ovrlpck(r(1,j),r(1,k),dummy,dummy2,pointr,npick,rms)
15       continue

        else if (styl.eq.'PATH') then
         rewind urcrd
         do 2 i=1,igrid
          call rpath(urcrd,i)
          j = (i-1)*npt
          do 2 k=1,npt
           l = j + k
           r(1,l)  = coor(1,k)
           r(2,l)  = coor(2,k)
           r(3,l)  = coor(3,k)
2        continue
         do i=1,igrid-1
           j = (i-1)*npt + 1
           k = i*npt + 1
           call ovrlpck(r(1,j),r(1,k),dummy,dummy2,pointr,npick,rms)
         end do

         else if (styl.eq.'INIT') then
          if (urcrd.ne.stdi) rewind urcrd
c getting unit number of reactants (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
           u1 = of()
          else
           level = 1
           call alert(name,namel,'Missing file name',17,level)
          end if
c getting unit number of products (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
           u2 = of()
          else
           level = 1
           call alert(name,namel,'Missing file name',17,level)
          end if
c reading reactants
          call getcrd(u1,'CHARM')
          do 3 i=1,npt
           r(1,i)   = coor(1,i)
           r(2,i)   = coor(2,i)
           r(3,i)   = coor(3,i)
3         continue
c reading products
          call getcrd(u2,'CHARM')
          j = (igrid-1)*npt
          do 4 i=1,npt
           l = j + i
           r(1,l)    = coor(1,i)
           r(2,l)    = coor(2,i)
           r(3,l)    = coor(3,i)
4         continue
c overlapping products with respect to reactants
          call ovrlpck(r(1,1),r(1,j+1),dummy,dummy2,pointr,npick,rms)
c calculating step size.
          div = dble(1.d0/(igrid-1))
          do 5 i=1,npt
           l = j + i
           coor(1,i) = (r(1,l)-r(1,i))*div
           coor(2,i) = (r(2,l)-r(2,i))*div
           coor(3,i) = (r(3,l)-r(3,i))*div
5         continue
c calculate intermediate structures
          do i=1,igrid-1
           k = i*npt
           l = (i-1)*npt
           do j=1,npt
            i1 = k + j
            i2 = l + j
            r(1,i1)  = r(1,i2) + coor(1,j)
            r(2,i1)  = r(2,i2) + coor(2,j)
            r(3,i1)  = r(3,i2) + coor(3,j)
           end do
          end do
        else if (styl.eq.'INTR') then
         do 7 i=1,lgrid
          istru(i) = 0
7        continue
         if (urcrd.ne.stdi)  rewind urcrd
         read(urcrd,*,err=999)nstru
         read(urcrd,*,err=999)(istru(i),i=1,nstru)
         if (istru(1).ne.1) then
          level = 1
          call alert(name,namel,'Missing reactants',17,level)
         else if (istru(nstru).ne.igrid) then
          level = 1
          call alert(name,namel,'Missing products',16,level)
         end if
         call rline(name,namel,urcrd)
         if (find('file')) then
          u1 = of()
         else
          level = 1
          call alert(name,namel,'Missing file name',17,level)
         end if
         do 8 i=1,nstru
          call rpath(u1,i)
          j = (istru(i)-1)*npt
          do 8 k=1,npt
           l = j + k
           r(1,l)  = coor(1,k)
           r(2,l)  = coor(2,k)
           r(3,l)  = coor(3,k)
8        continue
         do 11 i=1,nstru-1
          j = (istru(i)-1)*npt  + 1
          k = (istru(i+1)-1)*npt + 1
          l = istru(i)*npt+1
          call ovrlpck(r(1,j),r(1,k),dummy,dummy2,pointr,npick,rms)
          do 9 i1=0,npt-1
           r(1,l+i1) = (r(1,j+i1)+r(1,k+i1))*0.5d0
           r(2,l+i1) = (r(2,j+i1)+r(2,k+i1))*0.5d0
           r(3,l+i1) = (r(3,j+i1)+r(3,k+i1))*0.5d0
9         continue
11        continue
         end if
         return
999      continue
         level = 1
         call alert(name,namel,'Error reading nstru data',24,level)
         return
         end
