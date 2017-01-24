      subroutine preDEE()
c
c Eliminates rotamers that have a very high energy with the backbone, 
c such as in cases of clashes. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c
c      
      character*6 name
      integer namel,level
c
      integer i,j,removed,removed_local,rotlastreal
      integer rot1,rot1real,rotfirst,rotlast
      integer nmonenhstart
      double precision Eibackmin,Eib
c
      name = 'preDEE'
      namel = 6
      level=1
c
c.....SUMMARY BEFORE DEE
      write (6,*) 
      write (6,*) 'Before DEE:'
      write (6,*) '  #rotamers ',nmonenhleft
      write (6,*) '  residues:'
      do j=1,nposenh
         rotfirst=poirotenhaux(j-1)+1
         rotlast=poirotenhaux(j)
         write (6,40) j,rotlast-rotfirst+1
      end do
 40   format(2x,'position #',i4,' (',i4,' rotamers)')
c
c
c
      nmonenhstart=nmonenhleft
c
      removed=0
c
c.....loop over all positions and eliminate clashes with the backbone
      do i=1,nposenh
c     
c........first pick out of this position rotamers with strong clashes
c........with the backbone (some concern about too large negative
c........energies -- possible in OPLS/AMBER force field due to zero 
c........vdW radius of some hidrogens). Also evaluate the lowest Eiback 
c........energy
c
         removed_local=0
         rotfirst=poirotenhaux(i-1)+1
         rotlast=poirotenhaux(i)
         Eibackmin=1.0d30
c
         do rot1=rotfirst,rotlast
c
            rot1real=rotaux(rot1)
            Eib=Eiback(rot1real)
c
            if (abs(Eib).ge.hugeE) then
c
               kept(rot1real)=.false.
               removed = removed+1
               removed_local = removed_local+1
c
            else
c
               Eibackmin=dmin1(Eib,Eibackmin)
c               
            end if
c
         end do
c
         if (removed_local.eq.rotlast-rotfirst+1) then
            write (6,*) 'Problems with clashes with the backbone at ',
     &           'position', i
            call alert(name,namel,'all rotamers discarded ',23,level)
         end if
c
c
c........ACTUALLY ELIMINATE THE ROTAMERS
         if (removed_local.gt.0) then
            call removerot(i,rotfirst,rotlast)
         end if
c
c
c........second pick out of this position rotamers with energy beyond 
c........a predefined window above the lowest Eiback.
c
         removed_local=0
         rotfirst=poirotenhaux(i-1)+1
         rotlast=poirotenhaux(i)
         Eibackmin=Eibackmin+Eibackmax
c
         do rot1=rotfirst,rotlast
c
            rot1real=rotaux(rot1)
c
            if (Eiback(rot1real).gt.Eibackmin) then
c
               kept(rot1real)=.false.
               removed = removed+1
               removed_local = removed_local+1
c
            end if
c
         end do
c
         if (removed_local.eq.rotlast-rotfirst+1) then
            write (6,*) 'Problems with clashes with the backbone at ',
     &           'position', i
            call alert(name,namel,'all rotamers discarded ',23,level)
         end if
c
c
c........ACTUALLY ELIMINATE THE ROTAMERS
         if (removed_local.gt.0) then
            call removerot(i,rotfirst,rotlast)
         end if
c
      end do
c
c
c
      if (removed.gt.0) then
c
         write (6,*) 
         write (6,*) 'Clashes with the backbone:'
         write (6,*) '  #rotamers start',nmonenhstart
         write (6,*) '  #rotamers end  ',nmonenhleft
         do j=1,nposenh
            rotfirst=poirotenhaux(j-1)+1
            rotlast=poirotenhaux(j)
            write (6,40) j,rotlast-rotfirst+1
         end do
c
      else
c
         write (6,*)
         write (6,*)'no rotamers removed because of clashes with ',
     &        'the backbone'
         write (6,*)
c
      end if
c
      return
      end
c     
c







