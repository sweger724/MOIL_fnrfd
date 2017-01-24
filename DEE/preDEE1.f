      subroutine preDEE1()
c
c Eliminates rotamers that have clashes with all rotamers of some 
c position. 
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

      
      character*3 name
      integer namel,level
c
      double precision e1jmin,Eij1
      integer i,j,removed,removed_local,ind1,niterations,rotlastreal
      integer rot1,rotfirst,rotlast,rot1real
      integer rotT,rotfirstT,rotlastT,rotTreal
      integer nmonenhstart
c
      name = 'DEE'
      namel = 3
      level=1
c
c
      removed=999
      niterations=0
      nmonenhstart=nmonenhleft
c
      write (6,*) 
      write (6,*) 'Side Chain-Side Chain clashes:'
c
c
      do while (removed.gt.0)
c
         removed=0
         niterations=niterations+1
c
c........loop over all positions
         do i=1,nposenh
c     
            rotfirst=poirotenhaux(i-1)+1
            rotlast=poirotenhaux(i)
c
            do j=1,i-1
c
               rotfirstT=poirotenhaux(j-1)+1
               rotlastT=poirotenhaux(j)
               removed_local=0
c
               do rot1=rotfirst,rotlast
c
                  rot1real=rotaux(rot1)
                  e1jmin=1.0d30
c
                  do rotT=rotfirstT,rotlastT
                     rotTreal=rotaux(rotT)
                     ind1=pointEij(rotTreal-1)+rot1real-
     &                    poirotenh(j)
                     Eij1=abs(Eij(ind1))
c
                     e1jmin=dmin1(e1jmin,Eij1)
c
                  end do
c     
                  if (e1jmin.gt.Eibackmax) then
                     kept(rot1real)=.false.
                     removed = removed+1
                     removed_local=removed_local+1
                  end if
c                     
               end do
c     
               if (removed_local.eq.rotlast-rotfirst+1) then
                  write (6,*) 'Problems with clashes with the ',
     &                 'backbone at position', i
                  call alert(name,namel,'all rotamers discarded ',
     &                 23,level)
               end if
c
c..............ACTUALLY ELIMINATE THE ROTAMERS
c
               if (removed_local.gt.0) then
                  call removerot(i,rotfirst,rotlast)
               end if
c
CDEB
CDEB               write (6,*) ' (#removed local',removed_local,')'
CDEB               write (6,*) ' (#removed ',removed,')',i,j
c
            end do
c
            do j=i+1,nposenh
c
               rotfirstT=poirotenhaux(j-1)+1
               rotlastT=poirotenhaux(j)
               removed_local=0
c
               do rot1=rotfirst,rotlast
c
                  rot1real=rotaux(rot1)     
                  e1jmin=1.0d30
c
                  do rotT=rotfirstT,rotlastT
                     rotTreal=rotaux(rotT)
                     ind1=pointEij(rot1real-1)+rotTreal-
     &                    poirotenh(i)
                     Eij1=abs(Eij(ind1))
c
                     e1jmin=dmin1(e1jmin,Eij1)
c
                  end do
c     
                  if (e1jmin.gt.Eibackmax) then
                     kept(rot1real)=.false.
                     removed = removed+1
                     removed_local=removed_local+1
                  end if
c                     
               end do
c     
               if (removed_local.eq.rotlast-rotfirst+1) then
                  write (6,*) 'Problems with clashes with the ',
     &                 'backbone at position', i
                  call alert(name,namel,'all rotamers discarded ',
     &                 23,level)
               end if
c
c..............ACTUALLY ELIMINATE THE ROTAMERS
c
               if (removed_local.gt.0) then
                  call removerot(i,rotfirst,rotlast)
               end if
c
CDEB
CDEB               write (6,*) ' (#removed local',removed_local,')'
CDEB               write (6,*) ' (#removed ',removed,')',i,j
c
            end do
c
         end do
c
         write (6,*) '  finished iteration',niterations,
     &        ' (#removed',removed,')'
c
      end do
c
      write (6,*) '  #iteration',niterations
      write (6,*) '  #rotamers start',nmonenhstart
      write (6,*) '  #rotamers end  ',nmonenhleft
      write (6,*) '  final residues:'
      do j=1,nposenh
         rotfirst=poirotenhaux(j-1)+1
         rotlast=poirotenhaux(j)
         write (6,40) j,rotlast-rotfirst+1
      end do
 40   format(2x,'position #',i4,' (',i4,' rotamers)')
c
      return
      end
c     
c







