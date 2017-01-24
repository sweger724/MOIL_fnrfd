      subroutine DEE1()
c
c Eliminates rotamers using the Dead-End-Elimination algorithm. 
c In this routine we use the inequality for rotamers i and j. 
c    (Goldstein, Biophysical J. 66,1335(1994)
c    E(g_i)-<E(h_i)>+SUM_{j}{ min_{f'} [ E(g_i,f'_j)-<E(h_i,f'_j)>]}
c    
c       <A>=SUM_{k.ne.g} {C_k * A} ,
c           where  SUM_{k.ne.g} {C_k} =1 ,C_k>=0
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

      
      character*4 name
      integer namel,level
c
      integer nrotlocalmax
      parameter (nrotlocalmax=(maxnmonenh/maxnposenh)*10)
c
      double precision dtotal,ediffmin,Eibackave,Eijave
      integer i,j,removed,ind1,ind2,ind2temp,niterations,rotlastreal
      integer rot1,rotfirst,rotlast,rot1real,rot2realtemp
      integer rot2,rot2real(nrotlocalmax)
      integer rotT,rotfirstT,rotlastT,rotTreal,nrotlocal
      integer nmonenhstart
c
      name = 'DEE1'
      namel = 4
      level=1
c
c
      removed=999
      niterations=0
      nmonenhstart=nmonenhleft
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
c...........check if there is more than one left
            if (rotlast.gt.rotfirst) then
c
               do rot1=rotfirst,rotlast
c
                  rot1real=rotaux(rot1)
                  nrotlocal=0
c
                  do rot2=rotfirst,rot1-1
c
                     rot2realtemp=rotaux(rot2)
                     if (kept(rot2realtemp)) then
                        nrotlocal=nrotlocal+1
                        if (nrotlocal.le.nrotlocalmax) then
                           rot2real(nrotlocal)=rot2realtemp
                        else
                           call alert(name,namel,
     &                          'nrotlocal.gt.nrotlocalmax',16,level)
                        end if
                     end if
c
                  end do
c
c
                  do rot2=rot1+1,rotlast
c
                     rot2realtemp=rotaux(rot2)
                     if (kept(rot2realtemp)) then
                        nrotlocal=nrotlocal+1
                        if (nrotlocal.le.nrotlocalmax) then
                           rot2real(nrotlocal)=rot2realtemp
                        else
                           call alert(name,namel,
     &                          'nrotlocal.gt.nrotlocalmax',16,level)
                        end if
                     end if
c
                  end do
c
                  if (nrotlocal.gt.0) then
c
                     Eibackave=0.d0
                     do rot2=1,nrotlocal
                        Eibackave=Eibackave+Eiback(rot2real(rot2))
                     end do
                     Eibackave=Eibackave/dble(nrotlocal)
c
                     dtotal=Eiback(rot1real)-Eibackave
c
                     do j=1,i-1
c     
                        ediffmin=1.0d30
c     
                        rotfirstT=poirotenhaux(j-1)+1
                        rotlastT=poirotenhaux(j)
c
                        do rotT=rotfirstT,rotlastT
c
                           rotTreal=rotaux(rotT)
                           ind1=pointEij(rotTreal-1)+rot1real-
     &                          poirotenh(j)
c
                           ind2temp=pointEij(rotTreal-1)-poirotenh(j)
c
                           Eijave=0.d0
                           do rot2=1,nrotlocal
                              ind2=ind2temp+rot2real(rot2)
                              Eijave=Eijave+Eij(ind2)
                           end do
                           Eijave=Eijave/nrotlocal
c
                           ediffmin=dmin1(ediffmin,Eij(ind1)-Eijave)
c
                        end do
c     
                        dtotal=dtotal+ediffmin
c                     
                     end do
c     
                     do j=i+1,nposenh
c     
                        ediffmin=1.0d30
c     
                        rotfirstT=poirotenhaux(j-1)+1
                        rotlastT=poirotenhaux(j)
c     
                        do rotT=rotfirstT,rotlastT
c
                           rotTreal=rotaux(rotT)
                           ind1=pointEij(rot1real-1)+rotTreal-
     &                                poirotenh(i)
c
                           ind2temp=rotTreal-poirotenh(i)
c
                           Eijave=0.d0
                           do rot2=1,nrotlocal
                              ind2=pointEij(rot2real(rot2)-1)+ind2temp
                              Eijave=Eijave+Eij(ind2)
                           end do
                           Eijave=Eijave/nrotlocal
c
                           ediffmin=dmin1(ediffmin,Eij(ind1)-Eijave)
c
                        end do
c     
                        dtotal=dtotal+ediffmin
c                     
                     end do
c
                  end if
c
                  if (dtotal.gt.0) then
                     kept(rot1real)=.false.
                     removed = removed+1
                  end if
c
c
               end do
c
c
c..............ACTUALLY ELIMINATE THE ROTAMERS
c
               if (removed.gt.0) then
                  call removerot(i,rotfirst,rotlast)
               end if
c
c
            end if
c
         end do
c
      end do
c
      write (6,*) 
      write (6,*) 'DEE procedure 2:'
      write (6,*) '  #iteration',niterations
      write (6,*) '  #rotamers start',nmonenhstart
      write (6,*) '  #rotamers end  ',nmonenhleft
      write (6,*) '  final residues:'
      do j=1,nposenh
         rotfirst=poirotenhaux(j-1)+1
         rotlast=poirotenhaux(j)
         write (6,40) j,rotlast-rotfirst+1
         write (6,50)(rotaux(i),monames4(typerotenh(rotaux(i))),
     &        i=rotfirst,rotlast)
      end do
CDEB
CDEB      write (6,*) 'all residues flags'
CDEB      write (6,*)(kept(i),i=1,nmonenh)
c
 40   format(2x,'position #',i4,' (',i4,' rotamers)')
 50   format(8(1x,i4,1x,a4))
c
c
      return
      end
c     
c







