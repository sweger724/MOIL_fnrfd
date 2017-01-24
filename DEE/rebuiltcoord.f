      subroutine rebuiltcoord()
c
c Subroutine to rebuilt the coordinates of the rotamer library after 
c DEE. This is a simplification of the builtcoord() and addrotlib() 
c used before. It adds coordinates of the rotamer library at the end c
c of the list of backbone atoms. The coordinates are obtained from the 
c internal coordinates rotamer library.
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
      
      character*12 name
      integer namel,level

      integer i,j,k,irot,indrot,typerot,BBindex(maxptmonrotint)
      integer BBmonoindex,iBB
c
      name = 'rebuiltcoord'
      namel = 12
c
      npt=poiatrotenh(0)
      nb=nbback
c
c
      do i=1,nposenh
c........find actual indeces of backbone atoms N, CA and C 
c........(BBindex(1,3,4), only ones needed to built atoms beyond CB.
         BBmonoindex=indexposenh(i)
         j=poipt(BBmonoindex-1)
         if (standardAA(BBmonoindex))then
            do k=1,4
               BBindex(k)=j+k
            end do
         else
            BBindex(1)=j+1
            BBindex(3)=j+2
            BBindex(4)=j+3
c...........put the CN (or ME1) instead of H. 
            BBindex(2)=j+5
         end if
c
c........find the monomer type and internal index, and then 
c........built the coordinates of the rotamer.
         do irot=poirotenhaux(i-1)+1,poirotenhaux(i)
            indrot=intindrotenh(rotaux(irot))
            typerot=typerotenh(rotaux(irot))
            iBB=i
            call re_addrotamer(BBindex,typerot,indrot,iBB)
c
            poiatrotenh(rotaux(irot))=npt
c
         end do
      end do
c
c              
      return
      end

c
c-------------------------------------------------------------------
c
      subroutine re_addrotamer(BBindex,typerot,indrot,iBB)
c
c Add the coordinates of a given rotamer at the end of the list of 
c atoms. The coordinates are obtained from the internal coordinates 
c rotamer library.
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      
      character*13 name
      integer namel,level

      integer i,j,k,indrot,typerot,BBindex(*),iBB,i1,i2,i3,i4
      integer jmax,kmax,globalrotindex,localind
      double precision req12,angleq123,chi1234
      logical found
c
      name = 're_addrotamer'
      namel = 13
      level = 1
c
c.....start building from the CB (which comes after the backbone atoms)
c.....BBindex(6)
      npt=npt+1
c
      do k=1,3
         coor(k,npt)=coorcb(k,iBB)
      end do
      BBindex(6)=npt
c
c.....add bond CA-CB
      nb=nb+1
      ib1(nb)=BBindex(3)
      ib2(nb)=BBindex(6)
c
      globalrotindex=poityrotint(typerot-1)+indrot
c
c.....set properties of the atom
      localind=poipts4(typerot-1)+6
      call set_properties(localind,typerot,iBB)
c     
c.....then the next atoms
      k=poiptrotint(globalrotindex-1)
      kmax=poiptrotint(globalrotindex)
c
c.....check for overflow in ib*s4 matrix
      if((nb+kmax-k-5).gt.maxbond)
     1     call alert(name,namel,'nb exceeded maxbond',17,level)
c
      do i=k+7,kmax
c
         npt=npt+1
         localind=localind+1
         call set_properties(localind,typerot,iBB)
c         
c........find the standard equilibrium distance between pt1,pt2
         j=poinbs4(typerot-1)+1
         jmax=poinbs4(typerot)
         found=.false.
         do while ((.not.found).and.(j.le.jmax))
            if((ib1s4(j).eq.ptrotint1(i) .and. ib2s4(j).eq.ptrotint2(i))
     1           .or. (ib1s4(j).eq.ptrotint2(i) .and. 
     1           ib2s4(j).eq.ptrotint1(i))) then
               req12 = reqs4(j) 
               found=.true.
            else
               j=j+1                
            end if
         end do
         if (j.gt.jmax)
     1        call alert(name,namel,'could not find bond',19,level)
c
c........find the  standard equilibrium angle between pt1,pt2,pt3
         j=poinangls4(typerot-1)+1
         jmax=poinangls4(typerot)
         found=.false.
         do while ((.not.found).and.(j.le.jmax))
            if((iangl1s4(j).eq.ptrotint1(i) .and. 
     1           iangl2s4(j).eq.ptrotint2(i) .and. 
     1           iangl3s4(j).eq.ptrotint3(i)) .or. 
     1           (iangl3s4(j).eq.ptrotint1(i) .and. 
     1           iangl2s4(j).eq.ptrotint2(i) .and. 
     1           iangl1s4(j).eq.ptrotint3(i))) then
               angleq123 = angleqs4(j)
               found=.true.
            else
               j=j+1                
            end if
         end do
         if (j.gt.jmax)
     1        call alert(name,namel,'could not find angle',20,level)
c
c........put the atom at the correct position.
         chi1234=chirotint(i)
         BBindex(ptrotint1(i))=npt
         i1=BBindex(ptrotint1(i))
         i2=BBindex(ptrotint2(i))
         i3=BBindex(ptrotint3(i))
         i4=BBindex(ptrotint4(i))
         call poschi(i1,i2,i3,i4,req12,angleq123,chi1234)
c
c........add to the bond matrix
         nb=nb+1
         ib1(nb)=BBindex(ptrotint2(i))
         ib2(nb)=BBindex(ptrotint1(i))
c
c
      end do
c
c              
      return
      end
c






