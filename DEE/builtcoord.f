      subroutine builtcoord()
c
c Subroutine to built the coordinates of the backbone and the rotamer 
c library. It uses the connectivities of the backbone and of the 
c different aminoacids, as well as the coordinates of the backbone and 
c the internall coordinates for each aminoacid. 
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      
      character*10 name
      integer namel,level
c
c
c
      integer nptold
      integer oldindexcb(maxnposenh)
      integer newindex(maxpt)
c
      name = 'builtcoord'
      namel = 10
c
      call removeSC(oldindexcb,newindex,nptold)
c
      call fixindexBB(oldindexcb,newindex,nptold)
c
      call addrotlib()
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine removeSC(oldindexcb,newindex,nptold)
c
c Subroutine to remove the original side chains from the residues 
c chosen for enhancement. Checks for bridges and N-methil residues
c (using atoms type ME1 and CN). Uses the fact that the first atoms 
c are the backbone atom and the last ones are the bridge or 
c N-methilated, when appropriate.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c      
      character*8 name
      integer namel,level

      integer i,j,k,ii,jj,nptold
      double precision coor_temp(3,maxpt)
      integer oldindexcb(*),newindex(*)
      integer poipt_new(0:maxmono),iposenh,imono
      character*4 ptnm_temp(7),ptnmcb
      data ptnm_temp/'N   ','H   ','CA  ','C   ','O   ','ME1 ','CN  '/
      data ptnmcb/'CB  '/
c
      name = 'removeSC'
      namel = 8
c
      do i=1,npt
         do j=1,3
            coor_temp(j,i)=coor(j,i)
         end do
         newindex(i)=0
      end do
c
      nptold=npt
      npt=0
      iposenh=0
c
c
      do i=1,totmon
         if (ipickm(i).eq.1)then
c
            iposenh=iposenh+1
c
c...........get position of the H.
            j=poipt(i-1)+2
            if (ptnm(j).eq.ptnm_temp(2))then
c..............standard AA, with H at the second position.Get only the 
c..............five atoms at the backbone. 
               standardAA(i)=.true.
               do j=1,5
                  jj=poipt(i-1)+j
                  npt=npt+1
                  newindex(jj)=npt
                  if(ptnm(jj).ne.ptnm_temp(j))
     &                 call alert(name,namel,'wrong atom name',15,1)
                  do k=1,3
                     coor(k,npt)=coor_temp(k,jj)
                  end do
               end do
c..............get the position of the CB
               jj=poipt(i-1)+6
               do k=1,3
                  coorcb(k,iposenh)=coor_temp(k,jj)
               end do
c..............get the new index of CA and the old index of CB
               oldindexcb(iposenh)=jj
            else
c..............non-standard AA, (N-alquilated ones)
               standardAA(i)=.false.
               j=1
               jj=poipt(i-1)+j
               npt=npt+1
               newindex(jj)=npt
               if(ptnm(jj).ne.ptnm_temp(j))
     &              call alert(name,namel,'wrong atom name',15,1)
               do k=1,3
                  coor(k,npt)=coor_temp(k,jj)
               end do
               do j=2,4
                  jj=poipt(i-1)+j
                  npt=npt+1
                  newindex(jj)=npt
                  if(ptnm(jj).ne.ptnm_temp(j+1))
     &                 call alert(name,namel,'wrong atom name',15,1)
                  do k=1,3
                     coor(k,npt)=coor_temp(k,jj)
                  end do
               end do
c..............get the position of the CB
               jj=poipt(i-1)+5
               do k=1,3
                  coorcb(k,iposenh)=coor_temp(k,jj)
               end do
c..............get the new index of CA and the old index of CB
               oldindexcb(iposenh)=jj
c..............get the rest of the molecule (N-methyl or Bridge)
               do while ((jj.le.poipt(i)).and.
     &              (ptnm(jj).ne.ptnm_temp(6)).and.
     &              (ptnm(jj).ne.ptnm_temp(7)))
                  jj=jj+1
               end do
               if (jj.gt.poipt(i))
     &                 call alert(name,namel,'could not find CN',17,1)
               do j=jj,poipt(i)
                  npt=npt+1
                  newindex(j)=npt
                  do k=1,3
                     coor(k,npt)=coor_temp(k,j)
                  end do
               end do
            end if
         else
c...........just copy            
            standardAA(i)=.true.
            do j=poipt(i-1)+1,poipt(i)
               npt=npt+1
               newindex(j)=npt
               do k=1,3
                  coor(k,npt)=coor_temp(k,j)
               end do
            end do
         end if
c
         poipt_new(i)=npt
c
      end do
c
c.....Add to the exclusion list of CB the atoms before it that have 
c.....this CB in their exclusion list. Do it before fixing the atomic 
c.....indeces.
      exc1s5(0)=0
      totexs5=0
      do i=1,nposenh
c
c........scan over the atoms from the previous monomer to CB-1
c.......(the first and the last monomers are assumed to never be picked)
         imono=indexposenh(i)
         do j=poipt(imono-2)+1,oldindexcb(i)-1
            do k=exc1(j-1)+1,exc1(j)
               if (exc2(k).eq.oldindexcb(i))then
                  totexs5=totexs5+1
                  exc2s5(totexs5)=j
               end if
            end do
         end do
c
c........add the atoms from the exclusion list of CB
         j=oldindexcb(i)
         do k=exc1(j-1)+1,exc1(j)
            totexs5=totexs5+1
            exc2s5(totexs5)=exc2(k)
         end do
c         
         exc1s5(i)=totexs5
c
      end do      
c
c
      do i=1,totmon     
         poipt(i)=poipt_new(i)
      end do
c     
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine fixindexBB(oldindexcb,newindex,nptold)
c
c Subroutine to fix the bonded indeces of the backbone and to find out 
c all atoms in the backbone that interact with the CB (including atoms 
c from neighbor residues)
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      
      character*10 name
      integer namel,level

      integer i,j,k,inew,nptold,nb_new,totex_new,imono
      integer totexs5_new,exc1_new(0:maxpt),exc1s5i,jnew
      integer oldindexcb(*),newindex(*)
c
      name = 'fixindexBB'
      namel = 10
c
c.....fix according to the new backbone the indeces: poimon, ptid,
c.....ptnm, ptms, ptchg, epsgm6, epsgm12, flagchr, poipt
c
      do i=1,nptold
         inew=newindex(i)
         if (inew.ne.0) then
            poimon(inew)=poimon(i)
            ptid(inew)=ptid(i)
            ptnm(inew)=ptnm(i)
            ptms(inew)=ptms(i)
            ptchg(inew)=ptchg(i)
            epsgm6(inew)= epsgm6(i)
            epsgm12(inew)=epsgm12(i)
            flagchr(inew)=flagchr(i)
         end if
      end do
c
c.....fix according to the new backbone the indeces: ib1 ib2 
      nb_new=0
      do i=1,nb
         if((newindex(ib1(i)).ne.0).and.(newindex(ib2(i)).ne.0))then
            nb_new=nb_new+1
            ib1(nb_new)=newindex(ib1(i))
            ib2(nb_new)=newindex(ib2(i))
         end if
      end do
c
      nb=nb_new      
c
c.....fix according to the new backbone the indeces: exc1,exc2
c.....maybe exc1 and exc2 are not necessary anymore.
      totex_new=0
      do i=1,nptold
         inew=newindex(i)
         if (inew.ne.0) then
            do j=exc1(i-1)+1,exc1(i)
               jnew=newindex(exc2(j))
               if (jnew.ne.0)then
                  totex_new=totex_new+1
                  exc2(totex_new)=jnew
               end if
            end do
            exc1_new(inew)=totex_new
         end if
      end do
c
      totex=totex_new
c
      do i=1,npt
         exc1(i)=exc1_new(i)
      end do
c
c
c.....fix according to the new backbone the indeces: exc1s5, exc2s5
      totexs5=0
      exc1s5i=0
      do i=1,nposenh
         do j=exc1s5i+1,exc1s5(i)
            jnew=newindex(exc2s5(j))
            if (jnew.ne.0) then
               totexs5=totexs5+1
               exc2s5(totexs5)=jnew
            end if
         end do
         exc1s5i=exc1s5(i)
         exc1s5(i)=totexs5
      end do
c
c.....keep the number of bonds in the backbone
      nbback=nb
c
c              
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine addrotlib()
c
c Add the coordinates of the rotamer library at the end of the list of 
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
      
      character*9 name
      integer namel,level

      integer i,j,k,irot,indrot,typerot,BBindex(maxptmonrotint)
      integer BBmonoindex,iBB
c
      name = 'addrotlib'
      namel = 9
c
c
      poiatrotenh(0)=npt
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
         do irot=poirotenh(i-1)+1,poirotenh(i)
            indrot=intindrotenh(irot)
            typerot=typerotenh(irot)
            iBB=i
            call addrotamer(BBindex,typerot,indrot,iBB)
c
            poiatrotenh(irot)=npt
c
c...........built the exclusion list for the atoms of this rotamer
c...........(it has only side-chain to backbone exclusions. 
            call addrotexclusion(BBindex,irot,iBB,typerot)
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
      subroutine addrotamer(BBindex,typerot,indrot,iBB)
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
      
      character*10 name
      integer namel,level

      integer i,j,k,indrot,typerot,BBindex(*),iBB,i1,i2,i3,i4
      integer jmax,kmax,globalrotindex,localind
      double precision req12,angleq123,chi1234
      logical found
c
      name = 'addrotamer'
      namel = 10
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
c-------------------------------------------------------------------
c
      subroutine set_properties(localind,typerot,iBB)
c
c Fill the properties of the atom from the internal coordinates 
c connectivity file.
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      
      character*14 name
      integer namel,level

      integer localind,typerot,iBB,iBB1,iBB0,newmono
      character*1 numberAA(0:9)
      data numberAA/'0','1','2','3','4','5','6','7','8','9'/
c
      name = 'set_properties'
      namel = 14
c
      newmono=totmon+iBB
      poimon(npt)=newmono
      ptid(npt)=ptids4(localind)
      ptnm(npt)=ptnms4(localind)
      ptms(npt)=ptmss4(localind)
      ptchg(npt)=ptchgs4(localind)
      epsgm6(npt)=epsgm6s4(localind)
      epsgm12(npt)=epsgm12s4(localind)
      flagchr(npt)=flagchrs4(localind)
c
      if (iBB.lt.10)then
         moname(newmono)='ROT'//numberAA(iBB)
      else
         iBB1=iBB/10
         iBB0=mod(iBB,10)
         moname(newmono)='RO'//numberAA(iBB1)//numberAA(iBB0)
      end if
c              
      return
      end
c
c
c-----------------------------------------------------------------------
c
      subroutine addrotexclusion(BBindex,irot,iBB,typerot)
c
c Subroutine to built the exclusion list for the atoms of rotamer irot,
c including only side-chain to backbone exclusions. 
c (CB=exc1s5, CG+=exc1s4)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      
      character*15 name
      integer namel,level
c
      integer i,ii,j,jj,irot,BBindex(*),iCB,typerot
      integer iBB,totexold
c
      name  = 'addrotexclusion'
      namel = 15
      level =1
c
c.....first for the CB
      iCB = poiatrotenh(irot-1)+1
      totexold = totex 
      totex = totex + exc1s5(iBB) - exc1s5(iBB-1)
c.....check for overflow in exc2 matrix
      if(totex.gt.maxex)
     1     call alert(name,namel,'totex exceeded maxex',20,level)
c
      exc1(iCB) = totex 
c
      jj=exc1s5(iBB-1)
      do j = totexold+1,totex
         jj=jj+1
         exc2(j)=exc2s5(jj)
      end do
c
c.....then the rest
      ii=poipts4(typerot-1)+6
      do i = iCB+1,npt
         ii=ii+1
         totexold = totex 
         totex = totex + exc1s4(ii) - exc1s4(ii-1)
         exc1(i) = totex 
         jj=exc1s4(ii-1)
         do j = totexold+1,totex
            jj=jj+1
            exc2(j)=BBindex(exc2s4(jj))
         end do
      end do
c              
      return
      end
c
c






