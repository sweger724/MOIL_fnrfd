      subroutine getenergies()
c
c Subroutine to calculate the side-chain to backbone and side-chain to 
c side-chain energies. There is no need to calculate non-bonded lists 
c since the energy is calculated only once in the program. 
c The cut-off distances are determined per rotamer using the 
c geometric center of the side-chain residues. 
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
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'

      
      character*11 name
      integer namel,level
c
      double precision rot_cent(3,maxnmonenh),epstmp,pick
c
c
      name = 'getenergies'
      namel = 11
      level=1
c
      call init_energies(epstmp,pick)
c
      call getcentroids(rot_cent)
c
      if (evdyes.or.eelyes) then
c
         call getEiback(rot_cent,epstmp,pick)
c
         call getEij(rot_cent,epstmp,pick)
c
CDEB         call print_ener()
c
      end if
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine init_energies(epstmp,pick)
c
c Initialize the non-bonded interactions (as in the cdie subroutine). 
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'

c    
      character*13 name
      integer namel,level
c
      double precision pick,epstmp
      integer i,j
c
c
      name = 'init_energies'
      namel = 13
c
      if (eelyes) then
         epstmp = kofdie/eps
CDEB
CDEB         epstmp = 1.d0
      else
         epstmp = 0.0d0
      end if
      if (nocut) then
         cutvdw2 = 10000.d0
         cutele2 = 10000.d0
      end if
      
      if (evdyes) then
         pick = 1.0d0
      else
         pick = 0.0d0
      end if
c
      pointEij(0)=0
c     
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine getcentroids(rot_cent)
c
c Subroutine to calculate the geometric center of the different 
c rotamers.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
c    
      character*12 name
      integer namel,level
c
      double precision rot_cent(3,*),invnatm
      integer i,j,atbeg,atend,natm,iatm
c
c
      name = 'getcentroids'
      namel = 12
c
      do i=1,nmonenh
c
         atbeg=poiatrotenh(i-1)+1
         atend=poiatrotenh(i)
         natm=atend-atbeg+1
         invnatm=1.d0/dble(natm)
c
c........initialization
         do j=1,3
            rot_cent(j,i)=0
         end do
c
         do iatm=atbeg,atend
            do j=1,3
               rot_cent(j,i)=rot_cent(j,i)+coor(j,iatm)
            end do
         end do
c
         do j=1,3
            rot_cent(j,i)=rot_cent(j,i)*invnatm
         end do
c
      end do
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine getEiback(rot_cent,epstmp,pick)
c
c Subroutine to calculate the side-chain to backbone energies. No need 
c for non-bonded lists since the energy is calculated only once.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*9 name
      integer namel,level
c
      double precision epstmp,pick,rot_cent(3,*),maxdist,ric(3),r2
      double precision Eibacktemp
      integer i,j,k,atbeg,atend,exclmin,exclmax
      logical exclude
c
      name = 'getEiback'
      namel = 9
c
      maxdist=1.5d0*1.5d0*max(cutvdw2,cutele2)
c
      do i=1,nmonenh
c
c........first (CB) and last atom of the monomer.
         atbeg=poiatrotenh(i-1)+1
         atend=poiatrotenh(i)
c........smallest and largest BB atom index in the exclusion list of 
c........the monomer (from the CB exc. list).
         exclmin=npt+999
         exclmax=0
         do j=exc1(atbeg-1)+1,exc1(atbeg)
            exclmin=min(exclmin,exc2(j))
            exclmax=max(exclmax,exc2(j))
         end do             
c
c........initialization
         Eiback(i)=0
         e_vdw = 0.d0
         e_el  = 0.d0
c
         do j=1,poiatrotenh(0)
c
c...........check if it is necessary to check the exclusion list.
            if ((j.gt.exclmax).or.(j.lt.exclmin)) then
               exclude=.false.
            else
               exclude=.true.
            end if
c
c...........get centroid-atom distance
            do k=1,3
               ric(k)=rot_cent(k,i)-coor(k,j)
            end do
            r2 = ric(1)*ric(1)+ric(2)*ric(2)+ric(3)*ric(3)
c
            if (r2.lt.maxdist)then
               call energy_atom_side(j,atbeg,atend,exclude,
     &              epstmp,pick)
CDEB
CDEB            else
CDEB               write(*,*)'atom',j,'far away from rotamer',i
CDEB
            end if
c
         end do
c
         Eibacktemp=e_vdw+e_el
         if (abs(Eibacktemp).gt.hugeE) then
            Eiback(i)=hugeE
         else
            Eiback(i)=Eibacktemp
         end if
c
      end do
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine energy_atom_side(atback,atbeg,atend,exclude,
     &              epstmp,pick)
c
c Subroutine to calculate the atom-atom interactions.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*16 name
      integer namel,level
c
      double precision smallestdist
      parameter (smallestdist=1.d-6)
c
      double precision epstmp,pick,a,b,q,aback,bback,qback,
     &     rij(3),r2,s2,s,e2,s6,e1
      integer i,j,jmax,k,atback,atbeg,atend
      logical exclude,found,list1,list2,list3,charged
c
c
      name = 'energy_atom_side'
      namel= 16
      aback  = epsgm12(atback)*pick
      bback  = epsgm6(atback)*pick
      qback  = ptchg(atback)*epstmp
c
      do i=atbeg,atend
c
c........initialize the lists
         list1=.false.
         list2=.false.
         list3=.false.
         charged=(flagchr(i).and.flagchr(atback))
c........check the exclusion list for atback
         if (exclude) then
            found=.false.
            j=exc1(i-1)+1
            jmax=exc1(i)
            do while ((.not.found).and.(j.le.jmax))
               if(exc2(j).eq.atback)then
                  found=.true.
               else
                  j=j+1
               end if
            end do
         end if
         if (.not.found) exclude=.false.
c
         if (.not.exclude) then
c
c...........get atom-atom distance
            do k=1,3
               rij(k)=coor(k,i)-coor(k,atback)
            end do
            r2 = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
c
c...........avoid NAN's.
            if (r2.lt.smallestdist) then
               e_vdw=hugeE
               e_el =hugeE
            else
c
c..............make now a decision to which of the three lists this pair
c..............belongs:list1 includes both van der Waals and 
c..............          electrostatic and the max. distance is cutvdw2.
c..............        list2 uses the same cutoff (cutvdw2) but intends
c..............          to uncharged particle only.
c..............        list3 is for for distances larger than cutvdw2 
c..............          and smaller than cutele2 and includes only 
c..............          electrostatic.
c
c..............NOTE THAT cutvdw2 AND cutele2 INSTEAD OF cut-big2 AS IN
c..............THE NON-BONDED LIST.
c
               if (r2.lt.cutele2) then
c
c.................check for hydrogens. Hydrogens may have zero van-der 
c.................Waals radius if arith=.false. (default) and should be 
c.................placed in the third (electrostatic only) list 
c.................(if charged, otherwise exclude)
                  if (((epsgm12(i).lt.1.d-4).or.
     &                 (epsgm12(atback).lt.1.d-4)).and. charged) then
c
                     list3=.true.
c
                  else if (r2.le.cutvdw2) then
c
                     if (charged) then
c
                        list1=.true.
c
                     else
c
                        list2=.true.
c
                     end if
c
                  else if ((r2.gt.cutvdw2).and. charged) then
c
                     list3=.true.
c
                  end if
c
               end if
c
            end if 
c
         end if
c
c
c........continue the calculations. In case exclude=.true., rij was 
c........not calculated but it is OK since  all lists are .false.
c
         if (list1) then
c
            s2=1.0d0/r2
            s = dsqrt(s2)
            q = qback*ptchg(i)
            a = aback*epsgm12(i)
            b = bback*epsgm6(i)
            if (ctrue) then
               e2 = q*s
            else
               e2 = q*s2
            end if              
            s6=s2*s2*s2
            a = a*s6*s6
            b = b*s6
            e1 = a - b
CDEB
CDEB            e1 = pick
CDEB            e2 = epstmp
CDEB
            e_vdw = e_vdw + e1 
            e_el = e_el + e2
c
         else if (list2) then
c
            s2=1.0d0/r2
            a = aback*epsgm12(i)
            b = bback*epsgm6(i)
            s6=s2*s2*s2
            a = a*s6*s6
            b = b*s6
            e1 = a - b
CDEB
CDEB            e1 = pick
CDEB            e2 = epstmp
CDEB
            e_vdw = e_vdw + e1 
c
         else if (list3) then
c
            s2=1.0d0/r2
            s = dsqrt(s2)
            q = qback*ptchg(i)
            if (ctrue) then
               e2 = q*s
            else
               e2 = q*s2
            end if              
            s6=s2*s2*s2
CDEB
CDEB            e1 = pick
CDEB            e2 = epstmp
CDEB
            e_el = e_el + e2
c
         end if            
c            
      end do
c
c
      return
      end
c
c
c-------------------------------------------------------------------
c
      subroutine getEij(rot_cent,epstmp,pick)
c
c Subroutine to calculate the backbone to backbone energies. No need 
c for non-bonded lists since the energy is calculated only once.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*6 name
      integer namel,level
c
      double precision epstmp,pick,rot_cent(3,*),rij(3),r2,maxdist
      double precision Eij1
      integer i,j,k,jinit,atbegi,atendi,atbegj,atendj,indexint
c
      name = 'getEij'
      namel = 6
      indexint=0
      pointEij(0)=0
      maxdist=1.5d0*1.5d0*max(cutvdw2,cutele2)
c
      do i=1,nmonenh
c
c........first (CB) and last atom of the monomer.
         atbegi=poiatrotenh(i-1)+1
         atendi=poiatrotenh(i)
c
c........first rotamer of the next position
         jinit=poirotenh( posrotenh(i) )+1
c
         do j=jinit,nmonenh
c
            indexint=indexint+1
c
c...........initialization
            e_vdw = 0.d0
            e_el  = 0.d0
c
c...........first (CB) and last atom of the monomer.
            atbegj=poiatrotenh(j-1)+1
            atendj=poiatrotenh(j)
c
c
c...........get centroid-centroid distance
            do k=1,3
               rij(k)=rot_cent(k,i)-rot_cent(k,j)
            end do
            r2 = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
c
            if (r2.lt.maxdist)then
               call energy_side_side(atbegi,atendi,atbegj,atendj,
     &              epstmp,pick)
            end if
c
            Eij1=e_vdw+e_el 
c
c...........mark if i and j have a clash 
            if (abs(Eij1).gt.Eijmax) then
               DEPij(indexint)=.true.
               if (abs(Eij1).gt.hugeE) Eij1=hugeE
            else
               DEPij(indexint)=.false.
            end if
c
            Eij(indexint)=Eij1
c
         end do
c
         pointEij(i)=indexint
c
c
      end do
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine energy_side_side(atbegi,atendi,atbegj,atendj,
     &              epstmp,pick)
c
c Subroutine to calculate the atom-atom interactions.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*16 name
      integer namel,level
c
      double precision smallestdist
      parameter (smallestdist=1.d-6)
c
      double precision epstmp,pick,a,b,q,ai,bi,qi,
     &     rij(3),r2,s2,s,e2,s6,e1
      integer i,j,k,atbegi,atendi,atbegj,atendj
      logical found,list1,list2,list3,charged
c
c
      name = 'energy_side_side'
      namel= 16
c
      do i=atbegi,atendi
c
         ai  = epsgm12(i)*pick
         bi  = epsgm6(i)*pick
         qi  = ptchg(i)*epstmp
c
         do j=atbegj,atendj
c
c...........initialize the lists
            list1=.false.
            list2=.false.
            list3=.false.
            charged=(flagchr(i).and.flagchr(j))
c
c...........get atom-atom distance
            do k=1,3
               rij(k)=coor(k,i)-coor(k,j)
            end do
            r2 = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
c
c...........avoid NAN's.
            if (r2.lt.smallestdist) then
               e_vdw=hugeE
               e_el =hugeE
            else
c
c..............make now a decision to which of the three lists this pair
c..............belongs:list1 includes both van der Waals and 
c..............          electrostatic and the max. distance is cutvdw2.
c..............        list2 uses the same cutoff (cutvdw2) but intends
c..............          to uncharged particle only.
c..............        list3 is for for distances larger than cutvdw2 
c..............          and smaller than cutele2 and includes only 
c..............          electrostatic.
c
c..............NOTE THAT cutvdw2 AND cutele2 INSTEAD OF cut-big2 AS IN
c..............THE NON-BONDED LIST.
c
c
c
               if (r2.lt.cutele2) then
c
c.................check for hydrogens. Hydrogens may have zero van-der 
c.................Waals radius if arith=.false. (default) and should be 
c.................placed in the third (electrostatic only) list 
c.................(if charged, otherwise exclude)
                  if (((epsgm12(i).lt.1.d-4).or.(epsgm12(j).lt.1.d-4))
     &                 .and. charged) then
c
                     list3=.true.
c
                  else if (r2.le.cutvdw2) then
c
                     if (charged) then
c
                        list1=.true.
c
                     else
c
                        list2=.true.
c
                     end if
c
                  else if ((r2.gt.cutvdw2).and. charged) then
c
                     list3=.true.
c
                  end if
c
               end if
c
            end if
c
c..........continue the calculations. In case exclude=.true., rij was 
c..........not calculated but it is OK since  all lists are .false.
c
            if (list1) then
c
               s2=1.0d0/r2
               s = dsqrt(s2)
               q = qi*ptchg(j)
               a = ai*epsgm12(j)
               b = bi*epsgm6(j)
               if (ctrue) then
                  e2 = q*s
               else
                  e2 = q*s2
               end if              
               s6=s2*s2*s2
               a = a*s6*s6
               b = b*s6
               e1 = a - b
CDEB
CDEB               e1 = pick
CDEB               e2 = epstmp
CDEB
               e_vdw = e_vdw + e1 
               e_el = e_el + e2
c
            else if (list2) then
c
               s2=1.0d0/r2
               a = ai*epsgm12(j)
               b = bi*epsgm6(j)
               s6=s2*s2*s2
               a = a*s6*s6
               b = b*s6
               e1 = a - b
CDEB
CDEB               e1 = pick
CDEB               e2 = epstmp
CDEB
               e_vdw = e_vdw + e1 
c
            else if (list3) then
c
               s2=1.0d0/r2
               s = dsqrt(s2)
               q = qi*ptchg(j)
               if (ctrue) then
                  e2 = q*s
               else
                  e2 = q*s2
               end if              
CDEB
CDEB               e1 = pick
CDEB               e2 = epstmp
CDEB
               s6=s2*s2*s2
               e_el = e_el + e2
c
            end if            
c    
         end do
c
      end do
c
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine print_ener()
c
c Subroutine to print the energies of the different rotamers for 
c debugging purposes.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*10 name
      integer namel,level
c
      double precision zero(maxnmonenh)
      integer i,j,jinit1,ind1,ind2,nrotam1
c
      name = 'print_ener'
      namel = 10
c
      write (6,*) 'initial residues'
      write (6,50)(i,i=1,nmonenh)
 50   format(8(1x,i9))
c
      write (6,*)
      write (6,*) 'side-chain backbone energy'
c
      write (6,100)( Eiback(i),i=1,nmonenh)
 100  format(8(1x,f9.5))
c     
      do i=1,nmonenh
         zero(i)=0.d0
      end do
c     
      write (6,*)
      write (6,*) 'side-chain side-chain  energy'
c
      do i=1,nmonenh
c
c........first rotamer of the next position
         jinit1=poirotenh( posrotenh(i) )
         nrotam1=nmonenh-jinit1
         ind1=pointEij(i-1)
         ind2=pointEij(i)
c
         if ((ind2-ind1).ne.nrotam1)
     &        call alert(name,namel,'something wrong',15,level)
c
         write (6,*)
         write (6,*)' rotamer', i
c
         write(6,100)(zero(j),j=1,jinit1),
     &        (Eij(j),j=ind1+1,ind1+nrotam1)
c     
c
      end do
c
      return
      end
c



