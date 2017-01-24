      subroutine getenergies0()
c
c Subroutine to calculate the side-chain to backbone and side-chain to 
c side-chain energies. There is no need to calculate non-bonded lists 
c since the energy is calculated only once in the program. 
c The cut-off distances are determined per rotamer using the 
c geometric center of the side-chain residues. 
c
c Simplified version that returns 0 if Eiback or Eij < Eibackmax
c                                 1 if Eiback or Eij > Eibackmax
c 
c (Some of the subroutines are in getenergies.f)
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
         call getEiback0(rot_cent,epstmp,pick)
c
         call getEij0(rot_cent,epstmp,pick)
c
CDEB         call print_ener0()
c
      end if
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine getEiback0(rot_cent,epstmp,pick)
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
      double precision epstmp,pick,rot_cent(3,*),maxdist,ric(3),r2,ener
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
         Eiback(i)=e_vdw+e_el 
c
         if (abs(Eiback(i)).ge.Eibackmax) then
            Eiback(i)=Eibackmax*100.d0
         else
            Eiback(i)=0.d0
         end if
c
c
      end do
c
      return
      end
c
c
c-------------------------------------------------------------------
c
      subroutine getEij0(rot_cent,epstmp,pick)
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
            Eij(indexint)=e_vdw+e_el 
c
            if (abs(Eij(indexint)).ge.Eibackmax) then
               Eij(indexint)=Eibackmax*100.d0
            else
               Eij(indexint)=0.d0
            end if
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



