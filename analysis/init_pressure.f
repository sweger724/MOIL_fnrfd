       subroutine init_pressure() 

       implicit none

c common blocks

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/LD.BLOCK'
        include 'COMMON/PT.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PRESS.BLOCK'

c declare the variables you need
c loops
        integer i,imol
        logical check

c initialization

        do i=1,maxpt
         patom(i) = 0
         pmol(i) = 0
        enddo

c create the list of molecules 
c       call dmolecule(npt,nb,ib1,ib2,ptnm,moname,poimon,nmol,patom,pmol)
        call create_mol_BFS(npt,nb,ib1,ib2,ptnm,moname,poimon,
     1                      nmol,patom,pmol,ordLISTatm)
        write (*,*) 'nmol=',nmol
        do imol=1,nmol
         write (*,*) 'patom(',imol,')=',patom(imol)
        enddo
        check =.false.
        do i=1,npt
         if (ordLISTatm(i).ne.i) check=.true.
        enddo
        if (check) then
        do i=1,npt
         write (*,*) 'ordLISTatm(',i,')=',ordLISTatm(i)
        enddo
        endif
        do i=1,npt
         write (*,*) 'pmol(',i,')=',pmol(i)
        enddo
c create the molecular masses array
        call compute_massmol(npt,ptms,patom,massmol,ordLISTatm)
        do imol=1,nmol
         write (*,*) 'massmol(',imol,')=',massmol(imol)
        enddo
      return 
      end


         subroutine create_mol_BFS(npt,nb,ib1,ib2,ptnm,moname,poimon,
     1                             nmol,patom,pmol,ordLISTatm)
         implicit none
        
c input 
         integer npt,nb,ib1(*),ib2(*),poimon(*)
         character*4 moname(*),ptnm(*)
c output
         integer nmol,patom(*),pmol(*),ordLISTatm(*)
c indexes
         integer i,j,ipt,iatm,imol
c layers and lists 
         integer Nlist,Nlay,Nnewlay,Nord
         integer list(npt),layer(npt),new_layer(npt)

c this subroutine creates the molecules through the Breadth-First Search algorithm
c the search goes in layers:
c first layer: the first atom of the molecule
c  - search on the bonded list all the atoms bound to the first atom
c second layer: all the bonded atoms to the first atom of the molecule
c  - search on the bonded list all the atoms bound to those in the first layer
c ...
c all the atoms found in this way belong to a molecule.
c the list of atoms is then ordered through a binary three algorithm.
c EXCEPTION: water molecules (SPCE or TIP3P) are create according to their name

         ipt = 1
         imol = 0
         do i=1,npt
         layer(i) = 0 
         new_layer(i) = 0
         list(i) = 0
         enddo
         Nlist = 0
         Nord = 0
         do i=1,npt
          pmol(i) = 0
         enddo

         DO WHILE (ipt.lt.npt)

c build exeption for water due to missing bonds 
c due to mshk use
           if  (moname(poimon(ipt)).eq.'TIP3' .or.
     1          moname(poimon(ipt)).eq.'SPCE') then
            Nlist = 3
            list(1) = ipt
            list(2) = ipt + 1
            list(3) = ipt + 2
            goto 201
           endif

         if (layer(1).eq.0) then
          layer(1) = ipt
          Nlay = 1
          list(1) = ipt
          Nlist = 1
         endif

c create a new layer (checking that the particles in the layer are not 
c                     already in the list)
200      continue
         !write (*,*) 'list',Nlist,(list(i),i=1,10)
         call create_layer(ib1,ib2,nb,layer,new_layer,Nlay,Nnewlay,
     1                     list,Nlist)
!         write (*,*) 'layer',Nnewlay,new_layer(1),
!     1   new_layer(2),new_layer(3)
c add to the list 
         if (Nnewlay.gt.0) then
ccccc update layer
          do i=1,Nnewlay
           list(Nlist+i) = new_layer(i)
           layer(i) = new_layer(i)
           new_layer(i) = 0
          enddo
          Nlist = Nlist + Nnewlay
          Nlay = Nnewlay
          Nnewlay = 0
          goto 200
         endif
         do i=1,Nlist
          !write (*,*) i,list(i)
         enddo
         !stop
c if you are here the molecule is completed
c order the number of particles belonging to the molecule using a binary tree
         call binary_tree(list,Nlist) 
201      continue
         imol = imol + 1
         do i=1,Nlist
          ordLISTatm(i+Nord) = list(i)
          pmol(list(i)) = imol
         enddo
         Nord = Nord+Nlist
         patom(imol) = Nord
c go to the next particle, which is the one with the smallest particle index
c that was not yet assigned to a molecule
         do i=ipt+1,npt
          if (pmol(i).eq.0) goto 202
         enddo
202      continue
         ipt = i
         layer(1) = 0
        ENDDO

        nmol = imol

        return
        end

       subroutine create_layer(ib1,ib2,nb,layer,new_layer,Nlay,Nnewlay,
     1 list,Nlist)
       implicit none
       integer ilay,ilist,ib
       integer ib1(*),ib2(*),nb,layer(*),new_layer(*),Nlay,Nnewlay
       integer list(*),Nlist

       Nnewlay = 0
       do ilay=1,Nlay
        do ib=1,nb
         if (ib1(ib).eq.layer(ilay)) then
c check if the particle is already in the list
          do ilist=1,Nlist
!           write (*,*)'inin1',list(ilist),ib2(ib),list(ilist).eq.ib2(ib)
           if (list(ilist).eq.ib2(ib)) goto 300
          enddo
          do ilist=1,Nlay
           if (layer(ilist).eq.ib2(ib)) goto 300
          enddo
          do ilist=1,Nnewlay
           if (new_layer(ilist).eq.ib2(ib)) goto 300
          enddo
           Nnewlay = Nnewlay + 1
           new_layer(Nnewlay) = ib2(ib)
         else if (ib2(ib).eq.layer(ilay)) then
c check if the particle is already in the list
          do ilist=1,Nlist
!           write (*,*)'inin2',list(ilist),ib1(ib),list(ilist).eq.ib1(ib)
           if (list(ilist).eq.ib1(ib)) goto 300
          enddo
          do ilist=1,Nlay
           if (layer(ilist).eq.ib1(ib)) goto 300
          enddo
          do ilist=1,Nnewlay
           if (new_layer(ilist).eq.ib1(ib)) goto 300
          enddo
           Nnewlay = Nnewlay + 1
           new_layer(Nnewlay) = ib1(ib)
         endif
300      continue
!         write (*,*) 'in',ib,Nnewlay,new_layer(Nnewlay)
        enddo
       enddo

       return
       end

       subroutine binary_tree(list,Nlist)
       implicit none
       integer i,itree,ilist,iord
       integer list(*),Nlist
       integer order(Nlist)
       integer parent(Nlist),lchild(Nlist),rchild(Nlist)

       do i=1,Nlist
        parent(i) = 0
        lchild(i) = 0
        rchild(i) = 0
       enddo

c first element of the list is the root
c the other elements follow...

       do i=2,Nlist 
        itree = 1 
400     continue
        if (list(i).lt.list(itree)) then
         if (lchild(itree).gt.0) then
          itree = lchild(itree)
          goto 400
         else
          parent(i) = itree
          lchild(itree) = i
          goto 401
         endif
        else if (list(i).gt.list(itree)) then
         if (rchild(itree).gt.0) then
          itree = rchild(itree)
          goto 400
         else
          parent(i) = itree
          rchild(itree) = i
          goto 401
         endif
        endif
401     continue
        enddo
c read the tree ...
        do i=1,Nlist
         !write (*,*) 'tree',list(i),parent(i),lchild(i),rchild(i)
        enddo
        iord = 0
        itree = 1
402     continue
        do while (lchild(itree).gt.0)
         itree = lchild(itree)
        enddo
        iord = iord + 1
        order(iord) = list(itree)
        !write (*,*) 'ord',iord,order(iord)
        if (parent(itree).gt.0) then
         lchild(parent(itree)) = rchild(itree)
         parent(rchild(itree)) = parent(itree)
         itree = parent(itree)
        else
         if (rchild(itree).gt.0) then
          itree = rchild(itree)
          parent(itree) = 0
         else 
          goto 403
         endif
        endif
        goto 402

403     continue
        do i=1,Nlist
         list(i) = order(i)
        enddo
        return
        end

         subroutine dmolecule(npt,nb,ib1,ib2,ptnm,moname,poimon,
     1                        nmol,patom,pmol)

       implicit none

c      Program takes the bond information and classify
c      the particles as molecules 
c      the atoms covalently connected to each other defines th emolecule
c      patom(i) - is the last atom of the i th molecule 
c      hwat - integer flag to explicitly take into account that H1 and H2
c             of water molecules are bound to the oxygen.

         integer nb,ib1(*),ib2(*),ir(npt)
         integer i,npt,imol,nmol
         integer patom(*),poimon(*),pmol(*)
         integer hwat
 
c new method
         integer ifirst,ilarge,j
         logical b1test1,b1test2,b2test1,b2test2
 
         character*4 moname(*),ptnm(*)

c ******* new method ********

          ifirst = 0
          ilarge = 0
          imol = 0
          DO WHILE (ifirst .lt. npt)
           ifirst = ilarge + 1 
c build exeption for water due to missing bonds 
c due to mshk use
           if  (moname(poimon(ifirst)).eq.'TIP3' .or.
     1          moname(poimon(ifirst)).eq.'SPCE') then
            ilarge = ifirst + 2
            goto 100
           endif
C look for large bond of first particle
           ilarge = ifirst
           DO i=1,nb
            if (ib1(i).eq.ifirst) then
             if (ib2(i).gt.ilarge) ilarge = ib2(i)
            else if (ib2(i).eq.ifirst) then
             if (ib1(i).gt.ilarge) ilarge = ib1(i)
            endif
           ENDDO
c if ilarge eq ifirst the monomer is over
           if (ilarge .eq. ifirst) goto 100
c check the bonds of ilarge. if they are larger than ilarge that is the new ilarge
90         continue
           DO i=1,nb
            if (ib1(i).eq.ilarge) then
             if (ib2(i).gt.ilarge) then 
              ilarge = ib2(i)
              goto 90
             endif
            else if (ib2(i).eq.ifirst) then
             if (ib1(i).gt.ilarge) then
              ilarge = ib1(i)
              goto 90
             endif
            endif
           ENDDO
c if you are here it means that you did not find anything larger than the
c large ilarge found. Let's look to all the bonds belonging to this connected
c segment and see if any other particle with larger index are around
           DO i=1,nb
            b1test1 = ib1(i).ge.ifirst
            b1test2 = ib1(i).le.ilarge
            if (b1test1 .and. b1test2) then
             if (ib2(i).gt.ilarge) then
              ilarge = ib2(i)
              goto 90
             endif
            endif
            b2test1 = ib2(i).ge.ifirst
            b2test2 = ib2(i).le.ilarge
            if (b2test1 .and. b2test2) then
             if (ib1(i).gt.ilarge) then
              ilarge = ib1(i)
              goto 90
             endif
            endif
          ENDDO
100       continue
          imol = imol + 1
c if you reached this point the molecule is completed! create the lists and
c go back at the beginning
          do i=ifirst,ilarge
           pmol(i) = imol
          enddo
          patom(imol) = ilarge
          ifirst = ilarge
          ENDDO

          nmol = imol

          return
          end


        subroutine compute_massmol(nmol,ptms,patom,massmol,ordLISTatm)

        implicit none

        integer iatm,imol,k,nmol,patom(*),ordLISTatm(*)
        double precision ptms(*),massmol(*)

        k = 1
        DO imol=1,nmol
         massmol(imol) = 0.d0
         DO iatm=k,patom(imol)
          massmol(imol) = massmol(imol) + ptms(ordLISTatm(iatm))
         ENDDO
         k = patom(imol)+1
        ENDDO

        return
        end

