      subroutine rrotaconn(urota)
c
c Modification of rconn so that it reads the connectivity of a series 
c aminoacids and organize them in a long array but keeps the intrinsic
c numbering. It uses a set of arrays different from the ones in the 
c CONNECT_SHORT BLOCK (*s4 arrays). Keeps 
c only the relevant information about the sequence. It also changes 
c the exclusion-list of the atoms of the side chains to include only 
c the atoms of the backbone.
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      
      integer urota
      logical notendoffile

      integer imono,i,j,k,kk,k1,iold,l,junk(30)
      integer npt1,nb1,nangl1,totex1,totspe1,ntors1,nimp1
c
      character*4 moname1
      character*9 name
      integer namel,level

      name = 'rrotaconn'
      namel = 9

      imono = 0
      totmons4 = 0   
      npts4= 0
      nbs4= 0
      nangls4= 0
      totexs4= 0
c
      poipts4(0) = 0
      poinbs4(0) = 0
      poinangls4(0) = 0
      poiexs4(0) = 0
c
      level = 1
c      
      notendoffile=.true.
c
      do while (notendoffile)
c
c........monomer name
         read(urota,1000,err=99,end=999)moname1
 1000    format(a4)
         imono= imono+1 
         totmons4 = imono
         monames4(imono)=moname1
c
         read(urota,*)
         read(urota,*)
         read(urota,*)
c........reads: npt nb nangl totex totspe irot totchi nchg
         read(urota,1003)npt1,nb1,nangl1,ntors1,nimp1,totex1,totspe1
 1003    format(7x,7i6)
         npts4 = npts4+npt1
         poipts4(imono) = npts4
         nbs4 = nbs4+nb1
         poinbs4(imono) = nbs4 
         nangls4 = nangls4+nangl1
         poinangls4(imono) = nangls4
         totexs4 = totexs4+totex1
         poiexs4(imono) = totexs4 
c
         do i=1,9
            read(urota,*)
         end do
c
c........Properties of particles list : 
c........pt ptid lesid  ptnm   ptms   ptchg   epsgm6 epsgm12
c
         k = poipts4(imono-1)
         do 10 i=k+1, npts4
            read(urota,1004) ptnumbers4(i),ptids4(i),ptnms4(i),ptmss4(i)
     1           ,ptchgs4(i),epsgm6s4(i),epsgm12s4(i)
c
            poimons4(i) = imono
 10      continue
 1004    format(1x,i5,7x,i3,5x,a4,1x,f7.2,2x,f9.5,1x,e12.5,e12.5)
c
c
c........Bonds list: 
c........ib1 ib2 kbond req 
         read(urota,*)
         read(urota,*)
         kk = poinbs4(imono-1)+1
         do 20 i= kk, nbs4
            read(urota,1006)ib1s4(i),ib2s4(i),reqs4(i)
cout            ib1s4(i) = ib1s4(i) + k 
cout            ib2s4(i) = ib2s4(i) + k
 20      continue
 1006    format(2(1x,i5),12x,f10.4)
c
c
c........Angles list:
c........iangl1 iangl2 iangl3 kangl angleq
         read(urota,*)
         read(urota,*)
         kk = poinangls4(imono-1)+1
         do 30 i= kk, nangls4
            read(urota,1008)iangl1s4(i),iangl2s4(i),iangl3s4(i),
     1           angleqs4(i)
cout            iangl1s4(i) =iangl1s4(i) + k  
cout            iangl2s4(i) =iangl2s4(i) + k  
cout            iangl3s4(i) =iangl3s4(i) + k  
c YP
c            angleqs4(i) =angleqs4(i)* pi180
            angleqs4(i) =angleqs4(i)/pi180
 30      continue
 1008    format(3(1x,i5),12x,f10.5)
c
c........skip torsions and improper torsions
         do i=1,ntors1+nimp1+6
            read(urota,*)
         end do
c
c........Exclusion list 1-2 1-3 1-4
c........atom number, number of exclusions and list
         kk = poiexs4(imono-1)
         k = poipts4(imono-1)
         i = poipts4(imono-1)
         do while (kk.lt.poiexs4(imono))
            read(urota,117)j,k1
 117        format(2(1x,i5))
            j = j+k
            i = i+1
            if (i.lt.j) then
               iold=i
               do i=iold,j-1
                  exc1s4(i) = kk
               end do
            else if (i.gt.j) then
               call alert(name,namel,'Fishy exclusion list!',21,level)
            end if
c            
            exc1s4(j) = kk + k1
            read(urota,1018)(exc2s4(l),l=kk+1,exc1s4(j))
cout            do l=kk+1,exc1s4(j)
cout               exc2s4(l) = exc2s4(l) + k   
cout            end do
            kk = kk + k1
         end do
         do j=i+1,poipts4(imono)
            exc1s4(j) = exc1s4(j-1)
         end do
 1018    format(1x,10i5)
c
c........skip special list 1-4
         do i=1, totspe1+3
            read(urota,*)
         end do
c
c........Charge flags. T = charged atom , F = uncharged')
         k  = poipts4(imono-1)
         read(urota,1017)(flagchrs4(i),i=k+1,poipts4(imono))
 1017    format(1x,15l2)
         read(urota,*)
c........LES copy labels (dump into junk array)
         read(urota,1019,err=99)(junk(i),i=1,npt1)
 1019    format(1x,15i4)
c
      end do
c
 99   continue
      call alert(name,namel,'problems with reading',19,1)
 999  continue
c
c
c.....deconvolute the exclusion list
      call deconvolutes4()
c
      return
      end
c
c
c
      subroutine deconvolutes4()
c
c Deconvolutes the exclusion list such that it includes only the  
c sidechain atoms exclusions with respect to the backbone atoms.
c It is done to save computation in the calculation of the Side-chain
c Backbone energy. 
c Later in the program the exclusion list of CB is recalculated to 
c include the atoms in neighbor monomers. PAY ATTENTION TO THE 
c DIFFERENT INDECING OF THE ATOMS IN THE LISTS exc1 and exc2
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
c
      integer exc1s4_temp(0:maxpttype),exc2s4_temp(maxextype)
      integer totexs4_temp,imono,i,j,k,k0,k1,j0,j1,iatom
c
c
c.....backup the matrices
      totexs4_temp=totexs4
      do i=1,totexs4
         exc2s4_temp(i)=exc2s4(i)
      end do
      do i=0,npts4
         exc1s4_temp(i)=exc1s4(i)
      end do
c
c.....reinitialize the exclusion matrix
      totexs4=0
c
      do imono=1,totmons4
c
         k0=poipts4(imono-1)
         k1=poipts4(imono)
c
c........remove the exclusion of the backbone atoms
         do i=k0+1,k0+5
            exc1s4(i)=totexs4
         end do
c
c........loop over the side-chain atoms.
         do iatom=k0+6,k1
c
            do i=k0+1,k0+5
               j0=exc1s4_temp(i-1)+1
               j1=exc1s4_temp(i)
               do j=j0,j1
                  if ((exc2s4_temp(j)+k0).eq.iatom) then
                     totexs4=totexs4+1
                     exc2s4(totexs4)=i-k0
                  end if
               end do
            end do
c
            exc1s4(iatom)=totexs4   
c
         end do
c
         poiexs4(imono)=totexs4
c
      end do
c
      return 
      end


      




