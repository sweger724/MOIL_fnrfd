      subroutine rrotalib(urlib)
c
c Reads the internal coordinates of a series of rotamers of various 
c amino-acids. The list of aminoacids need not to be in the same order 
c as in the connectivity file. The order of the later is used as 
c indexing. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      integer urlib
      logical notendoffile

      integer i,j,k,kk
      integer imono
c
      character*8 name
      integer namel,level
c
      character*4 rotname1,atomname1
      integer ind1,ind2,ind3,natomsrot,ntyrotint,ii,jj,nfirstatrot
      integer imono1
c
      name = 'rrotalib'
      namel = 8
c
      nptsrotint=0
      nrotint=0
      ntyrotint=0
c
      poiptrotint(0) = 0
      poityrotint(0) = 0
c      
      notendoffile=.true.
c
c read comment line
      read(urlib,*)
c
      do while (notendoffile)
c
c read monomer name and check the indeces.
         read(urlib,1000,err=99,end=999)rotname1,ind1,ind2,ind3
 1000    format(a4,3(1x,i5))
c
         nrotint=nrotint+1
         if (ind3.ne.nrotint) 
     &        call alert(name,namel,'wrong global index',18,1)
c
c........check if the type changed
         if (ind1.ne.ntyrotint) then
c
            if ( (ind1-1).ne.ntyrotint )
     &           call alert(name,namel,'wrong residue index',19,1)
            ntyrotint=ntyrotint+1
c
c...........find the WCON index of the monomer
            imono=1 
            do while ((monames4(imono).ne.rotname1).and.
     &           (imono.le.totmons4))
               imono=imono+1
            end do
            if (imono.gt.totmons4) 
     &           call alert(name,namel,'unknown residue type',20,1)
c
            natomsrot=ptnumbers4( poipts4(imono) )
            nfirstatrot= poipts4(imono-1)
         end if
c
c
         typerotint(nrotint)=imono
c
         k=poiptrotint(nrotint-1)
         poiptrotint(nrotint)=k+natomsrot
c         
         do i=1,natomsrot
            ii=i+k
            jj=i+nfirstatrot
            glbindrotint(ii)=nrotint
            read(urlib,1001,err=99,end=99)ptrotint1(ii),atomname1,
     &           ptrotint2(ii),ptrotint3(ii),ptrotint4(ii),chirotint(ii)
            if (atomname1.ne.ptnms4(jj)) then
               call alert(name,namel,'wrong order in the atoms',24,1)
            end if
c YP
c            chirotint(ii)=chirotint(ii)*pi180
            chirotint(ii)=chirotint(ii)/pi180
         end do
 1001    format(i3,1x,a4,3(1x,i5),1x,f7.2)
      end do
c
 99   continue
      call alert(name,namel,'problems with reading',19,1)
c
c
 999  continue
cout      rewind urlib
c
c.....get the number of rotamer types (using the monomer indeces of 
c.....the connectivity). Assumes (and checks) that the rotamers appear 
c.....in order (even if some monomers are missing)
      do i=1,totmons4
         poityrotint(i)=0
      end do
c
      imono1=0
      do i=1,nrotint
         imono=typerotint(i)
         if (imono.gt.imono1) then
            do j=imono1,imono-1
               poityrotint(j)=i-1
            end do
            imono1=imono
         else if (imono.lt.imono1) then
            call alert(name,namel,'residues out of order',21,1)
         end if
      end do
      do j=imono,totmons4
         poityrotint(j)=nrotint
      end do
c
      do i=1,totmons4
         nrotinttype(i)=poityrotint(i)-poityrotint(i-1)
      end do
c
      return
      end




