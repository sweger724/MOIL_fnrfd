      subroutine print_totener(printall)
c
c Subroutine to print the energy of all combinations with the 
c respective rotamer index.
c debugging purposes.
c 
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c    
      character*13 name
      integer namel,level
c
      double precision etotal,etotalmin,detotal
      integer i,j,jinit1,ind1,gloop(maxnposenh),
     &     gloopinit(maxnposenh),gloopend(maxnposenh)
      logical initpos,printall
c
      name = 'print_totener'
      namel = 13
c
      etotalmin=1.d30
      detotal=0.5d0
cout      detotal=5.0d0
c
      do i=1,nposenh
         gloopinit(i) = poirotenhaux(i-1)+1
         gloopend(i)  = poirotenhaux(i)
         gloop(i)     = gloopinit(i)
      end do
c
      write (6,*) 
      write (6,*) 'PRINTING ENERGIES OF THE CONFORMATIONS'
      write (6,*) 
c
      initpos=.true.
      do while(initpos)
c
         etotal=0.d0
c
         do i=1,nposenh
            etotal=etotal+Eiback(rotaux(gloop(i)))
         end do
c
         do i=1,nposenh-1
            do j=i+1,nposenh
               ind1=pointEij(rotaux(gloop(i))-1)+rotaux(gloop(j))-
     &              poirotenh(i)
               etotal=etotal+Eij(ind1)
            end do
         end do
c
         if (printall) then
            if (etotal.le.Eibackmax) 
     &           write(6,*)(intindrotenh(rotaux(gloop(i))),i=1,nposenh),
     &           etotal
         else
c...........print only the lowest energies.
            if (etotal.le.(etotalmin+detotal)) then
               write(6,*)(intindrotenh(rotaux(gloop(i))),i=1,nposenh),
     &              etotal
               if (etotal.le.etotalmin) etotalmin=etotal
            end if
         end if
c
         call newloop(gloopinit,gloopend,gloop,initpos)
c
      end do


      return
      end
c
c
c-----------------------------------------------------------------------
c
c
      subroutine newloop(gloopinit,gloopend,gloop,initpos)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c
      logical initpos
      integer i,gloop(maxnposenh),gloopinit(maxnposenh),
     &     gloopend(maxnposenh)
c
      character*7 name
      integer namel,level
c
      name = 'newloop'
      namel = 7
c
c
      do 10 i=1,nposenh
         gloop(i) = gloop(i) + 1
         if (gloop(i).le.gloopend(i)) then
            return
         else
            if (i.lt.nposenh) then
               gloop(i) = gloopinit(i)
            else
               initpos=.false.
               return
            end if
         end if
 10   continue
c
      return
      end
c     
c
