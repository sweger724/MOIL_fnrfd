      subroutine input_SCpredict(uwcrd,uwcon)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
      include 'COMMON/CCRD.BLOCK'
c
      integer geti,of
      double precision getd
      logical find,fopen
      character*15 name
      integer namel,level
c
      character*4 monatpos(maxnmontyp)
      integer ires,imono,npos,nmon,ntype,i,j,k,ii,addnpt,finalnpt
c
      integer urcrd,uback,urlib,urota,uwcrd,uwcon
      data urcrd,uback,urlib,urota /4*99/
c
      name = 'input_SCpredict'
      namel = 15
c
CDEB
CDEB      open (stdi,file='SC_predict.inp',status='old')
c
      uwcrd=99
      uwcon=99
      nposenh=0
      nmonenh=0      
      poirotenh(0)=0
      poiatrotenh(0)=0
      do i=1,totmon
         ipickm(i) = 0
      end do
c
      level = 1
      nocut = .false.
      ctrue = .true. 
c
 1    continue
c
         call rline(name,namel,stdi)
         if (find('debu')) debug = .true.
c
         if (find('file')) then
c
c...........get connectivity of the backbone
            if (find('back')) then
               uback  = of()
               call rbackconn(uback)
            end if
c
c...........get coordinate (using charm format)
            if (find('rcrd')) then
               urcrd = of()
               if (.not.fopen(uback)) then
                  call alert(name,namel,'uback not opened',16,level)
               end if
               call getcrd(urcrd,'CHARM')
            end if
c
c...........get connectivity for the monomers used in the rotamer 
c...........library 
            if (find('rota')) then
               urota = of()
               call rrotaconn(urota) 
            end if
c
c...........get rotational library side chains in internal coordinates.
            if (find('rlib')) then
               urlib = of()
               if ((.not.fopen(urota))) then
                  call alert(name,namel,'urota not opened',17,level)
               end if
               call rrotalib(urlib) 
            end if
c
            if (find('wcrd')) uwcrd = of()
            if (find('wcon')) uwcon = of()
         end if
c
c
c---------------------------------------------
c........get residues and rotamers involved in the rotamer library.
         if (find('ROTA')) then
c
c...........initial checks
            if (.not.fopen(urlib)) then
               call alert(name,namel,'urlib not opened',17,level)
            end if
c
            finalnpt=npt
            nposenh = geti('#var',0) 
            if(nposenh.gt.maxnposenh) then
               call alert(name,namel,'too many positions to enhance',
     &              29,level)
            end if               
c
c...........read input lines
            do i=1,nposenh
               call rline(name,namel,stdi)
               ires = geti('#res',0)
               ipickm(ires) = 1
               indexposenh(i) = ires
               nmon = geti('#mon',0)
               read(stdi,1000)(monatpos(j),j=1,nmon)
 1000          format(15(a4,1x))
               write(stdo,1000)(monatpos(j),j=1,nmon)
c
               do j=1,nmon
c
c.................find monomer type
                  imono=1 
                  do while ((monames4(imono).ne.monatpos(j)).and.
     &                 (imono.le.totmons4))
                     imono=imono+1
                  end do
                  if (imono.gt.totmons4) then
                     if(monatpos(j).eq.'    ')then
                        call alert(name,namel,'missing residue type',
     &                       20,level)
                     else
                        call alert(name,namel,'unknown residue type',
     &                       20,level)
                     end if
                  end if
c
c.................get number of rotamers of the given type
                  ntype=nrotinttype(imono)
                  if (ntype.eq.0) then
                     call alert(name,namel,
     &                    'monomer not in the inter. coord. library',
     &                    40,level)
                  end if
c
c.................check beforehand if matrix indeces will overflow
                  if(nmonenh+ntype.gt.maxnmonenh)then
                     call alert(name,namel,'too many monomers',
     &                    17,level)
                  end if               
c
                  do ii=1,ntype
                     nmonenh=nmonenh+1
                     typerotenh(nmonenh)=imono
                     posrotenh(nmonenh)=i
                     intindrotenh(nmonenh)=ii
                  end do
c
                  addnpt=poipts4(imono)-poipts4(imono-1)-5
                  finalnpt=finalnpt+addnpt
c
               end do
c
               poirotenh(i)=nmonenh
               if ((poirotenh(i)-poirotenh(i-1)).gt.maxnmonenhppos) then
                  call alert(name,namel,
     &                 'too many rotamers in a single position',
     &                 35,level)
               end if

c
            end do
c
c...........check beforehand if matrix indeces will overflow
            if(finalnpt.gt.maxpt)then
               call alert(name,namel,
     &              'npt with rotamers larger than maxpt',35,level)
            end if
c
         end if
c
c---------------------------------------------
c........read energy parameters
         cutvdw2  = getd('rvmx',cutvdw2)
         cutele2  = getd('relx',cutele2)
         rmax     = getd('rmax',rmax)
         eps      = getd('epsi',eps)
c     
         if (find('cdie')) ctrue = .true.
         if (find('nbfi')) efyes = .true.
         if (find('rdie')) ctrue = .false.
c
         if (evdyes .and. find('novd')) evdyes  = .false.
         if (eelyes .and. find('noel')) eelyes  = .false.         
c........not used
         if (nocut  .and. find('cute')) nocut   = .false.
c
c........DEE parameters
         Eibackmax = getd('eibm',Eibackmax)
         Eijmax    = getd('eijm',Eijmax)
c
         if (find('acti')) goto 2
c
         goto 1
 2    continue
c     
c.....check that required files were opened
      if (.not.fopen(urcrd)) then
         call alert(name,namel,'urcrd not opened',16,level)
      else if (.not.fopen(uwcrd)) then
         call alert(name,namel,'uwcrd not opened',16,level)
      else if ((.not.fopen(uwcon))) then
         call alert(name,namel,'uwcon not opened',16,level)
      end if
c
c.....rmax is maintained here for old input (with a single cutoff 
c.....to work).
      if (rmax.gt.0) then
         cutvdw2 = rmax
         cutele2 = rmax
      end if
c
c.....NOT AS IN THE DYNAMICS CALCULATION, the present calculation do 
c.....not use cutvbig2 and cutebig2, since we need only single 
c.....energy calculations. The non-bonded list is also skipped.
c
c
c
      write(stdo,100)uback,urcrd,(cutvdw2),(cutele2)
 100  format(1x,'Energy calculations ',/,
     3     ' Backbone connectivity file to be read from unit ',i5,/,
     4     ' Initial coordinates to be read from unit ',i5,/,
     6     ' The cutoff distance for van der Waal is ',f10.5,/,
     7     ' The cutoff distance for electrostatic is ',f10.5)
      if (ctrue) then
         write(stdo,*)' Constant dielectric constant is used'
      else
         write(stdo,*)' distance dependent dielectric constant is used'
      end if
c
c      
      cutvdw2  = cutvdw2*cutvdw2
      cutele2  = cutele2*cutele2
c
      if(nposenh.le.0) then
         call alert(name,namel,'no positions to enhance',23,level)
      end if   
c
c.....initialize DEE matrices.
c
      do i=1,nmonenh
         kept(i)=.true.
         rotaux(i)=i
      end do
c
      do i=0,nposenh
         poirotenhaux(i)=poirotenh(i)
      end do
c
      nmonenhleft=nmonenh
c
      hugeE=10000*dmax1(Eijmax,Eibackmax)
c      
c            
      write(stdo,*)' clash energy', Eibackmax
c
c.....close files that already finished 
      close(uback)
      close(urota)
      close(urlib)
      close(urcrd)
c
c
      return
      end 











