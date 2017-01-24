        program conn
        implicit none
c----------------------------------------------------------------
c LENGTH.BLOCK  - includes defintion for vector lengths
c CONNECT.BLOCK - includes the topology and energy term for the
c                       molecule
c PROPERT.BLOCK - includes definition for atom bond angle torsion
c                       and improper torsion properties
c MONOMERS.BLOCK- includes information on connectivity of individual
c                       monomers
c LINE.BLOCK    - used for line interpreter
c UNITS.BLOCK   - standard units used
c DEBUG.BLOCK   - To initiate debugging printout
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/PROPERT.BLOCK'
        include 'COMMON/MONOMERS.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/MUTA.BLOCK'
c----------------------------------------------------------------


c----------------------------------------------------------------
c of    - integer function return a unit number of file to be opened
c UNITS that we need
c upro  - unit number of properties file
c umono - unit number of connectivity file and particle description
c upoly - unit number with the polymer description (sequence + solvent)
c ubond - unit number with lisiting of bonds to be added to exsiting
c               molecules (can do without it).
c uedit - unit for editing the connectivity file
c uwcon - unit to write conn. and properties of the molecule (output)
        integer of
        integer upro,umono,upoly,ubond,uwcon,uedit

c       double precision getd
c----------------------------------------------------------------

c----------------------------------------------------------------
c local variables
c
        character*4 name
        logical find,fopen,hydro,shakm
        integer namel,i,level
        integer unit
        integer numon
        double precision getd
        data upro,umono,upoly,ubond,uedit,uwcon/6*-1/

c muta
        integer igroup,ipick(maxpt)

c hk add start
        integer flag_pair
c hk add end 

c--------------------------------------------------------------

c-------------------------------------------------------------------
c initialize file units for standard i/o
c Initialize some molecular properties (see CONNECT.BLOCK for details)
c and unit numbers (to -1)
c
        silent = .false.
        stdi = 5
        stdo = 6
        stderr = 0
        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        nmb    = 0
        nchg   = 0
        lestyp = 0
        eps    = 1.d0
        shakm  = .false.

c muta
        muta   =.FALSE.

c hk add start
        flag_pair = 0
c hk add end

c--------------------------------------------------------------------

        debug       = .false.
        hydro       = .false.
        arith       = .false.
        prll_on_off = .false.
        mdivyes     = .true.
c       mdivyes     = .false.
        hvdw0       = .true.
        name  = 'CONN'
        namel = 4

c----------------------------------------------------------------
c open junk file unit (25) and
c Input to connect is set to standard input (stdi)
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
        unit = stdi
c----------------------------------------------------------------
c 1 is loop on commands to conn BEFORE calling propert. etc
c intention: to open files
c
1       continue
        call rline(name,namel,unit)
        eps = getd('epsi',eps)
        !  force field detection
        if (find('mshk')) shakm       = .true.
        if (find('prll')) prll_on_off = .true.
        if (find('debu')) debug = .true.
        if (find('hydr')) hydro = .true.
        if (find('muta')) muta  = .true.        
        if (find('arit')) arith = .true.
        if (find('mdiv')) mdivyes  = .true.
        if (find('nomd')) mdivyes  = .false.
        if (find('hvdw')) hvdw0    = .false.
        if (find('file')) then
         if (find('prop')) upro  = of()
         if (find('mono')) umono = of()
         if (find('poly')) upoly = of()
         if (find('wcon')) uwcon = of()

c hk start
         if (find('wco2')) then
            flag_pair = 1
            uwcon = of()
         endif
c hk end


         if (find('ubon')) ubond = of()
         if (find('uedi')) uedit = of()
        end if
        if (find('acti')) then
         go to 2
        end if
        go to 1
2       continue
         
        eps  = 1.d0/eps

        if (debug) then
         write(stdo,*)' After opening files '
         write(stdo,*)' pro mono poly wcon '
         write(stdo,*)upro,umono,upoly,uwcon
        end if
c fopen - inquire if units are open
c
         if (.not.fopen(umono)) then
          level = 1
          call alert(name,namel,'umono not open',14,level)
         else if (.not.fopen(upro)) then
          level = 1
          call alert(name,namel,'upro not open',13,level)
         else if (.not.fopen(upoly)) then
          level = 1
          call alert(name,namel,'Upoly not open',14,level)
         else if (.not.fopen(uwcon)) then
          level = 1
          call alert(name,namel,'Uwcon not open',14,level)
         end if
c----------------------------------------------------------------

c----------------------------------------------------------------
c READ PARTICLE PROPERTIES: NAME, CHARGE, VDW, 
c                           BOND ANGLE TORSION & IMPROPER TORSIONS
c
c data base for atoms' properties
c data: atom_name atom_mass atom_charge atom_epsilon atom_sigma
c

c data base for bonds
c data: atom_name atom_name bond_force_constant bond_equilibrium_distance
c

c data base for angles
c data: atom_name atom_name atom_name angle_forc_cnst angle_equi_posit
c

c data base for torsions
c data: atom_name atom_name atom_name atom_name amplitude1
c                               amplitude2 amplitude3 period phase
c

c data base for improper torsions
c data: atom_name atom_name atom_name atom_name forc_cnst equi
c
        silent = .true.
        call  prprt(upro,basenm,basems,basechg,
     1  baseeps,basesgm,maxunq,totunq,
     2  bondp,basekb,basereq,maxubnd,totubnd,
     3  anglp,basekt,baseaeq,maxuagl,totuagl,
     4  torsp,basamp1,basamp2,basamp3,baseper,basephs1,
     5  basephs2,basephs3,
     5  maxutrs,totutrs,improp,baseimp,
     6  baseieq,maxuimp,totuimp,arith,el14,v14)
c----------------------------------------------------------------

      if (el14.lt.0 .or. v14.lt.0) then
       level = 1
       call alert(name,namel,'1-4 interactions not specified!',31,level)
      end if

c-------------------------------------------------------------------
c GETTING MONOMERS PROPERTIES
c uconn - unit file for defining monomers connectivity
c mono  - character*4 names of monomers
c nprt  - intg pointr to particle indices in the monomer
c               nprt(i) is the last particle of monomer i
c mono_pt - vectors of names of particles in a monomer mono_pt(nprt(k-1)+i)
c               is the unique atom name of the i-th particle in the
c               k-1 monomer
c nbnd    - integer pointer to bond indices in the monomer. nbnd(i)
c               is the last bond of monomer i
c basenm  - particle names according to property index
c basechg - particle charges according to property index
c monob1 monob2 - pointer to specific bonds monob[1-2](nbnd(k-1)+i) the
c               the i-th bond in the k-th monomer
c updown  - particle identifier
c               updown(nprt(k-1)+i) corresponds to an action to take
c               on the i-th particle of a the k-th glue monomer
c               updown = -1  -> particle from the monomer before
c               updown =  0  -> particle from the current monomer
c               updown =  1  -> particle from the next monomer
c               updown = 999 -> particle to be deleted from next monomer
c               updown =-999 -> particle to be deleted from previous monomer
c iden    - a pointer to particle index. iden(nprt(k-1)+i) is the property
c               index of the i-th particle in the k-th monomer
c maxmono - maximum number of monomers allowed (must be larger or equal
c               to the some of normal monomers (numno) and glue monomers
c               (numglu) found 
c totunq  - number of unique particles found (input)
c numon   - number of monomers found (output)
c 
        call monomers(umono,mono,nprt,mono_pt,
     1          nbnd,basenm,basechg,monob1,
     2          monob2,updown,iden,maxumon,
     3          totunq,numon,maxdiv,monodiv,mdivyes,m_charge,m_said)
        silent = .false.
        if (debug) write(stdo,*) ' numon = ',numon
c-------------------------------------------------------------------


c-------------------------------------------------------------------
c Build the current molecule
c
        call lego(upoly,mono,nprt,mono_pt,nbnd,
     1          monob1,monob2,updown,iden,numon,
     2          basenm,basechg,basesgm,baseeps,basems,totunq,
     4          bondp,basekb,basereq,m_charge,m_said)
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c if ubond exists, add bond to exisiting system via addbond
c
        if (ubond.gt.0) then
         call addbond(ubond)
        end if
c-------------------------------------------------------------------
c call angtor to make lists of angles torsions and improper
c torsions
        call angtor(shakm)

c-------------------------------------------------------------------
c if uedit exists, add connectivity file via edcon
c
        if (uedit.gt.0) then
         call edcon(uedit)
        end if

c-------------------------------------------------------------------
c create a sorted charge file for faster calculation of electrostatic
        call chsort()

c
c If hydrophobicity is requested call set_hydro
        if (hydro) call set_hydro()

c
c Count the number of water molecules
        call count_waters()

         do i=1,npt
          mutaid(i)=0
         enddo

         if (muta) then
          call rline(name,namel,unit)
          if (find('MUTA')) then
           call pick(ipick,igroup)
           call add_muta(ipick,igroup)
          else
           level=1
           call alert(name,namel,'expected muta not found',23,level)
          endif
         else
          do i=1,npt
           mutaid(i)=0
          enddo
         endif

c-------------------------------------------------------------------
c chen
        call add_les()
c-------------------------------------------------------------------
c Call wconn: write down current topolopgy and properies of molecular
c system
c

        call wconn(uwcon,v14,el14,flag_pair)

        write(stdo,*) 'prll=',prll_on_off

        write(stdo,*)'el14 is ',el14
        write(stdo,*)'v14 is ',v14
        do 3 i=1,NBULK
         write(stdo,100)BULK(i)
100      format(/,1x,'*** MOLECULE ',a4, ' ASSEMBLED ***'
     1     ,//)
3       continue
        stop
        end







