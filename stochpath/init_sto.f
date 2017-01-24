      subroutine init_sto
c
c.v1.0 (last changed 24/2/97)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/COMM_LOCAL.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/METAL.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/PVM_LOCAL.BLOCK'
      include 'COMMON/RESTART.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/VELOC.BLOCK'
c      
c.....UNITS
      stdi = 5
      stdo = 6
      stderr = 0
      uwcrd = 2
      uwpth =99
c
C ENERGY DEFAULT PARAMETER
cc
      call init_ef()
c
c.....CONNECT
      totmon = 0
      npt    = 0
      nb     = 0
      nangl  = 0
      ntors  = 0
      nimp   = 0
      lestyp = 0
      nshak  = 0
      nshakm = 0
      iorie  = 0
      ncnst  = 0
      nbulk  = 0
c.....DEBUG
      debug  = .false.
c
      nori   = .false.
      eqms   = .false.
      lcent  = .false.
      shakl  = .false.
      shakb  = .false.
      shakm  = .false.
      freeze = .false.
      boltz  = .true.
      nocut  = .false.
      sdyes  = .false.
      no_scaling = .false.
      hvdw0      = .true.
c.....DYNA
      nstep        = 1
      neqstep      = 1
      ninfo        = 1
      ncoor        = 0
      newv         = 0
      n_tether     = 0
      nrigi        = 10
      nvelo        = 0
      nlist        = 1
      nscalv       = 0
      ntemp        = 1
      start_dyna   = 1
      irand        = -23518284
      dt           = 1.0d-3
      tempi(1)     = 300.d0
      tempf(1)     = 300.d0
      tolcons      = 1.d-7
      epshak       = 1.d-7
      epshakv      = 1.d-7
      itershak     = 100
      fmax         = -1.d0
      symanneal= .false.
c.....ACTPARA
      crdstyl     = 'INIT'
      nstepopt    = 1
      npri        = 1
      nwcrd       = 500
      fixend      = .true.
      tolg        = 1.d-3
      estred      = 1.0d-4
      gamma       = 1.d0
      stoscale    = 1.d0
      stobeta     = 0.d0
      lap         = .false.
      nliststo    = 1
      smlp        = 1
      dtsto       = 1.0d0
      twodtsto    = 2.d0*dtsto
      bc1         = .true.
      bc2         = .false.
      damp        = .false.
      dtopt       = 1.0d-5
      tempanl     = 30.d0
      irandsto    = -3141
      ianlmax     = 1
      anlfac      = 0.9d0
      lanneal  = .false.
      allupdate = .true.
czva  no frozen by default (sepfast - obsolet)
      sepfast = .false.
      freeze = .false.
czva  (STO: slow_frz =.true.) and  (DYNA: slow_frz =.false.)
      slow_frz =.true.
c
      call inpstoline
c
czva  reset the list of STO to one (lsts=1). 
c
      nliststo=1
c
      return
      end
