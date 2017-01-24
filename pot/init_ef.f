        subroutine init_ef()
        implicit none
c
c subrotuine to initialize some reasonable values for energy/force flags
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/FREADY.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
c inserted for xcenter etc
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/STEER.BLOCK'
        integer n

        ctrue      = .true.
        shift      = .false.
        eenmyes    = .false.
        ebyes      = .true.
        ethyes     = .true.
        etoyes     = .true.
        eimyes     = .true.
        evdyes     = .true.
        eelyes     = .true.
        ecnyes     = .false.
        esymyes    = .false.
        eteth_yes  = .false.
        efyes      = .false.
        ehyes      = .false.
        eballyes   = .false.
        e14el_yes  = .true.
        e14v_yes   = .true.
        eCGyes     = .false.
        e_Kumbyes  = .false.
        qssbp      = .false.
        hvdw0      = .true.
        test_d     = .false.
        ewaldyes   = .false.
        vp_flag    = .false.
        metalyes   = .false.
        gbsabool   = .false.
        sgbbool    = .false.
        lcent      = .false.
        esteeryes  = .false.
c 
        xcenter = .true.
        ycenter = .true.
        zcenter = .true.
c
        lastViolated = .true.
        point1(0)  = 0
        point2(0)  = 0
        point3(0)  = 0

        kappa=0.10395d0
        dielectric=78.5d0

C       FREADY parameters
        CG_cutoff = 13.5d0
        HB_cutoff = 8.d0
        smooth_hardcore = .false.
        Fix_2nd_structure = .false.

C       ELASTIC NETWORK parameters
        Enm_cutoff = 7.d0
        enm_alpha = 0.d0
        enm_beta  = 1.d1
        enm_gamma = 1.d0


c repulsive wall
          rwall = .FALSE.
          nwalls = 0
          normw(1) = 0
          normw(2) = 0
          normw(3) = 0
          normw(4) = 0
          normw(5) = 0
          normw(6) = 0
          w0(1) = 0.d0
          w0(2) = 0.d0
          w0(3) = 0.d0
          w0(4) = 0.d0
          w0(5) = 0.d0
          w0(6) = 0.d0
c units Kcal/mol/A^6
          weps = 10.d0

c electric field
          efbox = 0.d0
          efield = 0.d0
          nefield = 3
          DVtype = 2
          DV = 0.d0
          efield_yes = .FALSE.

c scaling waters far away
          scale_away = .false.

c-----------------------------
      do 1 n=1,maxmorsb
        emyes(n)  = .false.
        repyes(n) = .false.
1      continue
      emyes0 = .false.
      specl  = .false.
      repyes0    = .false.
c-----------------------------

c update rate of  GB calculation
        gbsu = 0

        eps  = 1.d0
        cutvdw2 = 6.d0
        cutele2 = 8.d0
        cutvbig2 = -1.d0
        cutebig2 = -1.d0
        rmax     = -1.d0
        cutmono2 = -1.d0
c force fconstants for steer
        kxy = 0.d0
        kz  = 0.d0
        return
        end
