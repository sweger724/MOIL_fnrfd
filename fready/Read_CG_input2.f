       subroutine Read_CG_input2()

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/FREADY.BLOCK'
        include 'COMMON/ELASTIC.BLOCK'

        integer uCGpar, i, of, urtet,lc
        double precision getd
        logical find
        character*6 getchar
        character*5 crdtyp

c Parameter for fixing secondary structure in
C the coarse grained model
C         if (find('fix2')) Fix_2nd_structure = .TRUE.


c Plastic network parameters
        enm_cutoff =getd('cutE',Enm_cutoff)
        CG_cutoff = getd('cutC',CG_cutoff)
        HB_cutoff = getd('cutH',HB_cutoff)
        enm_alpha = getd('alFh',enm_alpha)
        enm_beta  = getd('bEta',enm_beta)
        enm_gamma = getd('gamE',enm_gamma) 
        
        if (find('ENM2')) then 
           eenmyes = .TRUE.

c energy default parameters
c
           nocut  = .true.
           ebyes  = .false.
           ethyes = .false.
           etoyes = .false.
           eimyes = .false.
           evdyes = .false.
           eelyes = .false.
           e14el_yes  = .false.
           e14v_yes   = .false.
           eCGyes = .false.
        end if

        if (find('DECO')) smooth_hardcore = .true.
        
      return

      end
