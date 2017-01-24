        subroutine init_con()
c
c initialization of some con and units variables
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/DYNA.BLOCK'

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

        return
        end
