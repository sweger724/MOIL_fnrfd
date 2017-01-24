      subroutine init_d2v
c
c (last modified in 16/1/96 by RO)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      
c.....UNITS
      stdi = 5
      stdo = 6
      stderr = 0
c.....CONNECT
      totmon = 0
      npt    = 0
      nb     = 0
      nmb    = 0
      nangl  = 0
      ntors  = 0
      nimp   = 0
      lestyp = 0
      nbulk  = 0
c.....DEBUG
      debug  = .false.
c.....ENERGY
      nocut  = .true.
      lcent = .false.
c
c open scratch file for line manipulation (check if still used).
c
c.....LINE
      jnkf = 25
      open (unit=jnkf,status='scratch')
c

      return
      end

