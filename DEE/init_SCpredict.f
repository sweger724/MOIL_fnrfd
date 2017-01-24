      subroutine init_SCpredict()

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
c
      integer iseed
      integer i
c
c open junk file for rline
c
      jnkf=25
      open(unit=jnkf,status='scratch')
c      call terra_open_scratch(jnkf)
c
      debug = .false.
      stdi = 5
      stdo = 6
      stderr = 0
c
c.....connectivity for backbone
      totmon = 0
      npt    = 0
      nb     = 0
      totex  = 0
      totspe = 0
      poipt(0) = 0
      exc1(0)= 0
c
c.....energy
      ctrue      = .true.
      shift      = .false.
      evdyes     = .true.
      eelyes     = .true.
c
      eps  = 1.d0
      cutvdw2 = 6.d0
      cutele2 = 8.d0
      cutvbig2 = 0.d0
      cutebig2 = 0.d0
      rmax     = -1.d0
cout      cutEir = 999.d0
cout      maxcutEir = 999.d0
c
      Eibackmax=10000.d0
      Eijmax   =10000.d0
c
      return
      end


