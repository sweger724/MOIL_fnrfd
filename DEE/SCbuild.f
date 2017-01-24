      program SCbuild
c
c  Create a rotamer library with several amino-acids, each at several 
c  rotamer conformations, and finds the best sequence for a given 
c  backbone. Based on the DEE code of Ora Schueler version 3.
c
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
      include 'COMMON/CCRD.BLOCK'
c
      integer uwcrd,uwcon
      
      character*7 name
      integer namel,level
c
      integer niterations
      logical printall
c           
      name = 'SCbuild'
      namel = 9

c.....initiate parameters
      call init_SCpredict()

c.....read parameters from input file (incl conn and CRD of constant 
c.....part)
      call input_SCbuild(uwcrd,uwcon)
      call flush(stdo)
c
c.....build the rotamer libraries according to the internal coordinates 
c.....general rotamer library
      call builtcoord()
c
c
      write(stdo,*) 'finished with the coordinate building'
      call flush(stdo)
c
c.....calculate the energy for side-chain to backbone and side-chain 
c.....to side-chain. 
      call getenergies()
      write(stdo,*) 'finished with the energy calculation'
      call flush(stdo)
c
c.....print out the rotamers combinations (left) with the respective 
c.....total energies
      printall=.false.
      call print_totener(printall)
c
c.....write the connectivity in terms of bonded pairs
      call buildmolec()
c
c.....write the coordinates
      call putcrd(uwcrd,'CHARM')
c
      end


