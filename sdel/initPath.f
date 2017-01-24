       subroutine initPath()
       implicit none
c
c calculate path using SDEL algorithm
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SDEL.BLOCK'
        include 'COMMON/SGB.BLOCK'

        integer i,l

C  compute momentum of the first stucture
      if (first) then
        p0_initial = 2 * (ENERGYC-e0(1))
        if (p0_initial.le.0) then
           write (6,*) 'Initial momentum is < 0',e_total
           stop
        end if
        p0_initial = sqrt(p0_initial)
      endif 


C  compute momentum of the last stucture
      if (last) then
        p0_final = 2 * (ENERGYC-e0(pseg+2))
        if (p0_final.le.0) then
           write (6,*) 'Final momentum is < 0',e_total
           stop
        end if
        p0_final = sqrt(p0_final)
      endif

      return
      end
