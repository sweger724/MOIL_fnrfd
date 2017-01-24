        subroutine init_wre(uwene)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer uwene
        
        write(uwene,100)
100     format(1x,///,' Parameters for energy calculation ')
        if (.not.ebyes) then
         write(uwene,101)
101      format(1x,' Bond energy will be skipped')
        end if
        if (.not.ethyes) then
         write(uwene,102)
102      format(1x,' Angle energy will be skipped')
        end if
        if (.not.etoyes) then
         write(uwene,103)
103      format(1x,' Torsion energy will be skipped')
        end if
        if (.not.eimyes) then
         write(uwene,104)
104      format(1x,' Improper torsion energy will be skipped')
        end if
        if (.not.evdyes) then
         write(uwene,105)
105      format(1x,' van der Waals energy will be skipped')
        end if
        if (.not.eelyes) then
         write(uwene,106)
106      format(1x,' Electrostatic energy will be skipped')
        end if
        write(uwene,107)
107     format(/)
        if (ctrue) then
         write(uwene,108)dsqrt(cutele2),dsqrt(cutvdw2)
108      format(1x,' Constant dielectric will be used.',
     1   ' elec. Cutoff=',f12.5,/,1x,' vdW cutoff ',f12.5)
        end if
        if (shift) then
         write(uwene,109)
109      format(1x,' Shift will be applied to electrostatic ')
        end if
        if (esymyes) then
         write(uwene,111)
111      format(1x,' Cubic symmetry is used to simulate'
     1    ,'periodic system')
        end if

        if (gbsabool) then
            if (gbobcbool) then
                write(uwene,113)
113             format(1x,' GB polarized solvation energy ',
     &                ' required (Onufriev)')
                write(uwene, 114) gbalpha, gbbeta, gbgamma
114             format(1x,' gbalpha=',f10.3,3x,'gbbeta=',f10.3,3x, 
     &                     'gbgamma=',f10.3) 
           
           else
                write(uwene,112)
112             format(1x,' GB polarized solvation energy', 
     &                ' required (Hawkins)')
           end if 
        end if

        if (gbnpbool) then
           write(uwene,115)
115        format(1x,' GB non-polarized solvation energy', 
     &            ' required (LCPO)')
           write(uwene, 116) surften
116        format(1x,' surface tension =',f8.4,' kcal/mol/A**2')
        end if

        write(uwene,117)
117     format(///)
        return
        end
