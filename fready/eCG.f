      subroutine eCG()
        implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/FREADY.BLOCK'

        integer i,l

        if (.not. Fix_2nd_structure) then

           call eCG_theta()
           call eCG_bond()

        endif
        
        call eCG_tors()
        call eCG_backboneHB()
        if (.not. efyes)  call eCG_LJ()

      e_CG = e_vdw + e_tors + e_theta + e_bond + e_el + e_el14

        return
        end
