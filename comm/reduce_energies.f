        subroutine reduce_energies()
        implicit none
c
c this routine accumulate all the partial force vectors from
c all of the pe-s and add them up
c

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

c# terra include file
czva       include '../include/Tcommlib_f.h'
czva    include '../include/Trtlib_f.h'
c
c      include 'mpif.h'
c


        integer i
      
        if (energy_reduced) return
        energy_reduced = .true.
       
        call reduce_1(e_total)
        if (ebyes) call reduce_1(e_bond)
        if (ethyes) call reduce_1(e_theta)
        if (etoyes) call reduce_1(e_tors)
        if (eimyes) call reduce_1(e_imp)
        if (evdyes) then
                call reduce_1(e_vdw)
                call reduce_1(e_vdw14)
        end if
        if (eelyes) then
                call reduce_1(e_el)
                call reduce_1(e_el14)
        end if
        if (ewaldyes) then
                call reduce_1(e_corr)
                call reduce_1(e_ew_receip)
                call reduce_1(e_dir)
        endif
        if (ecnyes) then
                if (ncnst.eq.0) e_cnst = 0.d0
                call reduce_1(e_cnst)
        end if
        if (lcent)  call reduce_1(e_cent)
        if (esymyes) then
                call reduce_1(e_vsym)
                call reduce_1(e_lsym)
        end if
        do 1 i=1,nmb    
                if (emyes(i)) then
                        call reduce_1(e_morseb(i))
                end if
1       continue
        call reduce_1(tote_morsb)

        do 2 i=1,nmb
                if (repyes(i)) then
                        call reduce_1(e_repuls(i))
                end if
2       continue
        call reduce_1(tote_repuls)

        if (eteth_yes) call reduce_1(e_tether)
        if (ehyes) call reduce_1(e_hyd)

        return
        end



