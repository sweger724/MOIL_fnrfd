      subroutine eforce()
c     
c     a branch routine to accumulate the different energies and forces
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/METAL.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/EBALL.BLOCK'

      integer namel,i,level,n
      double precision e_tmp
      character*6 name

        integer j

      data name/'eforce'/
      data namel/6/

C@
c
c       do i=1,npt
c               dpot(1,i) = 1.d0
c               dpot(2,i) = 1.d0
c               dpot(3,i) = 1.d0
c       end do
c       return
C@
      e_total = 0.d0
      energy_reduced = .false.
c     NOTE - for fast modes required for OM it is assumed that only
c     bonds and angles (waters and counterions) as well as 
c     non-bonded interaction are present

c     jm & zva 
c     e_* should be initialized also here for the sto !!!!
c     double initialization is not so BAD (zva)
c     
      e_bond= 0.d0
      e_theta= 0.d0
      e_tors= 0.d0
      e_imp= 0.d0
      e_cnst= 0.d0
      e_el14=0.0d0 
      e_vdw14=0.0d0
c     
      e_vdw = 0.d0
      e_el  = 0.d0
      e_corr = 0.d0
      e_dir =0.d0 
czva
      e_vsym = 0.d0
      e_lsym = 0.d0   
c     jm
      e_ball = 0.d0
c re

      do 3 i=1,npt
         dpot(1,i) = 0.d0
         dpot(2,i) = 0.d0
         dpot(3,i) = 0.d0
 3    continue

cLuca
        virvdw1 = 0.d0  ! watwat vdw virial
        virvdw2 = 0.d0  ! watwat symm vdw vir
        virvdw3 = 0.d0  ! solute vdw vir
        virvdw4 = 0.d0  ! solute symm vdw vir
        virvdw5 = 0.d0  ! 1-4 vdw vir
        vircent = 0.d0  ! cent pot vir, zeroed again in ecent
        vircon  = 0.d0  ! constr vir


      if (gbsabool) then
       sndbool= 0
       e_gbsa=0.0d0

      call egb_f(npt,ptchg,dielectric,
     $  sndbool,0,0,debugsgb)
      e_total=e_total + e_gbsa
      endif

c     vp
c     fix positions of virtual prts before energy calculations
      if (vp_flag) call vp_locate()
c     end of vp


      if (ebyes  .and. nb.gt.0) then
         call ebond()
         e_total = e_total + e_bond
      end if
      if (ethyes .and. nangl.gt.0) then
         call etheta()
         e_total = e_total + e_theta
      end if
      if (etoyes .and. ntors.gt.0) then
         call etors()
         e_total = e_total + e_tors
      end if
      if (eimyes .and. nimp .gt.0) then
         call eimphi()
         e_total = e_total + e_imp
      end if
      if (ecnyes .and. ncnst .gt.0) then
         call cnstrn()
         e_total = e_total + e_cnst
      end if

      if (evdyes.or.eelyes) then
         if(test_d) call test_dis()
         if (ctrue) then
            if (shift) then
               level = 1
               call alert(name,namel,'Shift not supported',19,level)
            else if (.not.efyes) then
                  if (ewaldyes) then
                     call cdie_ewald()
                     if (.not.metalyes) then 
                       if (my_pe .eq. num_pes-1) then 
                         call ewald_receip()
                       else
                         e_ew_receip = 0.d0
                       end if
                     end if
                     e_total = e_total + e_ew_receip + e_corr + e_self
                  else
                     if (.not. 3*nwaters.eq.npt) then
                        call cdie()
                     end if
                  end if
            else if (efyes) then
               call nbfinit()
            end if
         else
            call rdie()
         end if          

         if (nwaters.gt.1) then
            e_tmp = e_vdw + e_el
            if (ewaldyes) then
               call watwat_ewald()
            else
               call watwat()
            end if
         end if
         e_total = e_total + e_vdw + e_el
         e_dir = e_dir + e_el
      end if


      if ((evdyes.or.eelyes) .and. totspe.gt.0 .and. (.not.efyes)) then
         if (e14el_yes.or.e14v_yes) then
            call ener14()
            e_total = e_total + e_el14 + e_vdw14
         endif
      end if
      if (metalyes) then
         e_lctd = 0.d0
         e_wall = 0.d0
         call metal_wall()
         e_total = e_total + e_wall + e_lctd
      end if

c     only cdie is supported for symmetry related operations
c     
      if (esymyes) then
         e_vsym = 0.d0
         e_lsym = 0.d0
         if (metalyes) e_lmet = 0.d0

         if (ewaldyes) then

            if (iblock1(symnum).gt.0 .or.
     1           iblock2(symnum).gt.0 .or. iblock3(symnum).gt.0)
     2           call symcdie_ewald()
            if (iblckwt1(symnum).gt.0 .or. iblckwt2(symnum).gt.0)
     1           call symwat_ewald()

         else

            if (iblock1(symnum).gt.0 .or.
     1           iblock2(symnum).gt.0 .or. iblock3(symnum).gt.0)
     2           call symcdie()
            if (iblckwt1(symnum).gt.0 .or. iblckwt2(symnum).gt.0)
     1           then
               e_tmp = e_vsym + e_lsym
               call symwat()
            end if
         end if
c     @ RE - moved into the loop.
         e_total = e_total + e_vsym + e_lsym
         if (metalyes) e_total = e_total + e_lmet
         e_dir = e_dir + e_lsym

      end if

      if (lcent.and.icenter.gt.0) then
         call ecent()
         e_total = e_total + e_cent
      end if
c      if (eballyes.and.nwaters.gt.0) then
       if (eballyes) then
         call eball()
         e_total = e_total + e_ball
      end if

c     ----------------------------------
c     repuls - exponential repuslion between atoms
c     special- Landau- Zener model for curve crossing
c     emyes  - Morse potential energy
c     
      if (nmb.gt.0) then
         if (debug) write(*,*) 'number of morse bonds ',nmb
         do 990 n=1,nmb
            if ( repyes(n)) then
               emyes(n) = .false.
               call repuls(n)
               e_total = e_total + e_repuls(n)
            else if (emyes(n)) then
               call emorse(n)
               e_total = e_total + e_morseb(n)
            end if
            if (debug) write(*,*)'Morse energy = ',e_morseb(n)
 990     continue
         if (specl) then
            call special()
         end if
      end if
c     ---------------------------------

      if (eteth_yes) then
         call etether()
         e_total = e_total + e_tether
      end if

      if (ehyes) then
         call ehydro()
         e_total = e_total + e_hyd
      end if

      if (eenmyes) then
        call eenm()
        e_total = e_total + e_enm
      end if

      if (prll_on_off) then
         call eforce_reduce()
      end if

      if (qssbp) then
         call ssbp1()
         e_total = e_total + enssbp
      endif


      if ( torscstr ) then
         call ctors()
         e_total = e_total + e_ctors
      end if

      if ( alabar ) then 
         call alabarrier()
         e_total = e_total + e_alabar
      end if

      if ( rgcst ) then
         call rgyrconst()
         e_total = e_total + e_rgyr
      end if



c     redistribute forces acting on virtual prts among its relatives
      if (vp_flag) call vp_fdistrib()
c     end of vp


c$$$      if (dpot_ave_yes) then
c$$$         do 2 i=1,npt
c$$$            dpot_ave(1,i) = dpot_ave(1,i) + dpot(1,i)
c$$$            dpot_ave(2,i) = dpot_ave(2,i) + dpot(2,i)
c$$$            dpot_ave(3,i) = dpot_ave(3,i) + dpot(3,i)
c$$$ 2       continue
c$$$      end if


C     @
c       do 1 i=1,npt
c       do 1 j=1,3
c       if (dabs(dpot(j,i)).gt.50.d0) then
c       write(*,*)'i dpot ',i,dpot(j,i)
c       end if
c1      continue


      return
      end
