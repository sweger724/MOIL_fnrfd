        program dynamics
        implicit none
c
c program for burning CPU(s)
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/PT.BLOCK'
        include 'COMMON/REPWALL.BLOCK'
        include 'COMMON/EFIELD.BLOCK'
        include 'COMMON/PRESS.BLOCK'

        character*8 name
        integer namel 

        logical failure

        integer k,l,iat1,iat2
	integer k1,k2
        integer ndegf_save
        integer j, istep
        integer kk,mat_idx(10*maxshak)

        double precision tmp,gradf
        double precision mat_val(10*maxshak)

        integer iii,i
        double precision rr_wfly,rx_wfly,ry_wfly,rz_wfly
        logical w_stop(maxmono), PT_accept

        integer wbin

c initialize variables
c
        name = 'dynamics'
        silent = .false.
        namel= 8

        stdi   = 5
        stdo   = 6
        stderr = 0
        wbin = 1        

        call init_paral()
        call init_io()
        call init_dyna()
        call init_press()

c printout energy parameters for the record
c
        call dump_param(stdo)

c print dynamics parameters for the record
c
        call dump_dyna(stdo)

        if (matshak) then
           call build_matrix(mat_idx,ishak1,ishak2,nshak)
        endif

       if (wfly) then
        do i=1,nwaters
            w_stop(i) = .false.
        end do
       end if

c if parallel call twice for nbondm (the same sequence below) 
c  due to load balancing

       if (prll_on_off) then
        if (esymyes) call squeeze(a,b,c)
        call nbondm()
        if (esymyes) call syminit()
        if (specl) call nbondm_spcl()
       end if


c start the dynamics
c DO GLOBAL LOOP 
c nlist is the number of steps between the update of the
c non-bonded lists
c
       do 6 istep=start_dyna,nstep,nlist

        if ( mod(istep,nsynch) .eq. (nsynch-1) ) then
          call broadcast_state()
        end if
c
c for annealing box size squeeze will be passed the current size
c
         if (symanneal.and.esymyes .and. (istep.le.neqstep)) then
            a=a+dxtra
            b=b+dytra
            c=c+dztra
            if (my_pe.eq.0) then
             if (my_pe .eq. 0) write (stdo,110)a,b,c
110          format(1x,'Current symmetry box size (x,y,z):',3(1x,f10.3))
            end if
         endif
c
c call list to get non-bonded list
c
c squeeze must be called before nbond or the list is wrong
c for squeezed particles.
c

         if (esymyes) call squeeze()
         call nbondm()
         if (esymyes) call syminit()
         if (specl) call nbondm_spcl()

c call eforce to get energy and forces
c forces are initializing after each update of nonbonded
c and symmetry lists
c

c  This tests force derivative
c         call deriv_test(nstep)
c         stop

         call eforce()
         if ( istep.eq. start_dyna) call wener(stdo)      
        
        if (sdyes) then
c If force is too large try to fix things by short minimization
c
                call force_norm(gradf)
                if (gradf.gt.fmax) then
                        call stdc(failure,velo)
                        if (failure) then
                         call alert(name,namel,'STDC FAILED!!!',14,1)
                        end if
                        call velinit(tempi,ntemp,tpo)
                end if
        end if
        if ((.not.esymyes) .and. (.not.nori) .and. (.not.qssbp)
     1          .and. (.not.eteth_yes) .and. (.not.freeze)) then
                call ovrlpck(coor2,coor,dpot,velo,jpick,iorie,rms)
                call prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)
        end if

c Prepare a reference coordinate for shake
c
          if (shakl.or.shakb .and. nshak.gt.0) then
           do 1 k=1,nshak
            iat1 = ishak1(k)
            iat2 = ishak2(k)
            cooref(1,k) = coor(1,iat1) - coor(1,iat2)
            cooref(2,k) = coor(2,iat1) - coor(2,iat2)
            cooref(3,k) = coor(3,iat1) - coor(3,iat2)
1          continue
           if (matshak) then
              call build_matrix_from_idx(mat_val,mat_idx)
           end if
        endif

C HERE START THE DYNAMICS LOOP
C
         !call info(stdo,0,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
         do 16 j=istep,min(nstep,istep+nlist-1)


c If Landau-Zener switching algorithm is on, check for a switch
C ### (lan_zen is not working, yet, in parallel)
c
          if (switch) call lan_zen()
        
          ! calculate kinetic energy for print out
          if (ninfo.ne.0 .and. mod(j,ninfo).eq.0) then
              call calculate_kinetic_energy()
          endif 


c Make a step for unfrozen particles (velocity is used to store the step)
c dividing finally by dt will recover the magnitude required for velocity
c calculations
          do  k=1,inofrz
            l = nofreez(k)
            velo(1,l) =  velo(1,l)*dt - factor1(k)*dpot(1,l)
            velo(2,l) =  velo(2,l)*dt - factor1(k)*dpot(2,l)
            velo(3,l) =  velo(3,l)*dt - factor1(k)*dpot(3,l)
          end do


c if (shake) call shake to fix position of shaked particles
          if (shakl.or.shakb .and. nshak.gt.0) then
             if (matshak) then
                call conjugate_grad_shakept (mat_idx,mat_val)
             else
                call shakept(itershak)
             endif
        end if

c call matrix shake if (shakm)
          if (shakm) then
                call mshakpt(dt,dt2)
          end if

c check if water fly away
        if (wfly) then
            do 17 i=1,nwaters
                iii = dpoipt(idxtip3(i))-2
                rx_wfly = coor(1,iii)
                ry_wfly = coor(2,iii)
                rz_wfly = coor(3,iii)
                rr_wfly = rx_wfly**2 + ry_wfly**2 + rz_wfly**2
                if (rr_wfly.gt.3600) then
                  velo(1,iii) = 0.d0
                  velo(2,iii) = 0.d0
                  velo(3,iii) = 0.d0
                  velo(1,iii+1) = 0.d0
                  velo(2,iii+1) = 0.d0
                  velo(3,iii+1) = 0.d0
                  velo(1,iii+2) = 0.d0
                  velo(2,iii+2) = 0.d0
                  velo(3,iii+2) = 0.d0
                  if (.NOT.w_stop(i)) then
               if (my_pe .eq. 0) write(6,*)' H2O number ', idxtip3(i),
     1                   ' was stopped! ',rr_wfly
                   w_stop(i) = .true.
                end if
             end if
17            continue
        end if

        if (mod(j,ninfo).eq.0 .and. ninfo.ne.0) then
            call write_press()    
        endif

c Make a coordinate step and prepare the velocity calculation
          tmp = 1.d0/dt
          do 3 k=1,inofrz
           l = nofreez(k)
           coor(1,l) = coor(1,l) + velo(1,l)
           coor(2,l) = coor(2,l) + velo(2,l)
           coor(3,l) = coor(3,l) + velo(3,l)
           velo(1,l) = velo(1,l)*tmp
           velo(2,l) = velo(2,l)*tmp
           velo(3,l) = velo(3,l)*tmp
3         continue

c calculate new forces
          call eforce()
C          if (e_total.gt.1.d6) then
C            if (my_pe .eq. 0) write(stdo,*) "Simulation is unstable!"
C            stop
C          end if
          
          if (sdyes) then
                call force_norm(gradf)
                if (gradf.gt.fmax) then
                   call stdc(failure,velo)
                   if (failure) then
                        call alert(name,namel,'STDC FAILED!!!',14,1)
                   end if
                   call velinit(tempi,ntemp,tpo)
                end if
           end if

           if (shakl.or.shakb .and. nshak.gt.0) then
              do 5 k=1,nshak
                 iat1 = ishak1(k)
                 iat2 = ishak2(k)
                 cooref(1,k) = coor(1,iat1) - coor(1,iat2)
                 cooref(2,k) = coor(2,iat1) - coor(2,iat2)
                 cooref(3,k) = coor(3,iat1) - coor(3,iat2)
 5            continue
              if (matshak) then
                 call build_matrix_from_idx(mat_val,mat_idx)
              endif
           end if

c calculate new velocities 
          do 4 k=1,inofrz
           l = nofreez(k)
           velo(1,l) = velo(1,l) - factor2(k)*dpot(1,l)
           velo(2,l) = velo(2,l) - factor2(k)*dpot(2,l)
           velo(3,l) = velo(3,l) - factor2(k)*dpot(3,l)
4         continue

        if (shakm) then
           call mshakvl(dt)
        endif

c if shake is on call shakevl
          if (shakl.or.shakb .and. nshak.gt.0) then
             if (matshak) then
                call conjugate_grad_shakevl (mat_idx,mat_val,epshakv)
             else
                call shakevl(epshakv,itershak)
             endif
         end if
c Are we still in equilibration period
          if (j.lt.neqstep) then

C *** EQUILIBRIUM
c the current temperature should be the average between tempi & tempf
c linear heating is assumed
c
           do k=1,ntemp
             tempi(k) = tempi(k) + dtemp(k)
           end do
c
c calculate current temperature(s), and scale velocities
c to obtain desired one(s)
c
           call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
           if (mod(j,ninfo).eq.0 .and. my_pe.eq.0) then
            write(stdo,100)(curtemp(k),k=1,ntemp)
            write(stdo,101)(tempi(k),k=1,ntemp)
100         format(1x,' Before scaling... temp ',5(1x,f14.4))
101         format(1x,'             requ. temp ',5(1x,f14.4),/)
           end if
           do 14 k=1,inofrz
            l = nofreez(k)
            tmp = dsqrt(tempi(tpo(l))/curtemp(tpo(l)))
            !write(6,*)"VVV:",velo(1,l)
            velo(1,l) = velo(1,l)*tmp
            velo(2,l) = velo(2,l)*tmp
            velo(3,l) = velo(3,l)*tmp
14         continue
          end if
C*** END OF EQUILIBRIUM CALCULATION

c time to print out?
c ninfo   - number of steps between information dump
c ncoor   - number of steps between writing coordinates
c nvelo   - number of steps between writing velocities
c nscalv  - temperature window width that specifies when rescaling happens
c ---------------------------------------------
       if (mod(j,ninfo).eq.0 .and. ninfo.ne.0) then
           call info(stdo,j,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
           if (my_pe.eq.0 .and. urst_crd .ne. -1) call up_res(j)
       end if
c ----------------------------------------
c
c time for new velocity assignment?
c
          if (newv.ne.0) then
             if (mod(j,newv).eq.0) then
              if (my_pe.eq.0)
     1          write(stdo,*)' New velocities are assigned at step ',j
                call velinit(tempi,ntemp,tpo)
           end if
          end if

C time to write coordinates?
C
          if (ncoor.ne.0 .and. my_pe.eq.0) then
           if (mod(j,ncoor).eq.0) then
           call wdyncrd(uwcrd,nstep/ncoor,j/ncoor,inofrz,nofreez,wbin)
           end if
          end if

C time to write velocities
C
          if (nvelo.ne.0 .and. my_pe.eq.0) then
           if (mod(j,nvelo).eq.0)
     1          call wdynvel(uwvel,nstep/nvelo,j/nvelo,inofrz,nofreez)
          end if


          if (Ltpick) then
c
c compute temperature for a subset of atoms
	 call multemp(velo,ptms,tpoint,tsele,tpo,tgroup,
     1         curtemp,ntemp)
	  write(*,*)' temperature selected for ',tsele,' prtc '
          write(*,*)'  is ',(curtemp(k)*tgroup(k)/(3*tsele),k=1,ntemp)
         end if

C check if time to rescale velocities
         if (.not.no_scaling) then
c
c calculate current temperature, and scale to obtain desired one
c
          if (.not.scale_away) then
           call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
           do 12 kk=1,ntemp
            if (dabs(curtemp(kk)-tempi(kk)).gt.nscalv) then
             if (my_pe.eq.0) then
              write(stdo,105)(curtemp(k),k=1,ntemp)
             end if
105          format(1x,' Scaling velocities at current temperatures '
     1              ,(1x,f10.4))
             do k=1,inofrz
              l = nofreez(k)
              tmp = dsqrt(tempi(tpo(l))/curtemp(tpo(l)))
              velo(1,l) = velo(1,l)*tmp
              velo(2,l) = velo(2,l)*tmp
              velo(3,l) = velo(3,l)*tmp
             end do
            end if
12	   continue
          else if (scale_away) then
c this scaling is specific to water molecules that are far away (no_of_away).
c Each water molcule has 3 atoms, so the ndegf should be 3*no_of_away.
c However each water molecule has three constraints. So the ttal ndegf for
c a water molecule is 6 or 2*no_of_atom_per_molecule, and hence the
c 2*no_of_away
c
             ndegf_save = tgroup(1)
             tgroup(1) = 2 * no_of_away
             call multemp(velo,ptms,far_water,no_of_away,tpo,tgroup,
     1          curtemp,ntemp)
             tmp = dsqrt(tempi(tpo(1))/curtemp(tpo(1)))
             tgroup(1) = ndegf_save
c
             write(*,*)" scale away : ",tempi(tpo(1)),curtemp(tpo(1))
             do k=1,no_of_away
               l = far_water(k)
               velo(1,l) = velo(1,l)*tmp
               velo(2,l) = velo(2,l)*tmp
               velo(3,l) = velo(3,l)*tmp
             end do
	  end if
         end if

C check if time to remove rigid body motion
      if ((mod(j,nrigi).eq.0) .and. (.not.freeze) .and.
     1          (.not.(esymyes.or.nori.or.eteth_yes))) then
C### to be parallelized later
            call ovrlpck(coor2,coor,dpot,velo,jpick,iorie,rms)
c the multemp below is debugging
C
            call prbm(velo,coor,ptms,grdlx,grdly,grdlz,npt)

            if (shakl.or.shakb .and. nshak.gt.0) then
               do 155 k=1,nshak
                  iat1 = ishak1(k)
                  iat2 = ishak2(k)
                  cooref(1,k) = coor(1,iat1) - coor(1,iat2)
                  cooref(2,k) = coor(2,iat1) - coor(2,iat2)
                  cooref(3,k) = coor(3,iat1) - coor(3,iat2)             
 155           continue
               if (matshak) then
                  call build_matrix_from_idx(mat_val,mat_idx)
               endif
            end if
         end if

        ! if we run parallel tempering ...
         if (PT_run .and. j.gt.neqstep.and.mod(j,nswap).eq.0) then
            call swap_conformations(j,PT_accept)
            if (PT_accept) goto 6            
         end if

16      continue         
6       continue
        
        call parallel_end()
        stop
        end
