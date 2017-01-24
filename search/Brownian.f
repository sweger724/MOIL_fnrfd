      subroutine Brownian(printEach, ireq,finished)

      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SEARCH.BLOCK'

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
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/TETHER.BLOCK'
        include 'COMMON/RESTART.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/EBALL.BLOCK'
        include 'COMMON/RGYRCST.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'
        include 'COMMON/LD.BLOCK'

        integer printEach

        character*8 name
        integer namel 

        logical failure

        integer my_pe,i,j,k,l,iat1,iat2,istep
        integer kk,mat_idx(10*maxshak)

        double precision tmp,gradf
        double precision mat_val(10*maxshak)

        integer iii
        double precision rr_wfly,rx_wfly,ry_wfly,rz_wfly
        double precision dr(3,maxpt), dv(3,maxpt), fff
        double precision coor3(3,maxpt),initTemp
        logical w_stop(maxmono)

        logical finished
        integer ireq
      
        integer wbin
 
        finished = .FALSE.
       
        name = 'Brownian'
        namel= 8
        nscalv=1

        wbin = 1

        if (matshak) then
           call build_matrix(mat_idx,ishak1,ishak2,nshak)
        endif

       if (wfly) then
        do 8 i=1,nwaters
            w_stop(i) = .false.
8       continue
       end if
       initTemp = tempi(1)
       call velinit(tempi,ntemp,tpo)

C  initialize langevin dynamics friction coefficient
        LDgamma = 0.3d0

        c0 = dexp(-LDgamma*dt)
        c1 = (1.d0 - c0)/(LDgamma*dt)
        c2 = (1.d0 - c1)/(LDgamma*dt)

        fff = 2.d0 + (- 3.d0 +
     &          (4.d0 - dexp(-LDgamma*dt))*dexp(-LDgamma*dt))
     &          / (LDgamma*dt)

        sigmaR = dsqrt (dt*kboltzmann*tempi(1)*fff /LDgamma)

        sigmaV= dsqrt(kboltzmann*tempi(1)*(1.d0 - dexp(-2.d0*LDgamma*dt)))


        crv = kboltzmann*tempi(1) * (1.d0 - dexp(-LDgamma*dt))**2
     &      / (LDgamma * sigmaR * sigmaV)

        do i=1, npt
          do l=1,3
            coor2(l,i)=coor(l,i)
          end do
        end do


c start the dynamics
c DO GLOBAL LOOP 
c nlist is the number of steps between the update of the
c non-bonded lists
c

        do 16 istep=1,nstep,nlist

c
c for annealing box size squeeze will be passed the current size
c
         if (symanneal.and.esymyes .and. (istep.gt.neqstep)) then
            a=a+dxtra
            b=b+dytra
            c=c+dztra
            if (my_pe.eq.0) then
             write (stdo,110)a,b,c
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

         call eforce()
C         call wener(stdo)

         call getRandomDisplacement(dr,dv)
C         if (istep.eq.1) call wener(stdo)       

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
         do 16 j=istep,min(nstep,istep+nlist-1)


c If Landau-Zener switching algorithm is on, check for a switch
C ### (lan_zen is not working, yet, in parallel)
c
          if (switch) call lan_zen()

c Make a step for unfrozen particles (velocity is used to store the step)
c dividing finally by dt will recover the magnitude required for velocity
c calculations

c@
c       write(*,*) ' In dyna after loop 2  '
c       do 33 k=1,inofrz
c       l = nofreez(k)
c       write(*,*)' velo',l,'=',(velo(kk,l),kk=1,3)
c33     continue
c@
c         if (prll_on_off) call gather_vel()

c if (shake) call shake to fix position of shaked particles
          if (shakl.or.shakb .and. nshak.gt.0) then
             if (matshak) then
                call conjugate_grad_shakept (mat_idx,mat_val)
             else
                call shakept(itershak)
             endif
        end if

c call matrix shake if (shakm)
c         write(*,*) 'shakm = ',shakm
          if (shakm) then
                call mshakpt(dt,dt2)
c@@             if (prll_on_off) call gather_wat()
          end if

c check if water fly away
c@
c       write(*,*) ' wfly = ',wfly
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
                 write(*,*)' H2O number ', idxtip3(i),
     1              ' was stopped! ',rr_wfly
                 w_stop(i) = .true.
                end if
             end if
17            continue
        end if

c Make a coordinate step and prepare the velocity calculation
C         write(*,*)'VVV: ',velo(1,1),velo(2,1),velo(3,1)
          do 3 k=1,inofrz
           l = nofreez(k)
           coor(1,l) = coor(1,l) + c1*dt*velo(1,l) - 
     &                 c2*dt**2*dpot(1,l)*invms(l) + dr(1,l)
           coor(2,l) = coor(2,l) + c1*dt*velo(2,l) -
     &                 c2*dt**2*dpot(2,l)*invms(l) + dr(2,l)   
           coor(3,l) = coor(3,l) + c1*dt*velo(3,l) -
     &                 c2*dt**2*dpot(3,l)*invms(l) + dr(3,l)

           velo(1,l) =  c0*velo(1,l) - (c1 -c2)*dt*
     &         dpot(1,l)*invms(l) + dv(1,l)
           velo(2,l) =  c0*velo(2,l) - (c1 -c2)*dt*
     &         dpot(2,l)*invms(l) + dv(2,l)
           velo(3,l) =  c0*velo(3,l) - (c1 -c2)*dt*
     &         dpot(3,l)*invms(l) + dv(3,l) 
3         continue

          call TestDone(ireq,finished)
          if (finished) then
C              write(6,*)"All structures collected."
              ! report to master that I realized that round is finished
              coor(1,1) = 0.d0
              coor(1,2) = 0.d0
              call Send_Double(coor(1,1),3*npt,0,100)
              tempi(1) =  initTemp 
              return
          end if

          if (mod(j,50).eq.0) then
            do i = 1, npt
              do l =1,3
                coor3(l,i) = coor(l,i)
              end do
            end do
C            call rmsd_weight(npt,coor2(1,1),coor3(1,1),rms,.false.,ptms)
            call Distance(coor3(1,1),coor2(1,1),rms)
C            write(stdo,*)"Distance: ",rms
            if (rms.gt.cutoff) then
C              write(6,*)"Simulation finished at step:",j
C              write(stdo,*)"Sending result."
              coor3(1,npt+1) = j
              call Send_Double(coor3(1,1),3*npt+1,0,100)
              tempi(1) = initTemp 
              return
            end if
          end if

c calculate new forces
          call eforce()
          call getRandomDisplacement(dr,dv)
          if (e_total.gt.1.d6) then
            write(stdo,*) "Simulation is unstable!"
            stop
          end if
          
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
           velo(1,l) = velo(1,l) - 1.d0*c2*dt*invms(l)*dpot(1,l)
           velo(2,l) = velo(2,l) - 1.d0*c2*dt*invms(l)*dpot(2,l)
           velo(3,l) = velo(3,l) - 1.d0*c2*dt*invms(l)*dpot(3,l)
4         continue


c
c if matrix shake on, correct velocities
        if (shakm) then
           call mshakvl(dt)
c@@        if (prll_on_off) call gather_wat()
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
          if (j.gt.neqstep) then

C *** EQUILIBRIUM
c the current temperature should be the average between tempi & tempf
c linear heating is assumed
c
           do k=1,ntemp
            tempi(k) = tempi(k) + dtemp(k)
           end do
          ! update LD parameters
           fff = 2.d0 + (- 3.d0 +
     &          (4.d0 - dexp(-LDgamma*dt))*dexp(-LDgamma*dt))
     &          / (LDgamma*dt)

        sigmaR = dsqrt (dt*kboltzmann*tempi(1)*fff /LDgamma)

        sigmaV= dsqrt(kboltzmann*tempi(1)*(1.d0 - dexp(-2.d0*LDgamma*dt)))


        crv = kboltzmann*tempi(1) * (1.d0 - dexp(-LDgamma*dt))**2
     &      / (LDgamma * sigmaR * sigmaV)
           ! end update        
         end if

C*** END OF EQUILIBRIUM CALCULATION

c time to print out?
c ninfo   - number of steps between information dump
c ncoor   - number of steps between writing coordinates
c nvelo   - number of steps between writing velocities
c nscalv  - number of steps between velocity scaling
c ---------------------------------------------
          if (mod(j,ninfo).eq.0 .and. ninfo.ne.0) then
           call info(stdo,j,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
           write(6,*)"Temp:", tempi(1)
           write(stdo,*)"Distance: ",rms
C           call up_res(j)
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

C check if time to remove rigid body motion
         if ((mod(j,nrigi).eq.0) .and. (.not.freeze) .and.
     1          (.not.(esymyes.or.nori.or.eteth_yes))) then
C### to be parallelized later
            call ovrlpck(coor2,coor,dpot,velo,jpick,iorie,rms)
c the multemp below is debugging
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

16      continue           
        rms = 0.d0
        ! give signal that simulation did not escape
        coor(1,1) = 0.d0
        call Send_Double(coor(1,1),3*npt,0,100)
        tempi(1) = initTemp 
        return
        end

          subroutine TestDone(ireq,finished)

          implicit none

          include "mpif.h"

          integer ireq,  status(MPI_STATUS_SIZE),rc
          logical finished

           call MPI_TEST(ireq,finished,status,rc)
        end

        subroutine ReceiveFinishSignal(ireq)

          implicit none

          include "COMMON/LENGTH.BLOCK"
          include "COMMON/SEARCH.BLOCK"
          include "mpif.h"

          integer ireq,rc

           call MPI_IRECV(doneB,1,MPI_Integer,0,155,
     &                    MPI_COMM_WORLD,ireq,rc)

        end
