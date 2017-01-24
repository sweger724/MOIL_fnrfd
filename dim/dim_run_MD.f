      subroutine dim_run_MD(reach,reached,ufpt,mm)

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

        character*8 name
        integer namel 

        logical failure

        integer my_pe,i,k,l,iat1,iat2
        integer*8 istep,j
        integer kk,mat_idx(10*maxshak)

        double precision tmp,gradf
        double precision mat_val(10*maxshak)

        integer iii, reach,ufpt
        double precision reached(*)
        double precision rr_wfly,rx_wfly,ry_wfly,rz_wfly
        double precision dr(3,maxpt), dv(3,maxpt), fff
        !double precision xxx, CGtorsion
        logical w_stop(maxmono)
        double precision vel_save(3,maxpt)
        double precision mm(*)

        name = 'Brownian'
        namel= 8
        reach = -2

        if (matshak) then
           call build_matrix(mat_idx,ishak1,ishak2,nshak)
        endif

       if (wfly) then
        do 8 i=1,nwaters
            w_stop(i) = .false.
8       continue
       end if

       !call velinit(tempi,ntemp,tpo)

       do i=1, npt
         do l=1,3
           coor2(l,i)=coor(l,i)
         end do
       end do

C  initialize langevin dynamics friction coefficient
        if (langevin) then
          write(6,*)"Friction coef:",LDgamma
          
          c0 = dexp(-LDgamma*dt)
          c1 = (1.d0 - c0)/(LDgamma*dt)
          c2 = (1.d0 - c1)/(LDgamma*dt)
        else
          c0 = 1.0d0
          c1 = 1.0d0
          c2 = 0.5d0
        end if


          fff = 2.d0 + (- 3.d0 +
     &          (4.d0 - dexp(-LDgamma*dt))*dexp(-LDgamma*dt))
     &          / (LDgamma*dt)

          sigmaR = dsqrt (dt*kboltzmann*tempi(1)*fff /LDgamma)

        sigmaV= dsqrt(kboltzmann*tempi(1)*(1.d0-dexp(-2.d0*LDgamma*dt)))


          crv = kboltzmann*tempi(1) * (1.d0 - dexp(-LDgamma*dt))**2
     &        / (LDgamma * sigmaR * sigmaV)

c start the dynamics
c DO GLOBAL LOOP 
c nlist is the number of steps between the update of the
c non-bonded lists
c
c        ! align all centers to coor
c        do i = 1, Ncells+1
c          l=(i-1)*npt+1
c          call rmsd_weight(npt,coor(1,1),centers(1,l),rms,.false.,mm)
c        end do

        !xxx = -CGtorsion(torsions(1,1),torsions(1,2),
      !&        torsions(1,3),torsions(1,4),coor)

        do 16 istep=1,nstep,nlist

c
c for annealing box size squeeze will be passed the current size
c
         if (symanneal.and.esymyes .and. (istep.le.neqstep)) then
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
         !do i=1,12
         !  write(6,*)"XXX",coor(1,i),velo(1,i)
         !end do
         call nbondm()
         if (esymyes) call syminit()
         if (specl) call nbondm_spcl()

c call eforce to get energy and forces
c forces are initializing after each update of nonbonded
c and symmetry lists
c

         call eforce()
         if (langevin) then
           call getRandomDisplacement(dr,dv)
         else
           do i = 1,npt
             do l=1,3
               dr(l,i) = 0.d0
               dv(l,i) = 0.d0
             end do
           end do 
         end if

c  This tests force derivative
c         call deriv_test(nstep)
c         stop
         
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
         
         if (j.eq.1) then
           !do i=1,12
           !  write(6,*)"YYY",coor(1,i),velo(1,i)
           !end do 
           call info(stdo,j,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
         end if

c If Landau-Zener switching algorithm is on, check for a switch
C ### (lan_zen is not working, yet, in parallel)
c
          if (switch) call lan_zen()

c     prepare velo & vel_save arrays
          do  k=1,inofrz
            l = nofreez(k)
            ! here velo is calculated only as delta coor
            ! it is used by SHAKE
            if (langevin) then
              vel_save(1,l) =  c0*velo(1,l) - (c1 -c2)*dt*
     &          dpot(1,l)*invms(l) + dv(1,l)
              vel_save(2,l) =  c0*velo(2,l) - (c1 -c2)*dt*
     &          dpot(2,l)*invms(l) + dv(2,l)
              vel_save(3,l) =  c0*velo(3,l) - (c1 -c2)*dt*
     &          dpot(3,l)*invms(l) + dv(3,l)
            end if

            velo(1,l) = dt*( c1*velo(1,l) -
     &                 c2*dt*dpot(1,l)*invms(l)) + dr(1,l)
            velo(2,l) = dt*( c1*velo(2,l) -
     &                 c2*dt*dpot(2,l)*invms(l)) + dr(2,l)
            velo(3,l) = dt*( c1*velo(3,l) -
     &                 c2*dt*dpot(3,l)*invms(l)) + dr(3,l)

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
                rr_wfly = rx_wfly**2 + ry_wfly**2 +  rz_wfly**2
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
          tmp = 1.d0/dt
          do 3 k=1,inofrz
           l = nofreez(k)
           coor(1,l) = coor(1,l) + velo(1,l)
           coor(2,l) = coor(2,l) + velo(2,l)
           coor(3,l) = coor(3,l) + velo(3,l)

           if (langevin) then
             velo(1,l) =  vel_save(1,l) 
             velo(2,l) =  vel_save(2,l) 
             velo(3,l) =  vel_save(3,l) 
           else
             velo(1,l) = velo(1,l) * tmp
             velo(2,l) = velo(2,l) * tmp
             velo(3,l) = velo(3,l) * tmp
           end if
3         continue 
         
          if (mod(j,1).eq. 0) call TestInterfaces(j,reach,reached,
     1       ufpt,mm)
          if (reach .ne. -2) then
            ! write new interface to output file
C change by alfredo, want to write in the output.pth 
C the final points of every trajectory
C so I am moving the writing of the file
c            if (reach .eq. 0) then
c              call Write2File(uwcrd,npt,coor(1,1),0.d0,0,0)
c            end if
            return
          end if

c calculate new forces
          call eforce()
          if (langevin) then
            call getRandomDisplacement(dr,dv)
          else
            do i = 1,npt
              do l=1,3
                dr(l,i) = 0.d0
                dv(l,i) = 0.d0
              end do
            end do
          end if
          if (e_total.gt.1.d6) then
            write(stdo,*) "Simulation is unstable!"
            return
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
            velo(1,l) = velo(1,l) - c2*dt*invms(l)*dpot(1,l)
            velo(2,l) = velo(2,l) - c2*dt*invms(l)*dpot(2,l)
            velo(3,l) = velo(3,l) - c2*dt*invms(l)*dpot(3,l)
4         continue

c
c if matrix shake on, correct velocities
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
           do 6 k=1,ntemp
            tempi(k) = tempi(k) + dtemp(k)
6          continue
           ! update LD parameters
           fff = 2.d0 + (- 3.d0 +
     &          (4.d0 - dexp(-LDgamma*dt))*dexp(-LDgamma*dt))
     &          / (LDgamma*dt)

        sigmaR = dsqrt (dt*kboltzmann*tempi(1)*fff /LDgamma)

        sigmaV= dsqrt(kboltzmann*tempi(1)*(1.d0-dexp(-2.d0*LDgamma*dt)))


        crv = kboltzmann*tempi(1) * (1.d0 - dexp(-LDgamma*dt))**2
     &      / (LDgamma * sigmaR * sigmaV)
           ! end update

          end if
C*** END OF EQUILIBRIUM CALCULATION

c time to print out?
c ninfo   - number of steps between information dump
c ncoor   - number of steps between writing coordinates
c nvelo   - number of steps between writing velocities
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

          if ( (.not. no_scaling) .and. (.not. langevin)) then
c
c calculate current temperature, and scale to obtain desired one
c
           call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
           do 12 kk=1,ntemp
            if (dabs(curtemp(kk)-tempi(kk)).gt.nscalv) then
            write(stdo,105)(curtemp(k),k=1,ntemp)
105        format(1x,' Scaling velocities at current temperatures '
     1              ,(1x,f10.4))

           do k=1,inofrz
            l = nofreez(k)
            tmp = dsqrt(tempi(tpo(l))/curtemp(tpo(l)))
             velo(1,l) = velo(1,l)*tmp
             velo(2,l) = velo(2,l)*tmp
             velo(3,l) = velo(3,l)*tmp
           end do
                go to 13
            end if
12         continue
13         continue
         end if

         if (andersen) then
           call Andersen_thermostat(dt*coupling)
         end if

C check if time to remove rigid body motion
         nrigi = 1
         if ((mod(j,nrigi).eq.0) .and. (.not.freeze) .and.
     1          (.not.(esymyes.or.nori.or.eteth_yes))) then
            call ovrlpck(coor2,coor,dpot,velo,jpick,iorie,rms)
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

c ---------------------------------------------
          if (mod(j,ninfo).eq.0 .and. ninfo.ne.0) then
           call info(stdo,j,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
           write(6,*)"Temp:", tempi(1)
C           write(6,*)"Distance:", myRms
          end if
c ----------------------------------------

16      continue           
        rms = 0.d0
        reach = -1
        write(6,*)"NO INTERFACE REACHED!!!"
        return
        end


        subroutine TestInterfaces(step,reach,reached,ufpt,mm)

          implicit none

          include "COMMON/LENGTH.BLOCK"
          include "COMMON/COORD.BLOCK"
          include "COMMON/CONNECT.BLOCK"
          include "COMMON/CONVERT.BLOCK"
          include "COMMON/PATH2.BLOCK"
          include "COMMON/SEARCH.BLOCK"
          include "COMMON/VELOC.BLOCK"
          include "COMMON/DYNA.BLOCK"

          integer i,j,k,l,step, reach, bb, bb2, bb3,ufpt,k2,bb1
          double precision reached(*), CGtorsion, DistAngle
          double precision D1,D2,D0, prod,vvv,mindx, coor3(3,maxpt)
          double precision mm(*)
          character*9 text

          !write(6,*) "WWW:", myCell,myCell2

          if (Nreduced.ne.0) then
          
            call GetReducedCoors()
            !write(6,*)"CCC",myReduced(1)
            D1 = DistAngle(ReducedCoor(1,myCell),myReduced(1),Nreduced)
            D2 = DistAngle(ReducedCoor(1,myCell2),myReduced(1),Nreduced)
          else

            do i = 1, npt
              do l = 1,3
                coor3(l,i) = coor(l,i)
              end do
            end do

            do i= 1,Ncells+1
              l = (i-1)*npt+1
            call rmsd_weight(npt,coor3(1,1),centers(1,l),D0,.false.,
     &      mm)
            end do

            k = (myCell-1)*npt+1
            k2 = (myCell2-1)*npt+1
            call Distance(coor3(1,1),centers(1,k),D1)
            call Distance(coor3(1,1),centers(1,k2),D2)

          end if
             
            mindx = 1000
            do i = 1, Ncells
              if (Nreduced.ne.0) then
              D0=DistAngle(ReducedCoor(1,cell(i)),myReduced(1),Nreduced)
              else
                k = (cell(i)-1)*npt +1
                call Distance(coor3(1,1),centers(1,k),D0)
              end if
C              if (cell(i) .eq. myCell2) write(6,*)"RRR:", rms - rms1
              if (D1 .le. (D2-dx)) then
              if ( D0 - D1 .lt. mindx) then
                 mindx = D0 - D1
                 bb = i
              end if
              else 
                           if ( D0 - D2 .lt. mindx) then
                 mindx = D0 - D2
                 bb = i
              end if
              endif
            
c            write(6,*) i,mindx,dx,cell(i)
            end do

            !write(6,*) "MMM:", mindx
            !D1 = CGtorsion(2,4,6,8,coor)*pi180
            !D2 = CGtorsion(4,6,8,10,coor)*pi180
            !write(6,*)"BBB:",D1
            !write(6,*)"AAA:",D2
            if (mindx .lt. -dx) then

              reach = bb
              if (D1 .le. (D2-dx)) then
              bb1 = myCell
              bb2 = myCell
              bb3 = myCell2
              else
              bb1 = myCell2
              bb2 = myCell
              bb3 = myCell2
              endif
             
              text = "reached "
              if ( reached(bb) .eq. 0 ) then
                text="Nreached "
                reach = 0
              end if
              reached(bb) = reached(bb) + 1
C this is for a dipeptide:
c              D1 = CGtorsion(2,4,6,8,coor)*pi180
c              D2 = CGtorsion(4,6,8,10,coor)*pi180
c              write(6,*)"BBB:",cell(bb),D1
c              write(6,*)"AAA:",cell(bb),D2
c              write(6,'(a,i4,1x,i4,1x,i4,a,i7,a,f9.4)') text,bb3,
c     &                 bb2,cell(bb)," after ",step,
c     &                     " steps ",weight
              write(ufpt,1005) text,bb3,
     &                 bb2,bb1,cell(bb)," after ",step,
     &                     " steps ",weight
1005          format(a,i4,1x,i4,1x,i4,1x,i4,a,i7,a,f9.4)
              call wener(6)
C  do the writting here -alfredo
              call Write2File(uwcrd,npt,coor(1,1),0.d0,0,0)
              return

            end if

        end
