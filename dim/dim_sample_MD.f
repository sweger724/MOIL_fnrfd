      subroutine dim_sample_MD(reached,utmpx)

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
        include 'COMMON/PARALLEL.BLOCK'

        integer utmpx

        character*8 name
        integer namel 

        logical failure, forward

        integer i,k,l,iat1,iat2
        integer*8 j, istep
        integer kk,mat_idx(10*maxshak)

        double precision tmp,gradf,coorR(3,maxpt)
        double precision mat_val(10*maxshak)

        integer iii, reached,js
        double precision dr(3,maxpt), dv(3,maxpt), fff
        logical w_stop(maxmono), OK, bad, sample, reverse
        double precision coor_save(3,maxpt), velo_save(3,maxpt)
        double precision vel_save(3,maxpt),CGtorsion,velo_my(3,maxpt)
        real ranx, maxlength

        name = 'Brownian'
        namel= 8
        reverse = .false.
        reached = -2
c       steps performed in sampling        
        js = 0    

        if (matshak) then
           call build_matrix(mat_idx,ishak1,ishak2,nshak)
        endif

       call velinit(tempi,ntemp,tpo)

       do i=1, npt
         do l=1,3
           coorR(l,i)=coor(l,i)
         end do
       end do

       total = nsave

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
     &      / (LDgamma * sigmaR * sigmaV)

c start the dynamics
c DO GLOBAL LOOP 
c nlist is the number of steps between the update of the
c non-bonded lists
c
c         turn on the interface constraining potential
        e_Kumbyes = .true.
        sample = .true.

        do 16 istep=1,nstep,nlist

          if ( mod(istep,nsynch) .eq. 1 ) then
            call broadcast_state()
          end if

c  restart
177      continue
c
c for annealing box size squeeze will be passed the current size
c
         if (symanneal.and.esymyes .and. (istep.le.neqstep)) then
            a=a+dxtra
            b=b+dytra
            c=c+dztra
             write (stdo,110)a,b,c
110          format(1x,'Current symmetry box size (x,y,z):',3(1x,f10.3))
         end if
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
c        call deriv_test(1,32)
c        stop


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
         weight = dexp(e_Kumb/(kboltzmann*tempi(1)))
         if (istep.eq.1) call wener(stdo)      

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
                call ovrlpck(coorR,coor,dpot,velo,jpick,iorie,rms)
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

c     prepare velo & velo_save arrays
          do  k=1,inofrz
            l = nofreez(k)
c             here velo is calculated only as delta coor
c             it is used by SHAKE
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

c calculate new forces
          call eforce()
          weight = dexp(e_Kumb/(kboltzmann*tempi(1)))
         
          if (e_total.gt.1.d5 .and. .false.) then
            write(stdo,*) "Simulation is unstable!", sample
            do i = 1, npt
              do l = 1,3
                coor(l,i) = coorR(l,i)
              end do
            end do
            call velinit(tempi,ntemp,tpo)
            sample = .true.
            e_Kumbyes = .true.
            goto 177
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
          if (j.le.neqstep/2) then

C *** EQUILIBRIUM
c the current temperature should be the average between tempi & tempf
c linear heating is assumed
c
           do 6 k=1,ntemp
            tempi(k) = tempi(k) + dtemp(k)
6          continue
c            update LD parameters
           fff = 2.d0 + (- 3.d0 +
     &          (4.d0 - dexp(-LDgamma*dt))*dexp(-LDgamma*dt))
     &          / (LDgamma*dt)

        sigmaR = dsqrt (dt*kboltzmann*tempi(1)*fff /LDgamma)

        sigmaV= dsqrt(kboltzmann*tempi(1)*(1.d0-dexp(-2.d0*LDgamma*dt)))


        crv = kboltzmann*tempi(1) * (1.d0 - dexp(-LDgamma*dt))**2
     &      / (LDgamma * sigmaR * sigmaV)
c            end update

          end if
C*** END OF EQUILIBRIUM CALCULATION

c time for new velocity assignment?
c
          if (newv.ne.0) then
             if (mod(j,newv).eq.0) then
                write(stdo,*)' New velocities are assigned at step ',j
                call velinit(tempi,ntemp,tpo)
           end if
          end if

C check if time to rescale velocities
         if ( (.not. no_scaling) .and. sample
     &        .and. (.not.langevin)           ) then
c
c calculate current temperature, and scale to obtain desired one
c
          call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
         if (j.lt. neqstep/2 .or. 
     1       dabs(curtemp(1)-tempi(1)).gt.nscalv ) then
           do 12 kk=1,ntemp
            if (dabs(curtemp(kk)-tempi(kk)).gt.nscalv ) 
     1       write(stdo,105)(curtemp(k),k=1,ntemp)
105        format(1x,' Scaling velocities at current temperatures '
     1              ,(1x,f10.4))

           do k=1,inofrz
            l = nofreez(k)
            tmp = dsqrt(tempi(tpo(l))/curtemp(tpo(l)))
             velo(1,l) = velo(1,l)*tmp
             velo(2,l) = velo(2,l)*tmp
             velo(3,l) = velo(3,l)*tmp
           end do
12         continue

c       if false        
        end if 
          if (andersen) then
            call Andersen_thermostat(dt*coupling)
          end if

         end if

          if (sample) then
            js = js + 1
          end if

          if (j .le. neqstep .and. mod(js,ncoor).eq.0) then
             if (verbose) write(6,*)"VVV:",weight, j,js,e_Kumb
          end if

          if (j .gt. neqstep .and. sample) then
           if (mod(js,ncoor).eq.0) then
              if (verbose) write(6,*)"WWW:",weight, j,js,e_Kumb
              call TestInterfaces(j,OK)
              if (OK) then
                call RANLUX(ranx,1)
                maxlength = 3 / ranx
                sample = .false.

                c0 = 1.d0
                c1 = 1.d0
                c2 = 0.5d0

                e_Kumbyes = .false.

                do i = 1,npt
                  do l = 1,3
                    coor_save(l,i) = coor(l,i)
                    velo_save(l,i) = velo(l,i)
                    velo(l,i) = -velo(l,i)
                  end do
                end do
c@
c Note clear why velocity is inititated afterit was reversed!
c the line below is therefore commented out. RE 7/15/2010
c
c                call velinit(tempi,ntemp,tpo)
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

                
                do i = 1,npt
                  do l = 1,3
                    velo_my(l,i) = velo(l,i)
                  end do
                end do

        
                forward=.false.
             
               if (verbose) write(6,*)"Propper interface point."
              end if
           end if

          end if
          
          if ( .not. sample .and. forward) then
c@
c   testing sample and forward. Not sure bad was initialized
c RE
            if (bad) then
              if (verbose) write(6,*) "NOT reactive."
            else
c@
c not clear what reverse is all about. Not defined clearly. RE
              call TestInterfaces2(j,OK,bad,reverse)
c@
            end if
            
            if (OK .and. (.not. bad)) then
              if (my_pe .eq. 0) then
                call Write2File(uwcrd,npt,Wcoor(1,1),weight,0,0)
                call Write2File(uwvel,npt,Wvelo(1,1),weight,0,0)
              end if
              !if (.not. reverse)
     *        !call info(stdo,j,nofreez,inofrz,tpo,tgroup,ntemp,curtemp)
              nsave = nsave - 1
              www(total-nsave) = weight
              write(6,*) "Phase point written!", weight, reverse
              
              if (nsave .eq. 0 ) return
              
            end if
            if (bad .or. OK) then
              sample = .true.
              e_Kumbyes = .true.
              forward = .false.

              if (langevin) then
                c0 = dexp(-LDgamma*dt)
                c1 = (1.d0 - c0)/(LDgamma*dt)
                c2 = (1.d0 - c1)/(LDgamma*dt)
              else
                c0 = 1.0d0
                c1 = 1.0d0
                c2 = 0.5d0
              end if

              do i = 1,npt
                do l = 1,3
                   coor(l,i) = coor_save(l,i)
                   velo(l,i) = velo_save(l,i)
                end do
              end do
                if (esymyes) call squeeze()
                call nbondm()
                if (esymyes) call syminit()
                if (specl) call nbondm_spcl()
                call eforce()
             
              if (shakl.or.shakb .and. nshak.gt.0) then
                do k=1,nshak
                  iat1 = ishak1(k)
                  iat2 = ishak2(k)
                  cooref(1,k) = coor(1,iat1) - coor(1,iat2)
                  cooref(2,k) = coor(2,iat1) - coor(2,iat2)
                  cooref(3,k) = coor(3,iat1) - coor(3,iat2)
                end do
                if (matshak) then
                  call build_matrix_from_idx(mat_val,mat_idx)
                endif
              end if

 
            end if
          end if 

          if (.not. sample .and. .not. forward) then
            call TestInterfaces3(j,OK,bad,reverse)
            if (OK .or. bad) then
              forward = .true.
              do i = 1,npt
                do l = 1,3
                   coor(l,i) = coor_save(l,i)
                   velo(l,i) = -velo_my(l,i)
                end do
              end do

               if (esymyes) call squeeze()
               call nbondm()
               if (esymyes) call syminit()
               if (specl) call nbondm_spcl()
               call eforce()
               
              if (shakl.or.shakb .and. nshak.gt.0) then
                do k=1,nshak
                  iat1 = ishak1(k)
                  iat2 = ishak2(k)
                  cooref(1,k) = coor(1,iat1) - coor(1,iat2)
                  cooref(2,k) = coor(2,iat1) - coor(2,iat2)
                  cooref(3,k) = coor(3,iat1) - coor(3,iat2)
                end do
                if (matshak) then
                 call build_matrix_from_idx(mat_val,mat_idx)
                endif
              end if

            end if
          end if
         
        if (sample .and. langevin) then
           call getRandomDisplacement(dr,dv)
        else
          do i = 1,npt
            do l=1,3
              dr(l,i) = 0.d0
              dv(l,i) = 0.d0
            end do
          end do
        end if 

C check if time to remove rigid body motion
         nrigi = 1
         if ((mod(j,nrigi).eq.0) .and. (.not.freeze) .and.
     1          (.not.(esymyes.or.nori.or.eteth_yes))) then
            call ovrlpck(coorR,coor,dpot,velo,jpick,iorie,rms)
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
          end if
c ----------------------------------------         

16      continue           
        rms = 0.d0
        reached = -1
        return
        end


        subroutine TestInterfaces(step,OK)

          implicit none

          include "COMMON/LENGTH.BLOCK"
          include "COMMON/COORD.BLOCK"
          include "COMMON/CONNECT.BLOCK"
          include "COMMON/PATH2.BLOCK"
          include "COMMON/SEARCH.BLOCK"
          include "COMMON/VELOC.BLOCK"
          include "COMMON/DYNA.BLOCK"
          include "COMMON/PARALLEL.BLOCK"
          include "COMMON/OVERLAP.BLOCK"

          integer i,l,step, ss, bb, reached,k1,k2,k
          double precision prod,size,coor3(3,maxpt)
          double precision D0, D1, D2, mindx

          double precision lower, upper, neigh,DistAngle
          
          logical OK
          logical firstT
          data firstT/.true./

          save firstT,lower,upper,neigh

          if (firstT) then
            lower = dx - dx_umbrella
            upper = dx
            neigh = dx
            firstT = .false.
          end if

          OK = .true.

          if (Nreduced .gt. 0) then 
            if (my_pe .ne. num_pes-1) then
              call GetReducedCoors()
            end if
  
            D1=DistAngle(ReducedCoor(1,myCell), myReduced(1),Nreduced)
            D2=DistAngle(ReducedCoor(1,myCell2), myReduced(1),Nreduced)
          else
            do i = 1, npt
              do l = 1,3
                coor3(l,i) = coor(l,i)
              end do
            end do

            do i = 1, Ncells+2
              l = (i-1)*npt+1
            call rmsd_weight(npt,coor3(1,1),centers(1,l),D0
     &          ,.false.,mmm)
            end do

            k1 = (myCell-1)*npt + 1
            k2 = (myCell2-1)*npt + 1
            call Distance(coor3(1,1),centers(1,k1),D1)
            call Distance(coor3(1,1),centers(1,k2),D2)
          end if
            
          mindx = D1 - D2
          if (mindx .gt. upper) then
            OK = .false.
            if (verbose) write (6,*) "too close to cell2", mindx,upper
            return
          end if

          if (mindx .lt. lower) then
            OK = .false.
            if (verbose) write (6,*) "too close to cell1", mindx,lower
            return
          end if
          
          do i = 1, Ncells
            if (Nreduced .gt. 0) then
            D0 = DistAngle(ReducedCoor(1,cell(i)),myReduced(1),Nreduced)
            else
              k = (cell(i)-1)*npt
              call Distance(coor3(1,1),centers(1,k+1),D0)
            end if
            
            if (D1-D0 .gt. neigh) OK = .false.
            if (.not. OK) then
              if (verbose) write(6,*) "in a neighbor!" ,D1-D0,D2-D0
              return
            end if
          end do
        end

        subroutine TestInterfaces2(j,OK,bad,reverse)
          
          implicit none

          include "COMMON/LENGTH.BLOCK"
          include "COMMON/COORD.BLOCK"
          include "COMMON/CONNECT.BLOCK"
          include "COMMON/PATH2.BLOCK"
          include "COMMON/SEARCH.BLOCK"
          include "COMMON/VELOC.BLOCK"
          include "COMMON/DYNA.BLOCK"
          include "COMMON/OVERLAP.BLOCK"

          integer i,j,l,step, ss, bb, reached,k1,k2,k
          double precision prod, size, coor3(3,maxpt),DistAngle
          double precision mindx, D1, D2, D0
          double precision lower, upper, neigh
          logical OK,bad, reverse
          logical firstT
          data firstT/.true./

          save firstT,lower,upper,neigh

          if (firstT) then
            lower = dx - dx_approximate
            upper = dx
            neigh = dx
            firstT = .false.
          end if

          OK = .false.
          bad = .false.

          if (Nreduced .gt. 0) then

            call GetReducedCoors()
           
            D1 = DistAngle(ReducedCoor(1,myCell),myReduced(1),Nreduced)
            D2 = DistAngle(ReducedCoor(1,myCell2),myReduced(1),Nreduced)
          else
            do i = 1, npt
              do l = 1,3
                coor3(l,i) = coor(l,i)
              end do
            end do

            do i = 1, Ncells+2
              l = (i-1)*npt+1
            call rmsd_weight(npt,coor3(1,1),centers(1,l),D0
     &          ,.false.,mmm)
            end do

            k1 = (myCell-1)*npt + 1
            k2 = (myCell2-1)*npt + 1
            call Distance(coor3(1,1),centers(1,k1),D1)
            call Distance(coor3(1,1),centers(1,k2),D2) 
          end if

          mindx = D1 - D2

          if (mindx .gt. upper .and. (.not. reverse)) then
            OK = .true.
            do i = 1,npt
              do l =1,3
                Wcoor(l,i) = coor(l,i)
                Wvelo(l,i) = velo(l,i)
              end do
            end do
          end if

          if (reverse .and. mindx .lt. lower) OK = .true.

          if ( .not. reverse .and. mindx .lt. lower) then
            bad = .true.
            if(verbose) write (6,*) "falling back to cell1",mindx,lower
            return
          end if

          if ( reverse .and. mindx .gt. upper)then
            bad = .true.
            if(verbose) write (6,*) "falling back to cell2",mindx,upper
            return
          end if

           do i = 1, Ncells
             if (Nreduced.gt.0) then
             D0 =DistAngle(ReducedCoor(1,cell(i)),myReduced(1),Nreduced)
             else
               k = (cell(i)-1)*npt
               call Distance(coor3(1,1),centers(1,k+1),D0)
             end if
             
             if (D1 - D0 .gt. neigh) bad = .true.
             if (bad) then
               if (verbose) write(6,*) "falling to a neighbor!"
               return
             end if
           end do
        end



        subroutine TestInterfaces3(j,OK,bad,reverse)

          implicit none

          include "COMMON/LENGTH.BLOCK"
          include "COMMON/COORD.BLOCK"
          include "COMMON/CONNECT.BLOCK"
          include "COMMON/PATH2.BLOCK"
          include "COMMON/SEARCH.BLOCK"
          include "COMMON/VELOC.BLOCK"
          include "COMMON/DYNA.BLOCK"

          integer i,j,l,step, ss, bb, reached,k1,k2,k
          double precision prod, size
          double precision D1,D2,D0, mindx, DistAngle,coor3(3,maxpt)
          double precision lower, upper, neigh
          logical OK,bad, reverse
          logical firstT
          data firstT/.true./

          save firstT,lower,upper,neigh

          if (firstT) then
            lower = dx - dx_approximate
            upper = dx
            neigh = dx
            firstT = .false.
          end if

          OK = .false.
          bad = .false.
          reverse = .false.
 
          if (Nreduced .gt. 0) then
            call GetReducedCoors()
          
            D1=DistAngle(ReducedCoor(1,myCell), myReduced(1),Nreduced)
            D2=DistAngle(ReducedCoor(1,myCell2), myReduced(1),Nreduced)
          else
            do i = 1, npt
              do l = 1,3
                coor3(l,i) = coor(l,i)
              end do
            end do

            do i = 1, Ncells+2
              l = (i-1)*npt+1
           call rmsd_weight(npt,coor3(1,1),centers(1,l),D0,
     &          .false.,mmm)
            end do

            k1 = (myCell-1)*npt + 1
            k2 = (myCell2-1)*npt + 1
            call Distance(coor3(1,1),centers(1,k1),D1)
            call Distance(coor3(1,1),centers(1,k2),D2)
          end if

          mindx = D1 - D2

          if (mindx .lt. lower) then
            OK = .true.
            return
          end if

          if (mindx .gt. upper) then
           
            OK = .true.
            reverse = .true.
            if (verbose) write(6,*)"reversing"
            do i=1,npt
              do l=1,3
                Wvelo(l,i) = velo(l,i)
                Wcoor(l,i) = coor(l,i)
              end do
            end do  
            return
          end if
        
          do i = 1, Ncells
            if (Nreduced.gt.0) then
            D0 = DistAngle(ReducedCoor(1,cell(i)),myReduced(1),Nreduced)
            else
             k = (cell(i)-1)*npt
             call Distance(coor3(1,1),centers(1,k+1),D0)
            end if

            if (D1 - D0 .gt. neigh) then
              bad = .true.
              return
            end if
          end do
          
        end
