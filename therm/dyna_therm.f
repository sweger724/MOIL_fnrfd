       subroutine dyna_therm(temp,tmpr,step,dt2,hlfdt,
     1  tfac,beta,partit,exp2,
     2  winda,windb,delta,lambda,npri,nstep,neqsvel,nsvel,
     3  dmass,divms,
     4  grdlx,grdly,grdlz,dpold,fact1,
     5  fact2,udata,ucrd,nwcrd,ntest,icol,newv,irand,debug,
     6  shakb,shakl,shakm,
     7  nlist,ndegf,massw)

        double precision temp,tmpr,step,tfac,beta,partit,exp2
        double precision winda,windb,delta,lambda
        double precision dt2,hlfdt
        integer iat1,iat2,iac1,iac2,k,ndegf,npri,nstep,neqsvel,nsvel,nlist
        integer ucrd,udata,nwcrd,ntest,icol,newv,irand
        logical debug,massw


        logical shakl,shakb,shakm
c
c temp   -  assigned temperature
c tmpr   -  actually calculated temperature (from kinetic energy)
c               using function hot_umbr
c step   -  time integration step (in AMKA units)
c tfac   -  conversion factor for time from PS to AMKA units
c winda    initial number along LAMBDA
c windb    final number  along LAMBDA
c delta    step along LAMBDA( user can calculate free energy
c           difference by TPM in the reverse direction )
c npri   -  print some useful(?) data each NPRI steps
c nstep  -  number of integration steps
c nsvel  -  number of steps between scaling of the velocities
c ucrd   -  unit number of file on which the coordinates (path format)
c               are written.
c udata  -  write info on the run on unit UDATA
c nwcrd  -  write coordinates on ucrd each NWCRD steps
c ntest  -  test the constraints each NTEST steps
c icol   -  use ICOL steps for equilibratio (before collecting data)
c newv   -  assign new velocities each NEWV steps
c irand  - a seed for a random nmber generator
c debug  - if .true. print a LOT of debugging info
c
c
        double precision dmass(*)
        double precision divms(*),grdlx(*),grdly(*),grdlz(*)
        double precision dpold(3,*),fact1(*),fact2(*)
c


        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONNECT2.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SYMM.BLOCK'
c
c local
        integer i,j,npt3,smlp,idyna,istp,lend
        double precision eorig,escond,enetot,hot_umbr
        double precision kinet,vfac
        double precision  diff,grad_cons,curpartit,curexp2

        double precision  cur_cons,par_cons,cur_energy,cur_entr
        double precision  par_energy,cur_grad,par_grad,cur_mix,par_mix

        double precision xforce1,yforce1,zforce1
        double precision xforce2,yforce2,zforce2
        double precision tmp,epsilon


        double precision xu(maxpt),yu(maxpt),zu(maxpt)
        double precision xc(maxpt),yc(maxpt),zc(maxpt)
        double precision xc_ori(maxpt),yc_ori(maxpt),zc_ori(maxpt)
        double precision xc_rem(maxpt),yc_rem(maxpt),zc_rem(maxpt)
        double precision xc_ini(maxpt),yc_ini(maxpt),zc_ini(maxpt)
        double precision dxc_ori(maxpt),dyc_ori(maxpt),dzc_ori(maxpt)
        double precision vxold(maxpt),vyold(maxpt),vzold(maxpt)
        double precision vxc_ori(maxpt),vyc_ori(maxpt),vzc_ori(maxpt)

c@
        epsilon = 1.d-5

        npt3  = npt*3
        smlp  = 1

        par_cons=0.d0
        par_energy=0.d0
        par_grad=0.d0
        par_mix=0.d0


c
c copy comp set of coordinates to current X,Y,Z
c

                do 2 i=1,npt
                coor(1,i)=coor2(1,i)
                coor(2,i)=coor2(2,i)
                coor(3,i)=coor2(3,i)
2               continue



          if (esymyes) call squeeze()
          call nbondm()
          if (esymyes) call syminit()

c
c calculate the energy and the derivatives
c for the original point along the lambda path
c

          call eforce()

          eorig=e_total

          do 340 i=1,npt
          dxc_ori(i)=dpot(1,i)
          dyc_ori(i)=dpot(2,i)
          dzc_ori(i)=dpot(3,i)
340       continue


          do 341 i=1,npt
          xc_ori(i)=coor(1,i)
          yc_ori(i)=coor(2,i)
          zc_ori(i)=coor(3,i)
341       continue



        write(udata,1001)  eorig
1001    format(1x,'first energy at the original point ',f10.4)
        write(udata,1002) lambda
1002    format(1x,'energy sampled at the position ',f10.4)





c
c calculate the energy for the remote point along the lambda path
c
          lambda=lambda+delta

          do 355 i=1,npt
          ptms(i)=(1.d0-lambda)*ptmsA(i)+ lambda*ptmsB(i)
         dmass(i) =(1.d0-lambda)*ptmsA(i)+lambda*ptmsB(i)
                divms(i) = 1.d0/dmass(i)
                fact1(i) = dt2*divms(i)
                fact2(i) = hlfdt*divms(i)
          ptwei(i)=(1.d0-lambda)*ptweiA(i)+ lambda*ptweiB(i)
          ptchg(i)=(1.d0-lambda)*ptchgA(i)+ lambda*ptchgB(i)
          epsgm6(i)=(1.d0-lambda)*epsgm6A(i)+
     1              lambda*epsgm6B(i)
          epsgm12(i)=(1.d0-lambda)*epsgm12A(i)+
     1              lambda*epsgm12B(i)
355       continue
          do 356 i=1,nb
          kbond(i)=(1.d0-lambda)*kbondA(i)+ lambda*kbondB(i)
          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)
356       continue
          do 357 i=1,nangl
          kangl(i)=(1.d0-lambda)*kanglA(i)+ lambda*kanglB(i)
          angleq(i)=(1.d0-lambda)*angleqA(i)+ lambda*angleqB(i)
357       continue
          do 358 i=1,ntors
          period(i)=(1.d0-lambda)*periodA(i)+lambda*periodB(i)
          ktors1(i)=(1.d0-lambda)*ktors1A(i)+ lambda*ktors1B(i)
          ktors2(i)=(1.d0-lambda)*ktors2A(i)+ lambda*ktors2B(i)
          ktors3(i)=(1.d0-lambda)*ktors3A(i)+ lambda*ktors3B(i)
          phase1(i) =(1.d0-lambda)*phase1A(i) + lambda*phase1B(i)
          phase2(i) =(1.d0-lambda)*phase2A(i) + lambda*phase2B(i)
          phase3(i) =(1.d0-lambda)*phase3A(i) + lambda*phase3B(i)
358       continue
          do 359 i=1,nimp
          kimp(i)=(1.d0-lambda)*kimpA(i)+ lambda*kimpB(i)
          impeq(i)=(1.d0-lambda)*impeqA(i)+ lambda*impeqB(i)
359       continue


c
c calculation energy at the remote point along lambda corrdinate
c


          call eforce()
          escond=e_total


          write(udata,1003) escond
1003      format(1x,' first energy at the remote point ',f10.4)
          write(udata,1004) lambda
1004      format(1x,' energy sampeld at the position ',f10.4)


          do 451 i=1,npt
          xc_rem(i)=coor(1,i)
          yc_rem(i)=coor(2,i)
          zc_rem(i)=coor(3,i)
451       continue





c
c go back to the original value of lambda
c

          lambda=lambda-delta


          do 955 i=1,npt
          ptms(i)=(1.d0-lambda)*ptmsA(i)+ lambda*ptmsB(i)
         dmass(i) =(1.d0-lambda)*ptmsA(i)+lambda*ptmsB(i)
                divms(i) = 1.d0/dmass(i)
                fact1(i) = dt2*divms(i)
                fact2(i) = hlfdt*divms(i)
          ptwei(i)=(1.d0-lambda)*ptweiA(i)+ lambda*ptweiB(i)
          ptchg(i)=(1.d0-lambda)*ptchgA(i)+ lambda*ptchgB(i)
          epsgm6(i)=(1.d0-lambda)*epsgm6A(i)+
     1              lambda*epsgm6B(i)
          epsgm12(i)=(1.d0-lambda)*epsgm12A(i)+
     1              lambda*epsgm12B(i)
955       continue
          do 956 i=1,nb
          kbond(i)=(1.d0-lambda)*kbondA(i)+ lambda*kbondB(i)
          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)
956       continue
          do 957 i=1,nangl
          kangl(i)=(1.d0-lambda)*kanglA(i)+ lambda*kanglB(i)
          angleq(i)=(1.d0-lambda)*angleqA(i)+ lambda*angleqB(i)
957       continue
          do 958 i=1,ntors
          period(i)=(1.d0-lambda)*periodA(i)+lambda*periodB(i)
          ktors1(i)=(1.d0-lambda)*ktors1A(i)+ lambda*ktors1B(i)
          ktors2(i)=(1.d0-lambda)*ktors2A(i)+ lambda*ktors2B(i)
          ktors3(i)=(1.d0-lambda)*ktors3A(i)+ lambda*ktors3B(i)
          phase1(i) =(1.d0-lambda)*phase1A(i) + lambda*phase1B(i)
          phase2(i) =(1.d0-lambda)*phase2A(i) + lambda*phase2B(i)
          phase3(i) =(1.d0-lambda)*phase3A(i) + lambda*phase3B(i)
958       continue
          do 959 i=1,nimp
          kimp(i)=(1.d0-lambda)*kimpA(i)+ lambda*kimpB(i)
          impeq(i)=(1.d0-lambda)*impeqA(i)+ lambda*impeqB(i)
959       continue

        diff=escond-eorig
        write(stdo,*)  diff
1005    format(1x,'first energy difference',f10.4)
        write(udata,*) '***************************************** '
        write(udata,*) '***  END OF TRAJECTORY INITIALIZATION *** '
        write(udata,*) '***************************************** '



        if (nlist .ne. 0 ) smlp=nlist
c
c Start integration loop
c

        do  100  idyna = 1,nstep,smlp
c
c do internal update each inbfrq for the non-bonded AND the images
c (if applicable). A single non-bonded list is generated according
c to the first structure.
c


         if (nlist.ne.0) then
          call nbondm()
         end if

         lend=min(nstep,smlp-1+idyna)
         do 100 istp = idyna,lend

         grad_cons=0.d0


c
c copy old potential derivatives to temporary data
c

                do 650  i=1,npt
                        dpold(1,i)=dxc_ori(i)
                        dpold(2,i)=dyc_ori(i)
                        dpold(3,i)=dzc_ori(i)
                        vxold(i)=velo(1,i)
                        vyold(i)=velo(2,i)
                        vzold(i)=velo(3,i)
650             continue


c ************************************************
c       if ( shakb .or. shakl ) then
c         nshak=0
c         call shakinit(shakl,shakb,shakm,epsilon)
c         end if

c ************************************************




          if (shakl.or.shakb .and. nshak.gt.0) then
           do 652 k=1,nshak
            iat1 = ishak1(k)
            iat2 = ishak2(k)
            cooref(1,k) = xc_ori(iat1) - xc_ori(iat2)
            cooref(2,k) = yc_ori(iat1) - yc_ori(iat2)
            cooref(3,k) = zc_ori(iat1) - zc_ori(iat2)
652        continue
          end if
c
c calculate a "free" (unconstrained) step
c
                 do 653 i=1,npt
                 xu(i) = xc_ori(i) + vxold(i)*step - dpold(1,i)*fact1(i)
                 yu(i) = yc_ori(i) + vyold(i)*step - dpold(2,i)*fact1(i)
                 zu(i) = zc_ori(i) + vzold(i)*step - dpold(3,i)*fact1(i)
653              continue

               do 564 i=1,npt
               xc_ini(i)=xc_ori(i)
               yc_ini(i)=yc_ori(i)
               zc_ini(i)=zc_ori(i)
564            continue



               do 353 i=1,npt

               velo(1,i) = vxold(i)*step-fact1(i)*dpold(1,i)
               velo(2,i) = vyold(i)*step-fact1(i)*dpold(2,i)
               velo(3,i) = vzold(i)*step-fact1(i)*dpold(3,i)


353            continue


          if (shakl.or.shakb .and. nshak.gt.0) then

           call shakept(itershak)

           end if

               tmp=1.d0/step

               do 354 i=1,npt

               coor(1,i) = xc_ori(i)+velo(1,i)
               coor(2,i) = yc_ori(i)+velo(2,i)
               coor(3,i) = zc_ori(i)+velo(3,i)

               velo(1,i)=velo(1,i)*tmp
               velo(2,i)=velo(2,i)*tmp
               velo(3,i)=velo(3,i)*tmp

354            continue





               do 655 i=1,npt
               xc_ori(i)=coor(1,i)
               yc_ori(i)=coor(2,i)
               zc_ori(i)=coor(3,i)
655            continue


               do 565 i=1,npt
               vxc_ori(i)=velo(1,i)
               vyc_ori(i)=velo(2,i)
               vzc_ori(i)=velo(3,i)
565            continue

               call eforce()
               eorig=e_total

               do 656 i=1,npt
               dxc_ori(i)=dpot(1,i)
               dyc_ori(i)=dpot(2,i)
               dzc_ori(i)=dpot(3,i)
656            continue



              do 3060 i=1,nb
              iac1=ib1(i)
              iac2=ib2(i)

          if(reqB(i).eq.reqA(i)) goto 3060


          coor(1,iac1)=xu(iac1)
          coor(2,iac1)=yu(iac1)
          coor(3,iac1)=zu(iac1)

          coor(1,iac2)=xu(iac2)
          coor(2,iac2)=yu(iac2)
          coor(3,iac2)=zu(iac2)

          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)

c ************************************************
c       if ( shakb .or. shakl ) then
c         nshak=0
c         call shakinit(shakl,shakb,shakm,epsilon)
c         end if

c ************************************************



          call shakecons(epshak,itershak,iac1,iac2)



          xforce1=
     1    (dmass(iac1)/(step*step))*(coor(1,iac1)-xu(iac1))

          yforce1=
     1    (dmass(iac1)/(step*step))*(coor(2,iac1)-yu(iac1))
          zforce1=
     1    (dmass(iac1)/(step*step))*(coor(3,iac1)-zu(iac1))

          xforce2=
     1    (dmass(iac2)/(step*step))*(coor(1,iac2)-xu(iac2))
          yforce2=
     1    (dmass(iac2)/(step*step))*(coor(2,iac2)-yu(iac2))
          zforce2=
     1    (dmass(iac2)/(step*step))*(coor(3,iac2)-zu(iac2))


          xc(iac1)=coor(1,iac1)
          yc(iac1)=coor(2,iac1)
          zc(iac1)=coor(3,iac1)


          xc(iac2)=coor(1,iac2)
          yc(iac2)=coor(2,iac2)
          zc(iac2)=coor(3,iac2)



          lambda=lambda+delta

          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)


c ************************************************
c       if ( shakb .or. shakl ) then
c         nshak=0
c         call shakinit(shakl,shakb,shakm,epsilon)
c         end if

c ************************************************


          coor(1,iac1)=xu(iac1)
          coor(2,iac1)=yu(iac1)
          coor(3,iac1)=zu(iac1)


          coor(1,iac2)=xu(iac2)
          coor(2,iac2)=yu(iac2)
          coor(3,iac2)=zu(iac2)



          call shakecons(epshak,itershak,iac1,iac2)


          xforce1=xforce1*(coor(1,iac1)-xc(iac1))/
     1          ((reqB(i)-reqA(i))*delta)
          yforce1=yforce1*(coor(2,iac1)-yc(iac1))/
     1          ((reqB(i)-reqA(i))*delta)
          zforce1=zforce1*(coor(3,iac1)-zc(iac1))/
     1          ((reqB(i)-reqA(i))*delta)


          xforce2=xforce2*(coor(1,iac2)-xc(iac2))/
     1          ((reqB(i)-reqA(i))*delta)
          yforce2=yforce2*(coor(2,iac2)-yc(iac2))/
     1          ((reqB(i)-reqA(i))*delta)
          zforce2=zforce2*(coor(3,iac2)-zc(iac2))/
     1          ((reqB(i)-reqA(i))*delta)

          grad_cons=grad_cons-(xforce1+xforce2+
     1     yforce1+yforce2+zforce1+zforce2)*(reqB(i)-reqA(i))

          lambda=lambda-delta


          coor(1,iac1)=xu(iac1)
          coor(2,iac1)=yu(iac1)
          coor(3,iac1)=zu(iac1)

          coor(1,iac2)=xu(iac2)
          coor(2,iac2)=yu(iac2)
          coor(3,iac2)=zu(iac2)

3060      continue





c
c find energy at the remote point along lambda coordinate
c
            lambda=lambda+delta

          do 25 i=1,npt
          ptms(i)=(1.d0-lambda)*ptmsA(i)+ lambda*ptmsB(i)
                dmass(i) =(1.d0-lambda)*ptmsA(i)+lambda*ptmsB(i)
                divms(i) = 1.d0/dmass(i)
                fact1(i) = dt2*divms(i)
                fact2(i) = hlfdt*divms(i)
          ptwei(i)=(1.d0-lambda)*ptweiA(i)+ lambda*ptweiB(i)
          ptchg(i)=(1.d0-lambda)*ptchgA(i)+ lambda*ptchgB(i)
          epsgm6(i)=(1.d0-lambda)*epsgm6A(i)+
     1              lambda*epsgm6B(i)
          epsgm12(i)=(1.d0-lambda)*epsgm12A(i)+
     1               lambda*epsgm12B(i)
25         continue
          do 26 i=1,nb
          kbond(i)=(1.d0-lambda)*kbondA(i)+ lambda*kbondB(i)
          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)
26         continue
          do 27 i=1,nangl
          kangl(i)=(1.d0-lambda)*kanglA(i)+ lambda*kanglB(i)
          angleq(i)=(1.d0-lambda)*angleqA(i)+ lambda*angleqB(i)
27         continue
          do 28 i=1,ntors
          period(i)=(1.d0-lambda)*periodA(i)+lambda*periodB(i)
          ktors1(i)=(1.d0-lambda)*ktors1A(i)+ lambda*ktors1B(i)
          ktors2(i)=(1.d0-lambda)*ktors2A(i)+ lambda*ktors2B(i)
          ktors3(i)=(1.d0-lambda)*ktors3A(i)+ lambda*ktors3B(i)
          phase1(i) =(1.d0-lambda)*phase1A(i) + lambda*phase1B(i)
          phase2(i) =(1.d0-lambda)*phase2A(i) + lambda*phase2B(i)
          phase3(i) =(1.d0-lambda)*phase3A(i) + lambda*phase3B(i)
28         continue
          do 29 i=1,nimp
          kimp(i)=(1.d0-lambda)*kimpA(i)+ lambda*kimpB(i)
          impeq(i)=(1.d0-lambda)*impeqA(i)+ lambda*impeqB(i)
29         continue



          do 995 i=1,npt
          coor(1,i)=xu(i)
          coor(2,i)=yu(i)
          coor(3,i)=zu(i)
995       continue


c ************************************************
c       if ( shakb .or. shakl ) then
c         nshak=0
c         call shakinit(shakl,shakb,shakm,epsilon)
c         end if

c ************************************************


               do 453 i=1,npt

               velo(1,i) = vxold(i)*step-fact1(i)*dpold(1,i)
               velo(2,i) = vyold(i)*step-fact1(i)*dpold(2,i)
               velo(3,i) = vzold(i)*step-fact1(i)*dpold(3,i)


453            continue


          if (shakl.or.shakb .and. nshak.gt.0) then

           call shakept(itershak)

           end if


               do 454 i=1,npt

               coor(1,i) = xc_ini(i)+velo(1,i)
               coor(2,i) = yc_ini(i)+velo(2,i)
               coor(3,i) = zc_ini(i)+velo(3,i)

               velo(1,i)=velo(1,i)*tmp
               velo(2,i)=velo(2,i)*tmp
               velo(3,i)=velo(3,i)*tmp



454            continue






               do 755 i=1,npt
               xc_rem(i)=coor(1,i)
               yc_rem(i)=coor(2,i)
               zc_rem(i)=coor(3,i)
755            continue

               call eforce()
               escond=e_total

c
c go back to the original point
c
          lambda=lambda-delta


          do 255 i=1,npt
          ptms(i)=(1.d0-lambda)*ptmsA(i)+ lambda*ptmsB(i)
         dmass(i) =(1.d0-lambda)*ptmsA(i)+lambda*ptmsB(i)
                divms(i) = 1.d0/dmass(i)
                fact1(i) = dt2*divms(i)
                fact2(i) = hlfdt*divms(i)
          ptwei(i)=(1.d0-lambda)*ptweiA(i)+ lambda*ptweiB(i)
          ptchg(i)=(1.d0-lambda)*ptchgA(i)+ lambda*ptchgB(i)
          epsgm6(i)=(1.d0-lambda)*epsgm6A(i)+
     1              lambda*epsgm6B(i)
          epsgm12(i)=(1.d0-lambda)*epsgm12A(i)+
     1              lambda*epsgm12B(i)
255       continue
          do 256 i=1,nb
          kbond(i)=(1.d0-lambda)*kbondA(i)+ lambda*kbondB(i)
          req(i)=(1.d0-lambda)*reqA(i)+ lambda*reqB(i)
256       continue
          do 257 i=1,nangl
          kangl(i)=(1.d0-lambda)*kanglA(i)+ lambda*kanglB(i)
          angleq(i)=(1.d0-lambda)*angleqA(i)+ lambda*angleqB(i)
257       continue
          do 258 i=1,ntors
          period(i)=(1.d0-lambda)*periodA(i)+lambda*periodB(i)
          ktors1(i)=(1.d0-lambda)*ktors1A(i)+ lambda*ktors1B(i)
          ktors2(i)=(1.d0-lambda)*ktors2A(i)+ lambda*ktors2B(i)
          ktors3(i)=(1.d0-lambda)*ktors3A(i)+ lambda*ktors3B(i)
          phase1(i) =(1.d0-lambda)*phase1A(i) + lambda*phase1B(i)
          phase2(i) =(1.d0-lambda)*phase2A(i) + lambda*phase2B(i)
          phase3(i) =(1.d0-lambda)*phase3A(i) + lambda*phase3B(i)
258       continue
          do 259 i=1,nimp
          kimp(i)=(1.d0-lambda)*kimpA(i)+ lambda*kimpB(i)
          impeq(i)=(1.d0-lambda)*impeqA(i)+ lambda*impeqB(i)
259       continue

c
c go back to the original x,y,z
c

          do 999 i=1,npt
          coor(1,i)=xc_ori(i)
          coor(2,i)=yc_ori(i)
          coor(3,i)=zc_ori(i)
999       continue


               do 997 i=1,npt
               velo(1,i)=vxc_ori(i)
               velo(2,i)=vyc_ori(i)
               velo(3,i)=vzc_ori(i)
997            continue

c
c calculate "free" velocities
c
                 do 675 i=1,npt
                  velo(1,i) = velo(1,i) - fact2(i)*dxc_ori(i)
                  velo(2,i) = velo(2,i) - fact2(i)*dyc_ori(i)
                  velo(3,i) = velo(3,i) - fact2(i)*dzc_ori(i)
675              continue



c
c correct the velocities to satisfy the constraints
c

          if (shakl.or.shakb .and. nshak.gt.0) then
           do 552 k=1,nshak
            iat1 = ishak1(k)
            iat2 = ishak2(k)
            cooref(1,k) = xc_ori(iat1) - xc_ori(iat2)
            cooref(2,k) = yc_ori(iat1) - yc_ori(iat2)
            cooref(3,k) = zc_ori(iat1) - zc_ori(iat2)
552        continue
          end if



          if (shakl.or.shakb .and. nshak.gt.0) then
           call shakevl(epshakv,itershak)
          end if





c
c check conservation of energy
c
        if (istp/ntest*ntest .eq. istp ) then
                kinet = 0.d0
                do 8 i = 1,npt
                        kinet = kinet + dmass(i)*(velo(1,i)*velo(1,i)+
     1          velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
8               continue
                kinet = kinet/2.d0
                enetot = eorig + kinet
        write(udata,*)'**** checking energy conservation **** '
        write(udata,1006)istp,enetot
1006    format(1x,' at step ',i7,' total energy is ',e10.2)
        write(udata,1007)eorig,kinet
1007    format(1x,' tot. pot. ',e10.2,' total kin. ',e10.2)
        end if
c
c FIRST OPTION: EQILIBRATION PERIOD
c if the run is normal MD, do the following:
c if the temperature is more than 20 degrees different from the
c assigned temperature: scale the velocities

         if((neqsvel.ne.0).and.(istp.le.icol))then
          if  (istp/neqsvel*neqsvel.eq.istp) then
                tmpr = hot_umbr(velo,dmass,npt,ndegf,massw)
                if ((dabs(tmpr-temp).gt.20.d0))then
                  vfac=dsqrt(temp/tmpr)
                  write(udata,1008)istp
1008              format(1x,' scaling velocities at step ',i7)
        write(udata,1009)tmpr,temp,vfac
1009    format(1x,' current temperature is ',f8.2,
     1   ' Desired temperature is ',f8.2,
     2   ' Scaing factor = ',f8.2)
                  do 9 i=1,npt
                        velo(1,i)=velo(1,i)*vfac
                        velo(2,i)=velo(2,i)*vfac
                        velo(3,i)=velo(3,i)*vfac
9                 continue
         end if
        end if
       end if
c SECOND  OPTION: AFTER EQILIBRATION PERIOD
c if the run is normal MD, do the following:
c if the temperature is more than 20 degrees different from the
c assigned temperature: scale the velocities
         if ((nsvel.ne.0).and.(istp.gt.icol)) then
          if  (istp/nsvel*nsvel.eq.istp) then
                tmpr = hot_umbr(velo,dmass,npt,ndegf,massw)
                if ((dabs(tmpr-temp).gt.20.d0))then
                  vfac=dsqrt(temp/tmpr)
                  write(udata,1008)istp
                  write(udata,1009)tmpr,temp,vfac
                  do 10 i=1,npt
                        velo(1,i)=velo(1,i)*vfac
                        velo(2,i)=velo(2,i)*vfac
                        velo(3,i)=velo(3,i)*vfac
10                continue
         end if
        end if
       end if
c
c the time to reassign velocities?!
c
        if (istp/newv*newv.eq.istp) then

         write(udata,1010)' *** assigning new velocities at step ',istp
1010     format(1x,' *** Assigning new vlocities at step ',i7)
c
c set velocities - initial conditions
c
         call velinit(temp,1,0)

        if (.not.esymyes) then
         call prbm(velo,coor,dmass,grdlx,grdly,grdlz,npt)
        end if






         if (debug) then
        write(udata,*)' reassigning velocities at ',istp,' irand ',irand
       write(udata,103)((velo(j,i),j=1,3),i=1,npt)
103       format(1x,8(f9.4,1x),/)
         end if
        end if
c
c if after thermalization collect data
c


        if (istp.gt.icol) then
         par_energy=par_energy+eorig
         diff=escond-eorig
         partit=partit+diff
         par_grad=par_grad+diff/delta
         par_mix=par_mix+eorig*(diff/delta)
         par_cons=par_cons+grad_cons
         exp2=exp2+diff*diff
        if (istp/npri*npri.eq.istp) then
         curpartit=partit/(istp-icol)
         cur_energy=par_energy/(istp-icol)
         cur_grad=par_grad/(istp-icol)
         cur_mix=par_mix/(istp-icol)
         cur_cons=par_cons/(istp-icol)
         cur_entr=cur_energy*cur_grad-cur_mix
         curexp2=exp2/(istp-icol)
         curexp2= dsqrt((curexp2-curpartit*curpartit)/(istp-icol))
         write(udata,1011)curpartit,cur_entr,cur_cons,curexp2
1011     format(1x,' partition f / entropy / cons / variance ',
     1      4(1x,e9.4))
        end if
        end if


c
c printing out some useful information
c
        if(istp/npri*npri.eq.istp) then
        kinet=0.d0
        do 88 i=1,npt
        kinet=kinet+dmass(i)*(velo(1,i)*velo(1,i)+
     1          velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
88      continue
        kinet=kinet/2.d0
        enetot=e_total+kinet
        write(udata,*) ' **** PRINTING  INFO  **** '
        write(udata,1012) escond,eorig
1012    format(1x,' Energy at he remote / original point ',
     1     2(1x,e10.3))
        write(udata,241) kinet,e_total,escond-eorig,enetot
241   format(1x,' Kinetic Energy ',e11.4,/,1x,' Potential Energy (1)'
     1  ,e11.4,/,1x,' Potential Energy Difference (2)-(1) '
     2  ,f9.3,/,1x,' Total Energy ' ,e11.4)
                tmpr = hot_umbr(velo,dmass,npt,ndegf,massw)
        write(udata,242) tmpr
242     format(1x, 'current temperature is ',f9.3)
        end if
c
c each nwcrd steps write coordinates in path format on ucrd
c
        if (debug) write(*,*) ' istp nwcrd icol ' , istp,nwcrd,icol
        if (istp/nwcrd*nwcrd.eq.istp) then
        write(ucrd)e_total,((coor(j,i),i=1,npt),j=1,3)
        end if
100     continue
c
c ************ END INTEGRATION LOOP ***********************
c

        partit=partit/(nstep-icol)
        exp2=exp2/(nstep-icol)
        exp2=dsqrt((exp2-partit*partit)/(nstep-icol))
c *********************************************************
c PARTITION FUCTION AND FREE ENERGY DIFFERNCE ARE DONE
c *********************************************************

        return
        end


