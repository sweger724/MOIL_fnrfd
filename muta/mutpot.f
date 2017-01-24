      subroutine mutpot(c_pr,De_lambda)
c routine for muta. if the mutant particle is labelled
c with 1, then its interaction is weighted with lambda,
c if it is labelled with 2 with 1-lambda. Double countig
c are solved in cdie.f. 
c if c_pr = 0 it scales for the first time the interactions
c if c_pr = 1 it rescales the interactions with new lambda
c if c_pr = 2 it counts U1-U2, and gives it back as De_lambda.


        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
C        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
C        include 'COMMON/SHAKE.BLOCK'
C        include 'COMMON/MSHAKE.BLOCK'
C        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
C        include 'COMMON/DYNA.BLOCK'
C        include 'COMMON/SWITCH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
C        include 'COMMON/SYMM.BLOCK'
C        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SPECL.BLOCK'
C        include 'COMMON/CONSTRAN.BLOCK'
C        include 'COMMON/TETHER.BLOCK'
C        include 'COMMON/RESTART.BLOCK'
C        include 'COMMON/SSBP.BLOCK'
C        include 'COMMON/EWALD.BLOCK'
C        include 'COMMON/METAL.BLOCK'
C        include 'COMMON/SGB.BLOCK'
C        include 'COMMON/EBALL.BLOCK'
C        include 'COMMON/RGYRCST.BLOCK'
C        include 'COMMON/MASSFAC.BLOCK'
         include 'COMMON/MUTA.BLOCK'
         include 'COMMON/PARALLEL.BLOCK'

        double precision evd,ebbe,eth,eto,ei,ene14

        double precision De_lambda,lambda_old,rde

        integer i,j,l,cntr_param
        integer c_pr,ii,k,h

        double precision  epsgm6_old(30000),epsgm12_old(30000)
        double precision  ptchg_old(30000)
        double precision  kbond_old(30000), kangl_old(30000)
        double precision  ktors1_old(30000), ktors2_old(30000)
        double precision  kimp_old(30000)


        De_lambda=0.0d0
        De_lambda2=0.0d0

        cntr_param=c_pr

        if (cntr_param.eq.0) then

C******************************************************
C* IF ONE OF THE PARTICLES IS INVOLVED IN MUTATIONS   *
C* THE EPS6, EPS12, PTCHG WILL BE SCALED WITH EITHER  *
C* L OR (1-L). DOUBLE COUNTING OF THE SCALING (INTERA-*
C* CTIONS BETWEEN TWO MUTANT PARTICLES) WILL BE SOLVED*
C* IN SUBROUTINE CDIE.F                               * 
C******************************************************

        do i=1,npt
         epsgm6_old(i)=epsgm6(i)
         epsgm12_old(i)=epsgm12(i)
         ptchg_old(i)=ptchg(i)
        enddo
 
         do j=1,npt
           if (mutaid(j).eq.1) then
           epsgm6(j)=epsgm6(j)*lambda
           epsgm12(j)=epsgm12(j)*lambda
           ptchg(j)=ptchg(j)*lambda
            else if (mutaid(j).eq.2) then
           epsgm6(j)=epsgm6(j)*(1.0d0-lambda)
           epsgm12(j)=epsgm12(j)*(1.0d0-lambda)
           ptchg(j)=ptchg(j)*(1.0d0-lambda)
           endif
          enddo

        do i=1,npt
!         write(udeb,*) i,mutaid(i),lambda,epsgm6_old(i),epsgm6(i)
!         write(udeb,*) i,mutaid(i),lambda,epsgm12_old(i),epsgm12(i)
!         write(udeb,*) i,mutaid(i),lambda,ptchg_old(i),ptchg(i)
        enddo

C CHANGE P14

c******************************************************
c* P14 ARE SCALED WITH LAMBDA OR (1-LAMBDA) DEPENDING *
c* ON MUTAID VALUE. IF TWO PARTICLES BELONG TO SAME   *
c* MUTANT SIDE CHAIN THE ENERGY MUST BE RESCALED ONLY *
c* BY L, NOT BY L^2                                   *
c******************************************************

       do i=1,totspe
          j = spec1(i)
          k = spec2(i)
!         write (udeb,*) 'P14 - i,j,k',i,j,k
!         write (udeb,*) 'lambda -',lambda
!         write (udeb,*) ' 1,2,3 ',p14(1,i),p14(2,i),p14(3,i)
          if (p14(1,i).gt.1.d-12) then
           if (mutaid(j) .eq. 1 .or. mutaid(k) .eq. 1) then
            p14(1,i)=p14(1,i)*lambda
            p14(2,i)=p14(2,i)*lambda
           else if (mutaid(j) .eq. 2 .or. mutaid(k) .eq. 2) then
            p14(1,i)=p14(1,i)*(1.0d0-lambda)
            p14(2,i)=p14(2,i)*(1.0d0-lambda)
           endif
          endif
          if (dabs(p14(3,i)).gt.1.d-12) then
           if (mutaid(j) .eq. 1 .or. mutaid(k) .eq. 1) then
            p14(3,i)=p14(3,i)*lambda
           else if (mutaid(j) .eq. 2 .or. mutaid(k) .eq. 2) then
            p14(3,i)=p14(3,i)*(1.0d0-lambda)
           endif
          endif
!         write (udeb,*) 'after l',p14(1,i),p14(2,i),p14(3,i)
        enddo   


C CHANGE BONDS!

c*******************************************************
c* IF ONE OF THE PARTICLES INVOLVED IN A BOND, OR BOTH *
c* ARE MUTANT PARTICLES, THEN RESCALE THEIR ENERGY WITH*
C* EITHER L OR (1-L). THE SAME IS DONE FOR ANGLES,     *
C* TORSIONS, IMPROPER TORSIONS                         *
c*******************************************************

!          write (udeb,*) 'bonds'
 

          do ii=1,nb
           i = ib1(ii)
           j = ib2(ii)
!           write (udeb,*) 'i=',i,'j=',j,'kbond=',kbond(ii)
           if (mutaid(i) .eq. 1 .or. mutaid(j) .eq. 1) then
              kbond(ii)=kbond(ii)*lambda
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'kbond=',kbond(ii)
           else if (mutaid(i) .eq. 2 .or. mutaid(j) .eq. 2) then
              kbond(ii)=kbond(ii)*(1.0d0-lambda)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'kbond=',kbond(ii)
           endif
          enddo

C CHANGE ANGLES!

!          write (udeb,*) 'angles'

          do ii=1,nangl
           i=iangl1(ii)
           j=iangl2(ii)
           k=iangl3(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 .or. 
     &         mutaid(k).eq.1)then
              kangl(ii)=kangl(ii)*lambda
!         write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2) then
              kangl(ii)=kangl(ii)*(1.0d0-lambda)
!         write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           endif
          enddo  
       
C CHANGE TORSIONS!

!          write (udeb,*) 'torsions'
          do ii=1,ntors
           i=itor1(ii)
           j=itor2(ii)
           k=itor3(ii)
           h=itor4(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 
     &                        .or. mutaid(k).eq.1 
     &                        .or. mutaid(h).eq.1 ) then
              ktors1(ii)=ktors1(ii)*lambda
              ktors2(ii)=ktors2(ii)*lambda
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2 .or. mutaid(h).eq.2) then
              ktors1(ii)=ktors1(ii)*(1.0d0-lambda)
              ktors2(ii)=ktors2(ii)*(1.0d0-lambda)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           endif
           enddo

C CHANGE IMPROPER TORSIONS!

!          write (udeb,*) 'Improper torsions'

          do ii=1,nimp
           i=iimp1(ii)
           j=iimp2(ii)
           k=iimp3(ii)
           h=iimp4(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 
     &                        .or. mutaid(k).eq.1 
     &                        .or. mutaid(h).eq.1 ) then
              kimp(ii)=kimp(ii)*lambda
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2 .or. mutaid(h).eq.2) then
              kimp(ii)=kimp(ii)*(1.0d0-lambda)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           endif
           enddo



 
         else  if (cntr_param.eq.1) then


          lambda_old=lambda-lambda_step

C******************************************************
C* IF ONE OF THE PARTICLES IS INVOLVED IN MUTATIONS   *
C* THE EPS6, EPS12, PTCHG WILL BE SCALED WITH EITHER  *
C* L/L_OLD OR (1-L)/(1-L_OLD). DOUBLE COUNTING OF THE *
C* SCALING (INTERACTIONS BETWEEN TWO MUTANT PARTICLES)* 
C* WILL BE SOLVED IN SUBROUTINE CDIE.F                * 
C******************************************************

         do j=1,npt
           if (mutaid(j).eq.1) then
           epsgm6(j)=epsgm6(j)*lambda/lambda_old
           epsgm12(j)=epsgm12(j)*lambda/lambda_old
           ptchg(j)=ptchg(j)*lambda/lambda_old
            else if (mutaid(j).eq.2) then
           epsgm6(j)=epsgm6(j)*(1.0d0-lambda)/
     &     (1.0d0-lambda_old)
           epsgm12(j)=epsgm12(j)*(1.0d0-lambda)/
     &     (1.0d0-lambda_old)
           ptchg(j)=ptchg(j)*(1.0d0-lambda)/
     &     (1.0d0-lambda_old)
           endif
          enddo 

C CHANGE P14

c******************************************************
c* P14 ARE SCALED WITH L/L_OLD OR (1-L)/(1-L_OLD)     *
C* DEPENDING ON MUTAID VALUE. IF TWO PARTICLES BELONG *
C* TO SAME MUTANT SIDE CHAIN THE ENERGY MUST BE       *
C* RESCALED ONLY BY L, NOT BY L^2                     *
c******************************************************

       do i=1,totspe
          j = spec1(i)
          k = spec2(i)
!         write (udeb,*) 'P14 - i,j,k',i,j,k
!         write (udeb,*) 'lambda -',lambda,lambda_old
!         write (udeb,*) ' 1,2,3 ',p14(1,i),p14(2,i),p14(3,i)
          if (p14(1,i).gt.1.d-12) then
           if (mutaid(j) .eq. 1 .or. mutaid(k) .eq. 1) then
            p14(1,i)=p14(1,i)*lambda/lambda_old
            p14(2,i)=p14(2,i)*lambda/lambda_old
           else if (mutaid(j) .eq. 2 .or. mutaid(k) .eq. 2) then
            p14(1,i)=p14(1,i)*(1.0d0-lambda)/(1.0d0-lambda_old)
            p14(2,i)=p14(2,i)*(1.0d0-lambda)/(1.0d0-lambda_old)
           endif
          endif
          if (dabs(p14(3,i)).gt.1.d-12) then
           if (mutaid(j) .eq. 1 .or. mutaid(k) .eq. 1) then
            p14(3,i)=p14(3,i)*lambda/lambda_old
           else if (mutaid(j) .eq. 2 .or. mutaid(k) .eq. 2) then
            p14(3,i)=p14(3,i)*(1.0d0-lambda)/(1.0d0-lambda_old)
           endif
          endif
!         write (udeb,*) 'after l',p14(1,i),p14(2,i),p14(3,i)
        enddo

!        stop

C CHANGE BONDS!

c*******************************************************
c* IF ONE OF THE PARTICLES INVOLVED IN A BOND, OR BOTH *
c* ARE MUTANT PARTICLES, THEN RESCALE THEIR ENERGY WITH*
C* EITHER L/L_OLD OR (1-L)/(1-L_OLD). THE SAME IS DONE * 
C* FOR ANGLES. TORSIONS, IMPROPER TORSIONS             *
c*******************************************************

!          write (udeb,*) 'bonds'

          do ii=1,nb
           i = ib1(ii)
           j = ib2(ii)
!           write (udeb,*) 'i=',i,'j=',j,'kbond=',kbond(ii)
           if (mutaid(i) .eq. 1 .or. mutaid(j) .eq. 1) then
              kbond(ii)=kbond(ii)*lambda/lambda_old
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'kbond=',kbond(ii)
           else if (mutaid(i) .eq. 2 .or. mutaid(j) .eq. 2) then
              kbond(ii)=kbond(ii)*(1.0d0-lambda)/(1.0d0-lambda_old)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'kbond=',kbond(ii)
           endif
          enddo

C CHANGE ANGLES!

C          write (udeb,*) 'angles'

          do ii=1,nangl
           i=iangl1(ii)
           j=iangl2(ii)
           k=iangl3(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 
     &                        .or. mutaid(k).eq.1)then
              kangl(ii)=kangl(ii)*lambda/lambda_old
!         write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2) then
              kangl(ii)=kangl(ii)*(1.0d0-lambda)/(1.0d0-lambda_old)
!         write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'kangl=',kangl(ii)
           endif
          enddo

C CHANGE TORSIONS!

!          write (udeb,*) 'torsions'
          do ii=1,ntors
           i=itor1(ii)
           j=itor2(ii)
           k=itor3(ii)
           h=itor4(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 
     &                        .or. mutaid(k).eq.1 
     &                        .or. mutaid(h).eq.1 ) then
              ktors1(ii)=ktors1(ii)*lambda/lambda_old
              ktors2(ii)=ktors2(ii)*lambda/lambda_old
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2 .or. mutaid(h).eq.2) then
              ktors1(ii)=ktors1(ii)*(1.0d0-lambda)/(1.0d0-lambda_old)
              ktors2(ii)=ktors2(ii)*(1.0d0-lambda)/(1.0d0-lambda_old)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'k1=',ktors1(ii),'k2=',ktors2(ii)
           endif
           enddo

C CHANGE IMPROPER TORSIONS!

          write (udeb,*) 'Improper torsions'

          do ii=1,nimp
           i=iimp1(ii)
           j=iimp2(ii)
           k=iimp3(ii)
           h=iimp4(ii)
!           write (udeb,*) 'i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           if (mutaid(i).eq.1 .or. mutaid(j).eq.1 
     &                        .or. mutaid(k).eq.1 
     &                        .or. mutaid(h).eq.1 ) then
              kimp(ii)=kimp(ii)*lambda/lambda_old
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           else if (mutaid(i).eq.2 .or. mutaid(j).eq.2 .or.
     &              mutaid(k).eq.2 .or. mutaid(h).eq.2) then
              kimp(ii)=kimp(ii)*(1.0d0-lambda)/(1.0d0-lambda_old)
!              write (udeb,*) 'NEW! - i=',i,'j=',j,'k=',k,'h=',h,
!     &                    'kimp=',kimp(ii)
           endif
           enddo

!          write (udeb,*) 'after'
          do j=1,cnt_l
!             write(udeb,*) 'j=',j,'c_pr=',c_pr
!             write(udeb,*) '   epsgm6',epsgm6(poi_l(j))
!             write(udeb,*) '   epsgm12',epsgm12(poi_l(j))
!             write(udeb,*) '   ptchg',ptchg(poi_l(j))
          enddo
!          write (udeb,*) '******'


        else if (cntr_param.eq.2) then

        tmp_e_lambda=0.0d0
        evd=0.0d0
        ebbe=0.0d0
        eth=0.0d0
        eto=0.0d0
        ei=0.0d0
        call cdie()
        write (udeb,*) '1',tmp_e_lambda
        evd=tmp_e_lambda
        edeb(1,i_lambda) = edeb(1,i_lambda) +evd
        edeb2(1,i_lambda)= edeb2(1,i_lambda)+evd**2        
        call ebond()
        ebbe=tmp_e_lambda-evd
        edeb(2,i_lambda) = edeb(2,i_lambda) +ebbe
        edeb2(2,i_lambda)= edeb2(2,i_lambda)+ebbe**2
        call etheta()
        eth=tmp_e_lambda-evd-ebbe
        edeb(3,i_lambda) = edeb(3,i_lambda) +eth
        edeb2(3,i_lambda)= edeb2(3,i_lambda)+eth**2
        call etors()
        eto=tmp_e_lambda-evd-ebbe-eth
        edeb(4,i_lambda) = edeb(4,i_lambda) +eto
        edeb2(4,i_lambda)= edeb2(4,i_lambda)+eto**2
        call eimphi()
        ei=tmp_e_lambda-evd-ebbe-eth-eto
        edeb(5,i_lambda) = edeb(5,i_lambda) +ei
        edeb2(5,i_lambda)= edeb2(5,i_lambda)+ei**2
        call ener14()
        ene14=tmp_e_lambda-evd-ebbe-eth-eto-ei
        edeb(6,i_lambda) = edeb(6,i_lambda) +ene14
        edeb2(6,i_lambda)= edeb2(6,i_lambda)+ene14**2
!        write (udeb,*) 'tmp_e_lambda,evd,ebbe,eth,eto,ei,ene14'
!        write (udeb,*) tmp_e_lambda,evd,ebbe,eth,eto,ei,ene14

!         write (udeb,*) 'k=',kbond(1)
!         write (udeb,*) 'C=',coor(1,1),coor(2,1),coor(3,1)
!         write (udeb,*) 'O=',coor(2,1),coor(2,2),coor(3,2)
         rde=((coor(1,1)-coor(1,2))**2+(coor(2,1)-coor(2,2))**2+
     &      (coor(3,1)-coor(3,2))**2)**(0.5)
!         write (udeb,*) 'r=',rde,'l=',i_lambda
!         write (udeb,*) 'k2=',kbond(1)
!         write (udeb,*) 'C2=',coor(1,1),coor(2,1),coor(3,1)
!         write (udeb,*) 'O2=',coor(2,1),coor(2,2),coor(3,2)
         rde=((coor(1,3)-coor(1,4))**2+(coor(2,3)-coor(2,4))**2+
     &      (coor(3,3)-coor(3,4))**2)**(0.5)
!         write (udeb,*) 'r2=',rde,'l=',i_lambda
!         write (udeb,*) 'req2=',req(2)
!         write (udeb,*) 'ebond=',ebbe,'l=',i_lambda

         De_lambda = tmp_e_lambda

         endif


         return
         end
                  
