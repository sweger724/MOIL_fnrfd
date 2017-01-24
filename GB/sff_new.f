c mysff_f.f
c This file also contains some debugging subroutines: 
c search 'egb_all_pair','egb_self', 'check_derivative' etc

      SUBROUTINE egb_calc_pair_f(coori,coorj,qi,qj,ri,rj,epol,dpoti,
     $     dpotj,moilvectorouti,moilvectoroutj,moilvectori,moilvectorj,
     $     massfaci,massfacj,diel_ext,ifsndbool,ii,jj)

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

C - - - arg types - - -

      integer  ifsndbool,ii,jj
      double precision 
     $     qi,qj,ri,rj,massfaci,massfacj, diel_ext, 
     &     coori(*),coorj(*),dpoti(*),
     $     dpotj(*),epol,moilvectorouti(*),moilvectoroutj(*),
     $     moilvectori(*),moilvectorj(*) 

C - - - local declarations - - -

      integer  kl,l,zero, i, j, k
      double precision
     $     xij,yij,zij,r2,rb2,efac,qiqj,fgbi,fgbk,expmkf,dielfac,e,
     $     de,dedx,dedy,dedz,temp4,temp5,temp6,lphi,lpsi,A1,A2,B1,B2,C1,
     $     C2,beta,betasqrt,beta1p5,beta2,beta2p5,kkonst,bit1,bit2,
     $     derivbit,deltabit,part1,term1,term2,term3,term4,rijk,rijl,
     $     mass,betaderivi,betaderivj,psiderivi,psiderivj,
     $     temp7, xi,yi,zi,qi2h, qid2h, ri1i, daix, daiy, daiz

      double precision KSCALE, BOFFSET
      parameter (KSCALE=0.73D0, BOFFSET=0.09D0)
      double precision multfact, q(MAXPT), 
     &                 dax(MAXPT), day(MAXPT), daz(MAXPT)      

      double precision rgbmaxpsmax2, dij1i, dij2i,dij3i, dij,
     &                 temp1,sj,sj2, datmp, rgbmax2, rgbmax1i,
     &                 rgbmax2i, thi, dumbo, tmpsd

	double precision eps_gb

C - - - begin - - -
   

      xij = coori(1) - coorj(1)
      yij = coori(2) - coorj(2)
      zij = coori(3) - coorj(3)

      r2 = xij*xij + yij*yij + zij*zij

      qiqj = qi * qj
      rb2 = ri * rj


      efac = DEXP(-r2/(4.0D0*rb2))
      fgbi = 1.0D0/DSQRT(r2 + rb2*efac)
      fgbk = -(kappa)* KSCALE /fgbi

      expmkf = DEXP( fgbk )/(diel_ext)
      dielfac = 1.0D0 - expmkf


      e = -qiqj * dielfac * fgbi
      epol = epol + e

      temp4 = fgbi*fgbi*fgbi
      temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf)
      de =  temp6 * (1.0d0 - 0.25d0 * efac)

      temp5 = 0.5d0 * efac * temp6 * (rb2 + 0.25d0 * r2)
      sumdeijda(ii) = sumdeijda(ii) + ri * temp5
      sumdeijda(jj) = sumdeijda(jj) + rj * temp5

c ---------------------------------------------------------------------
c Note: The derivative of Eij w.r.t. ri is actually equal to 
c 'temp5/ri'. Here, it is pre-multiplied by a factor (ri* ri),
c which is presented in the formula of the derivative of ri
c w.r.t. cartersian coordinates
c --------------------------------------------------------------------- 

c      write(*,*) "DEBUG> sumdeijda(i), (j), epol",
c     &       ii, jj,  sumdeijda(ii), sumdeijda(jj), epol

      dedx = de * xij
      dedy = de * yij
      dedz = de * zij

C Dont modify the derivatives if this is a second derivative call

      if (ifsndbool.eq.0) then
         dpoti(1) = dpoti(1) + dedx
         dpoti(2) = dpoti(2) + dedy
         dpoti(3) = dpoti(3) + dedz
         dpotj(1) = dpotj(1) - dedx
         dpotj(2) = dpotj(2) - dedy
         dpotj(3) = dpotj(3) - dedz
      end if


      IF (.not.(ifsndbool.eq.0)) THEN

         lphi = expmkf
         lpsi = efac
         A1 = -1.0D0/(2.0D0 *ri * rj) * lpsi
         B1 = 2.0D0 - 0.5D0 * lpsi
         beta = r2 + rb2*efac

         betasqrt = DSQRT(beta)
         beta1p5 = beta * betasqrt
         beta2 = beta * beta
         beta2p5 = beta2*betasqrt
         kkonst = (kappa)* 0.73D0
         eps_gb = dielfac

         bit1 = (eps_gb/beta1p5 - kkonst /beta * lphi)
         derivbit = + 0.5D0 * kkonst / beta2 * lphi  
     $        + kkonst / beta2 * lphi  
     $        - 1.5D0 * eps_gb / beta2p5  
     $        + 0.5D0 * kkonst * kkonst / beta1p5 * lphi


         bit2 = (1.0D0-0.25D0 * lpsi)
         part1 = qiqj*bit1 * bit2
         mass = massfaci * massfacj

         do  kl = 1,3
            rijk = coori(kl) - coorj(kl)
            betaderivi = B1 * rijk
            betaderivj = - B1 * rijk
            psiderivi = A1 * rijk
            psiderivj = -A1 * rijk
            
            do l = 1,3
               rijl = coori(l) - coorj(l)
               moilvectorouti(kl) = moilvectorouti(kl) +( qiqj )*  
     $              (derivbit * betaderivi * 
     $              (1.0D0 - 0.25D0 * lpsi) * rijl +  
     $              bit1 * psiderivi * (-0.25D0 * rijl)  
     $              ) * massfaci * massfaci * moilvectori(l)

               
               moilvectoroutj(kl) = moilvectoroutj(kl) +( qiqj )*  
     $              (derivbit * betaderivj * 
     $              (1.0D0 - 0.25D0 * lpsi) * rijl + 
     $              bit1 * psiderivj * (-0.25D0 * rijl)  
     $              ) * massfaci * massfacj * moilvectori(l)
               

               moilvectorouti(kl) = moilvectorouti(kl) -( qiqj )*  
     $              (derivbit * betaderivi * 
     $              (1.0D0 - 0.25D0 * lpsi) * rijl +  
     $              bit1 * psiderivi * (-0.25D0 * rijl)  
     $              ) * massfaci * massfacj * moilvectorj(l)

               moilvectoroutj(kl) = moilvectoroutj(kl) -( qiqj )*  
     $              (derivbit * betaderivj * 
     $              (1.0D0 - 0.25D0 * lpsi) * rijl +  
     $              bit1 * psiderivj * (-0.25D0 * rijl)  
     $              ) * massfacj*massfacj * moilvectorj(l)
               
               if (kl.eq.l) then 
                  moilvectorouti(kl) = moilvectorouti(kl) +
     $                 ( part1 * massfaci*massfaci )*  
     $                 moilvectori(l)

                  moilvectorouti(kl) = moilvectorouti(kl) - 
     $                 ( part1 * massfaci*massfacj )*  
     $                 moilvectorj(l)
                  
                  moilvectoroutj(kl) = moilvectoroutj(kl) +
     $                 ( -part1 * massfaci*massfacj )*  
     $                 moilvectori(l)
                  
                  moilvectoroutj(kl) = moilvectoroutj(kl) -
     $                 ( -part1 * massfacj*massfacj )*  
     $                 moilvectorj(l)
               end if  
            end do 
          end do 
        end if  
      return 
      END 

C --------------------------------------------------

      SUBROUTINE egb_deriv_rad(natom,qmoil,diel_ext,epol)

c YS     Computing derivatives of effective radius w.r.t coord
c        plus the "diagonal" term
c        Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
c 6/9/06
      include 'COMMON/CONNECT.BLOCK'

c Parameters
      integer natom
      double precision qmoil(*), diel_ext,epol

c Local variables
cRE
      integer i,j,k,jbeg1,jend1,jbeg3,jend3,jbegex,jendex

      double precision KSCALE, BOFFSET
      parameter (KSCALE=0.73D0, BOFFSET=0.09D0)
      double precision  dax(MAXPT), day(MAXPT), daz(MAXPT)      
cRE
	double precision sumda(MAXPT)
	double precision si, si2, rj
      double precision rgbmaxpsmax2, r2, dij1i, dij2i,dij3i, dij,
     &                 temp1, temp7, sj,sj2, datmp, rgbmax2, rgbmax1i,
     &                 rgbmax2i, dumbo, tmpsd,xij,yij,zij
      double precision  thi,multfact,expmkf,dielfac,qi2h,qid2h 
      double precision q(MAXPT),xi,yi,zi,ri,ri1i,daix,daiy,daiz


c YS FGB taylor coefficients follow from A to H;
c    1/3, 2/5, 3/7, 4/9, 5/11, 4/3, 12/5, 24/7, 40/9, 60/11
      double precision TA, TB, TC, TD, TDD, TE, TF, TG, TH, THH
      parameter (TA =0.33333333333333333333)
      parameter (TB =0.4D0)
      parameter (TC =0.42857142857142857143)
      parameter (TD =0.44444444444444444444)
      parameter (TDD=0.45454545454545454545)
      parameter (TE =1.33333333333333333333)
      parameter (TF =2.4D0)
      parameter (TG =3.42857142857142857143)
      parameter (TH =4.44444444444444444444)
      parameter (THH=5.45454545454545454545)
      
c Begin
      rgbmax2 = rgbmax*rgbmax
      rgbmax1i = 1.0 / rgbmax
      rgbmax2i = rgbmax1i*rgbmax1i
      rgbmaxpsmax2 = ( rgbmax + Fsmax ) * ( rgbmax + Fsmax )
    
      multfact = 18.220867158288597763891252209377D0
      do i = 1,natom
         q(i) = qmoil(i)*multfact
      end do 
c ---------------------------------------------------------------------
c Compute the "diagnal" energy that is a function of only the
c effective radius Ri. Also compute the contribution of 
c the diagnal energy term to the sum by which the derivative of
c Ri will be multiplied
c --------------------------------------------------------------------- 


      egb_self_ene = 0.0
	jbeg1 = 1
	jbeg3 = 1
      do i =1, natom 
	 expmkf = DEXP(-KSCALE * kappa * reff(i) )/ diel_ext
	 dielfac = 1.0D0 - expmkf
	 qi2h = 0.5 * q(i)*q(i)
	 if (lesid(i).ne.0) qi2h = qi2h/ptwei(i)
	 qid2h = qi2h * dielfac
         egb_self_ene = egb_self_ene - qid2h / reff(i)
	 sumdeijda(i) = sumdeijda(i) + qid2h - KSCALE*kappa*qi2h*
     &               expmkf*reff(i)
      end do
c      write(*,*) "Self-energy: ", egb_self_ene

      epol = epol + egb_self_ene

c RE: initializig dax/y/z(i) (not used in the previous version) to allow
c symmetric treatment of i and j.
c
c Comment: there are a lot of redundant calculations done here and
c in egb_f. Merging the two routines will make the code even messier, but 
c probably more efficient.
c
	do i = 1, natom
	 dax(i) = 0.d0
	 day(i) = 0.d0
	 daz(i) = 0.d0
	end do
	 


c --------------------------------------------------------------------- 
c Compute the derivatives of the effective radius Ri of atom i
c w.r.t. the cartersian coordinates of each atom j (datmp),
c --------------------------------------------------------------------- 
      do i = 1, natom

	 xi = coor(1, i)
	 yi = coor(2, i)
	 zi = coor(3, i)
	 ri = rborn(i) - BOFFSET 
	 si = fs(i)*ri
	 si2 = si * si
	 ri1i = 1.0D0/ri

         sumda(i) = sumdeijda(i)
         if (gbobcbool) then
             thi = DTANH( ( gbalpha - gbbeta * psi(i) + 
     &                      gbgamma * psi(i) * psi(i) ) * psi(i) )
             sumda(i) = sumda(i) * ( (gbalpha - 2.0 * gbbeta * psi(i) + 
     &               3.0 * gbgamma * psi(i) * psi(i) ) * ( 1.0 - 
     &               thi * thi ) * ri/rborn(i) )
         end if 

	 jend1 = point1(i)

c RE: the loop below is for list 1. 
	 do k = jbeg1,jend1
	     j = list1(k)
	     if (lesid(i).ne.0 .and. lesid(j).ne.0
     &		 .and. (lesid(i).eq.lesid(j) .and.
     &		  cplbl(i).ne.cplbl(j)) ) go to 1001 
	     xij = xi - coor(1, j)
	     yij = yi - coor(2, j)
	     zij = zi - coor(3, j)
	     r2=xij*xij + yij*yij+zij*zij
	     dij1i = 1.0 /DSQRT(r2)
	     dij2i = dij1i * dij1i
	     dij= r2*dij1i
	     rj = rborn(j) - BOFFSET
	     sj = fs(j) * rj
	     sj2 = sj * sj

c
c The following are the numerator of the first derivatives of the  
c effective redius Ri w.r.t. the interatomic distance Dij
c
c RE: The loooong if below is for the j particles. Similarly loooong
c if follows this sequence for the i particle
c
	     if (dij .gt. rgbmax + sj) go to 1001
	     if (dij .gt. rgbmax -sj ) then
		 temp1 = 1.0/(dij - sj)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + sj2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*sj) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  ri+sj) then 
		 temp1 = 1.0/(r2 - sj2)
		 datmp = temp1 * sj * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - sj) / (dij + sj))
	     else if ( dij .gt. DABS(ri-sj) ) then
		 temp1 = 1.0 /( dij + sj)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - ri * ri + sj2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(ri * temp1))
	     else if ( ri .lt. sj) then
		 temp1 = 1.0/(r2-sj2)
		 datmp = -0.5 * (sj * dij2i *temp1 - 2.0 * sj *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((sj - dij ) / (sj + dij)))
	     else
		 datmp = 0
	     end if

	     dax(i) = dax(i) + xij * datmp                  
	     day(i) = day(i) + yij * datmp                  
	     daz(i) = daz(i) + zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,j) = dpot(1,j) + xij * datmp/ptwei(i)
             	 dpot(2,j) = dpot(2,j) + yij * datmp/ptwei(i)
             	 dpot(3,j) = dpot(3,j) + zij * datmp/ptwei(i)
		else
             	 dpot(1,j) = dpot(1,j) + xij * datmp
             	 dpot(2,j) = dpot(2,j) + yij * datmp
             	 dpot(3,j) = dpot(3,j) + zij * datmp
		end if

1001         continue 
cRE: for symmetry compute the derivative of j born radius
c


	     if (dij .gt. rgbmax + si) go to 1002
	     if (dij .gt.  rgbmax -si ) then
		 temp1 = 1.0/(dij - si)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + si2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*si) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  rj+si) then 
		 temp1 = 1.0/(r2 - si2)
		 datmp = temp1 * si * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - si) / (dij + si))
	     else if ( dij .gt. DABS(rj-si) ) then
		 temp1 = 1.0 /( dij + si)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - rj * rj + si2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(rj * temp1))
	     else if ( rj .lt. si) then
		 temp1 = 1.0/(r2-si2)
		 datmp = -0.5 * (si * dij2i *temp1 - 2.0 * si *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((si - dij ) / (si + dij)))
	     else
		 datmp = 0
	     end if

	     dax(j) = dax(j) - xij * datmp                  
	     day(j) = day(j) - yij * datmp                  
	     daz(j) = daz(j) - zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,i) = dpot(1,i) - xij * datmp/ptwei(i)
             	 dpot(2,i) = dpot(2,i) - yij * datmp/ptwei(i)
             	 dpot(3,i) = dpot(3,i) - zij * datmp/ptwei(i)
		else
             	 dpot(1,i) = dpot(1,i) - xij * datmp
             	 dpot(2,i) = dpot(2,i) - yij * datmp
             	 dpot(3,i) = dpot(3,i) - zij * datmp
		end if

1002         continue 

c end do j
	 end do         
	 jbeg1 = jend1 + 1

cRE: get ready for list3

	 
	 jend3 = point3(i)
c RE: the loop below is for list 1. 
	 do k = jbeg3,jend3
	     j = list3(k)
	     if (lesid(i).ne.0 .and. lesid(j).ne.0
     &		 .and. (lesid(i).eq.lesid(j) .and.
     &		  cplbl(i).ne.cplbl(j)) ) go to 1001 
	     xij = xi - coor(1, j)
	     yij = yi - coor(2, j)
	     zij = zi - coor(3, j)
	     r2=xij*xij + yij*yij+zij*zij
	     dij1i = 1.0 /DSQRT(r2)
	     dij2i = dij1i * dij1i
	     dij= r2*dij1i
	     rj = rborn(j) - BOFFSET
	     sj = fs(j) * rj
	     sj2 = sj * sj

c
c The following are the numerator of the first derivatives of the  
c effective redius Ri w.r.t. the interatomic distance Dij
c
c RE: The loooong if below is for the j particles. Similarly loooong
c if follows this sequence for the i particle
c
	     if (dij .gt. rgbmax + sj) go to 1003
	     if (dij .gt. rgbmax -sj ) then
		 temp1 = 1.0/(dij - sj)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + sj2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*sj) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  ri+sj) then 
		 temp1 = 1.0/(r2 - sj2)
		 datmp = temp1 * sj * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - sj) / (dij + sj))
	     else if ( dij .gt. DABS(ri-sj) ) then
		 temp1 = 1.0 /( dij + sj)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - ri * ri + sj2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(ri * temp1))
	     else if ( ri .lt. sj) then
		 temp1 = 1.0/(r2-sj2)
		 datmp = -0.5 * (sj * dij2i *temp1 - 2.0 * sj *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((sj - dij ) / (sj + dij)))
	     else
		 datmp = 0
	     end if

	     dax(i) = dax(i) + xij * datmp                  
	     day(i) = day(i) + yij * datmp                  
	     daz(i) = daz(i) + zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,j) = dpot(1,j) + xij * datmp/ptwei(i)
             	 dpot(2,j) = dpot(2,j) + yij * datmp/ptwei(i)
             	 dpot(3,j) = dpot(3,j) + zij * datmp/ptwei(i)
		else
             	 dpot(1,j) = dpot(1,j) + xij * datmp
             	 dpot(2,j) = dpot(2,j) + yij * datmp
             	 dpot(3,j) = dpot(3,j) + zij * datmp
		end if

1003         continue 
cRE: for symmetry compute the derivative of j born radius
c


	     if (dij .gt. rgbmax + si) go to 1004
	     if (dij .gt.  rgbmax -si ) then
		 temp1 = 1.0/(dij - si)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + si2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*si) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  rj+si) then 
		 temp1 = 1.0/(r2 - si2)
		 datmp = temp1 * si * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - si) / (dij + si))
	     else if ( dij .gt. DABS(rj-si) ) then
		 temp1 = 1.0 /( dij + si)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - rj * rj + si2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(rj * temp1))
	     else if ( rj .lt. si) then
		 temp1 = 1.0/(r2-si2)
		 datmp = -0.5 * (si * dij2i *temp1 - 2.0 * si *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((si - dij ) / (si + dij)))
	     else
		 datmp = 0
	     end if

	     dax(j) = dax(j) - xij * datmp                  
	     day(j) = day(j) - yij * datmp                  
	     daz(j) = daz(j) - zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,i) = dpot(1,i) - xij * datmp/ptwei(i)
             	 dpot(2,i) = dpot(2,i) - yij * datmp/ptwei(i)
             	 dpot(3,i) = dpot(3,i) - zij * datmp/ptwei(i)
		else
             	 dpot(1,i) = dpot(1,i) - xij * datmp
             	 dpot(2,i) = dpot(2,i) - yij * datmp
             	 dpot(3,i) = dpot(3,i) - zij * datmp
		end if

1004         continue 

c end do j for list3
	 end do         
	 jbeg3 = jend3 + 1

c RE: now the exclusion list (which is part of the game now...)
         jbegex = exc1(i-1)+1
         jendex = exc1(i)

         do k = jbegex,jendex
	  j = exc2(k)
	     if (lesid(i).ne.0 .and. lesid(j).ne.0
     &		 .and. (lesid(i).eq.lesid(j) .and.
     &		  cplbl(i).ne.cplbl(j)) ) go to 1001 
	     xij = xi - coor(1, j)
	     yij = yi - coor(2, j)
	     zij = zi - coor(3, j)
	     r2=xij*xij + yij*yij+zij*zij
	     dij1i = 1.0 /DSQRT(r2)
	     dij2i = dij1i * dij1i
	     dij= r2*dij1i
	     rj = rborn(j) - BOFFSET
	     sj = fs(j) * rj
	     sj2 = sj * sj


c The following are the numerator of the first derivatives of the  
c effective redius Ri w.r.t. the interatomic distance Dij
c
c RE: The loooong if below is for the j particles. Similarly loooong
c if follows this sequence for the i particle
c
	     if (dij .gt. rgbmax + sj) go to 1005
	     if (dij .gt. rgbmax -sj ) then
		 temp1 = 1.0/(dij - sj)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + sj2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*sj) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  ri+sj) then 
		 temp1 = 1.0/(r2 - sj2)
		 datmp = temp1 * sj * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - sj) / (dij + sj))
	     else if ( dij .gt. DABS(ri-sj) ) then
		 temp1 = 1.0 /( dij + sj)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - ri * ri + sj2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(ri * temp1))
	     else if ( ri .lt. sj) then
		 temp1 = 1.0/(r2-sj2)
		 datmp = -0.5 * (sj * dij2i *temp1 - 2.0 * sj *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((sj - dij ) / (sj + dij)))
	     else
		 datmp = 0
	     end if

	     dax(i) = dax(i) + xij * datmp                  
	     day(i) = day(i) + yij * datmp                  
	     daz(i) = daz(i) + zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,j) = dpot(1,j) + xij * datmp/ptwei(i)
             	 dpot(2,j) = dpot(2,j) + yij * datmp/ptwei(i)
             	 dpot(3,j) = dpot(3,j) + zij * datmp/ptwei(i)
		else
             	 dpot(1,j) = dpot(1,j) + xij * datmp
             	 dpot(2,j) = dpot(2,j) + yij * datmp
             	 dpot(3,j) = dpot(3,j) + zij * datmp
		end if

1005         continue 
cRE: for symmetry compute the derivative of j born radius
c


	     if (dij .gt. rgbmax + si) go to 1006
	     if (dij .gt.  rgbmax -si ) then
		 temp1 = 1.0/(dij - si)
                 dij3i  = dij1i * dij2i
		 datmp = 0.125 * dij3i * (( r2 + si2) *
     &                       (temp1 * temp1 - rgbmax2i) -
     &                       2.0 * DLOG(rgbmax*temp1))

	     else if (dij .gt. 4.0*si) then
		 tmpsd = sj2 * dij2i
		 dumbo = TE + tmpsd * (TF+tmpsd*(TG+tmpsd*
     &                       (TH+tmpsd*THH)))
		 datmp = tmpsd*sj*dij2i*dij2i*dumbo
	     else if (dij .gt.  rj+si) then 
		 temp1 = 1.0/(r2 - si2)
		 datmp = temp1 * si * (-0.5 * dij2i + temp1) +
     &                       0.25 * dij1i * dij2i *
     &                       DLOG((dij - si) / (dij + si))
	     else if ( dij .gt. DABS(rj-si) ) then
		 temp1 = 1.0 /( dij + si)
		 dij3i = dij1i * dij1i * dij1i
		 datmp = -0.25 * (-0.5 * ( r2 - rj * rj + si2) *
     &                       dij3i * ri1i * ri1i + dij1i * temp1 *
     &                       (temp1 - dij1i) - dij3i * 
     &                       DLOG(rj * temp1))
	     else if ( rj .lt. si) then
		 temp1 = 1.0/(r2-si2)
		 datmp = -0.5 * (si * dij2i *temp1 - 2.0 * si *
     &                       temp1 * temp1 - 0.5 * dij2i * dij1i *
     &                       DLOG((si - dij ) / (si + dij)))
	     else
		 datmp = 0
	     end if

	     dax(j) = dax(j) - xij * datmp                  
	     day(j) = day(j) - yij * datmp                  
	     daz(j) = daz(j) - zij * datmp                   

             datmp = datmp * sumda(i)
	      if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
             	 dpot(1,i) = dpot(1,i) - xij * datmp/ptwei(i)
             	 dpot(2,i) = dpot(2,i) - yij * datmp/ptwei(i)
             	 dpot(3,i) = dpot(3,i) - zij * datmp/ptwei(i)
		else
             	 dpot(1,i) = dpot(1,i) - xij * datmp
             	 dpot(2,i) = dpot(2,i) - yij * datmp
             	 dpot(3,i) = dpot(3,i) - zij * datmp
		end if

1006         continue 

c end do j
	 end do         
cRE: now add all the accumulated forces in dax to dpot(i)
c end do i
	end do

	do i = 1,natom

	     if (lesid(i).ne.0 .and. (lesid(i).eq.lesid(j))) then
         	dpot(1,i) = dpot(1,i) - sumda(i) * dax(i)/ptwei(i)
         	dpot(2,i) = dpot(2,i) - sumda(i) * day(i)/ptwei(i)
         	dpot(3,i) = dpot(3,i) - sumda(i) * daz(i)/ptwei(i)
	     else
         	dpot(1,i) = dpot(1,i) - sumda(i) * dax(i)
         	dpot(2,i) = dpot(2,i) - sumda(i) * day(i)
         	dpot(3,i) = dpot(3,i) - sumda(i) * daz(i)
	     end if

c end do i    
      end do

      return 
      END

C --------------------------------------------------
      SUBROUTINE egb_f(natom,qmoil,diel_ext,
     $     ifsndbool,moilvector,moilvectorout,
     $     moildebug)


      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/MASSFAC.BLOCK'
C --------------------------------------------------



C - - - arg types - - -

      integer  natom,ifsndbool,moildebug
      double precision qmoil(*),moilvector(*),moilvectorout(*),
     $     diel_ext,epol

C     - - - local declarations - - -

      integer  i,j,k,size_da,ic,iaci,skip,iexcl,jexcl,jexcl_last,nexcl,
     $         atompick, natom3,iloop,jloop,jbeg1,jend1,jbeg3,jend3,
     &         jbegex, jendex,l ,m,n,kl,ll,il,ik,ii,jj,runstep

      double precision 
     &       e,dielfac,qi,qiqj,expmkf,q(MAXPT),
     $       xi,yi,zi,xij,yij,zij,dij1i,dij3i,
     $       temp,temp1,qi2h,qid2h,datmp,daix,daiy,daiz,theta,ri1i,
     &       dij2i,dij, eel,f6,f12,rinv,r2inv,r6inv,r2,ri,rj,
     &       sj,sj2,uij,multfact, Atermsave,ftermsnd,fterm1st,term,
     &       mass,expmkfi,expmkfj,dielfaci, dielfacj,qj2h,qjd2h,
     &       ftermsndi,ftermsndj,fterm,ftermsum,betaderivk,betaderivl,
     &       psiderivk,psiderivl,bit3,part2

C RE
	double precision sum(maxpt)
	double precision si2,xj,yj,zj,uji
C YS variables for Onufriev's GB model 
      double precision BOFFSET, KSCALE
      parameter ( BOFFSET = 0.09D0, KSCALE = 0.73D0 )

C YS LCPO stuff follows
      double precision si 
      double precision dumbo, tmpsd
      double precision rgbmax2,rgbmax1i,rgbmax2i,tmpcs, rgbmaxpsmax2 
      integer count, counta
c RE bug fixd for lES add wqi variable
	double precision wqi
C@	
	double precision dpot_tmp

C YS FGB taylor coefficients follow from A to H;
C    1/3, 2/5, 3/7, 4/9, 5/11, 4/3, 12/5, 24/7, 40/9, 60/11
      double precision TA, TB, TC, TD, TDD, TE, TF, TG, TH, THH
      parameter (TA =0.33333333333333333333)
      parameter (TB =0.4D0)
      parameter (TC =0.42857142857142857143)
      parameter (TD =0.44444444444444444444)
      parameter (TDD=0.45454545454545454545)
      parameter (TE =1.33333333333333333333)
      parameter (TF =2.4D0)
      parameter (TG =3.42857142857142857143)
      parameter (TH =4.44444444444444444444)
      parameter (THH=5.45454545454545454545)


c Begin execution

c YS   based on NAB code eff.c
      rgbmax2 = rgbmax*rgbmax
      rgbmax1i = 1.0/rgbmax
      rgbmax2i = rgbmax1i*rgbmax1i
      rgbmaxpsmax2 = ( rgbmax + Fsmax ) * ( rgbmax + Fsmax )


      natom3 = natom*3
      multfact = 18.220867158288597763891252209377D0
      do i = 1,natom
         q(i) = qmoil(i)*multfact
	 sum(i) = 0.d0
      end do 

c YS: haven't used 'runstep', what's the purpose for ?

      size_da = 3*(natom)
      if (.not.(gbsu.eq.0)) then 
         runstep = MOD(stepno,(gbsu))
      else 
         if (stepno.gt.0) then 
            runstep = 0
         else 
            runstep = 1
         end if 
      end if 

c -----------------------------------------------------------
c Algorithm Step 1
c Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275
c -----------------------------------------------------------

c YS  Get the "effective" Born radii via the approximate pairwise
c     method Use Eqs 9-11 of Hawkins, Cramer, Truhlar; J. Phys. 
c     Chem. 100:19824(1996). Also added a modified method based on
c     Onufriev, Bashford, Case, Proteins, 55:383(2004). eqs 6
c
      count = 1      
      counta = 1
      ineighbor_ptrs(1) = 1
      overlappair_ptrs(1) = 1
	jbeg1 = 1
	jbeg3 = 1
      do i = 1,natom
	  jend1 = point1(i)
          xi = coor(1, i)
          yi = coor(2, i)
          zi = coor(3, i)
          ri = rborn( i ) - BOFFSET 
          ri1i = 1.0D0/ri
	  si = fs(i)*ri
	  si2 = si*si

c YS: NAB code select j from the pair list, may be considered 
c in the furture
          do k = jbeg1,jend1
	      j = list1(k)
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then
c@
c		 if (lesid(i).ne.0) write(*,*) i,j,lesid(i),lesid(j)
		  xj  = coor(1, j)
		  yj  = coor(2, j)
		  zj  = coor(3, j)
		  rj  = rborn(j) - BOFFSET
                  xij = xi - xj
                  yij = yi - yj
                  zij = zi - zj
                  r2 = xij * xij + yij * yij + zij * zij
                  if ( r2 .gt. rgbmaxpsmax2 ) go to 100
                  dij1i = 1.0 / DSQRT(r2)
                  dij = r2 * dij1i
                  sj = fs(j) * rj
                  sj2 = sj * sj
 
c YS: from NAB code: 
c     followings are from the Appendix of Schaefer and Froemmel
c     JMB 216:1045-1066, 1990; Taylor series expansion for d >> s
c     is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
c     Andreas Svrcek-Seiler and Alexy Onufriev  

c RE
c the text below is now running twice, once for the derivative of the Born
c radius of i once for j
c
                  if ( dij .gt. rgbmax+sj ) go to 100
                  if ( dij .gt. rgbmax-sj ) then 
                     uij = 1.0 / ( dij - sj )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uij+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - sj2)+
     &                      2.0 * DLOG((dij - sj) * rgbmax1i))
                     sum(i) = sum(i) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = sj2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(i) = sum(i) - sj * tmpsd * dij2i * dumbo                
                  else if (dij .gt. ri + sj ) then 
		     temp =  0.5 *( sj/(r2-sj2) +
     $                    0.5 * dij1i * DLOG((dij-sj)/(dij+sj)))
                     sum(i) = sum(i) - temp
                  else if (dij .gt. DABS(ri-sj) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + ri*ri -sj2)
                     uij = 1.0 / (dij+sj)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uij + dij1i*DLOG(ri*uij))
                     sum(i) = sum(i) - temp
                  else if ( ri .lt. sj ) then
		     temp = 0.5 * (sj / (r2 - sj2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((sj - dij)/(sj + dij)))
                     sum(i) = sum(i) - temp

                  end if 
                  if ( dij .gt. rgbmax+si ) go to 100
                  if ( dij .gt. rgbmax-si ) then 
                     uji = 1.0 / ( dij - si )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uji+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - si2)+
     &                      2.0 * DLOG((dij - si) * rgbmax1i))
                     sum(j) = sum(j) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = si2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(j) = sum(j) - si * tmpsd * dij2i * dumbo                
                  else if (dij .gt. rj + si ) then 
		     temp =  0.5 *( si/(r2-si2) +
     $                    0.5 * dij1i * DLOG((dij-si)/(dij+si)))
                     sum(j) = sum(j) - temp
                  else if (dij .gt. DABS(rj-si) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + rj*rj -si2)
                     uji = 1.0 / (dij+si)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uji + dij1i*DLOG(rj*uji))
                     sum(j) = sum(j) - temp
                  else if ( rj .lt. si ) then
		     temp = 0.5 * (si / (r2 - si2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((si - dij)/(si + dij)))
                     sum(j) = sum(j) - temp

                  end if 
              end if

              if( gbnpbool ) then 
                  if( (P0(ptsaid(i)) + P0(ptsaid(j))) .gt. dij ) then
                      if ( (P0(ptsaid(i)) .gt. 2.5 ) .and. 
     &                     (P0(ptsaid(j)) .gt. 2.5 ) ) then
                           ineighbor(count) = j
                           count = count + 1
		           if ( j .gt. i ) then
			       overlappair(counta) = j
			       counta = counta + 1
                           end if
                      end if
                  end if
              end if 

100       continue  
          end do 
	  jbeg1 = jend1 + 1

	  jend3 = point3(i)

          do k = jbeg3,jend3
	      j = list3(k)
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then
c@
c		 if (lesid(i).ne.0) write(*,*) i,j,lesid(i),lesid(j)
		  xj  = coor(1, j)
		  yj  = coor(2, j)
		  zj  = coor(3, j)
		  rj  = rborn(j) - BOFFSET
                  xij = xi - xj
                  yij = yi - yj
                  zij = zi - zj
                  r2 = xij * xij + yij * yij + zij * zij
                  if ( r2 .gt. rgbmaxpsmax2 ) go to 101
                  dij1i = 1.0 / DSQRT(r2)
                  dij = r2 * dij1i
                  sj = fs(j) * rj
                  sj2 = sj * sj
 
c YS: from NAB code: 
c     followings are from the Appendix of Schaefer and Froemmel
c     JMB 216:1045-1066, 1990; Taylor series expansion for d >> s
c     is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
c     Andreas Svrcek-Seiler and Alexy Onufriev  

c RE
c the text below is now running twice, once for the derivative of the Born
c radius of i once for j
c
                  if ( dij .gt. rgbmax+sj ) go to 101
                  if ( dij .gt. rgbmax-sj ) then 
                     uij = 1.0 / ( dij - sj )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uij+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - sj2)+
     &                      2.0 * DLOG((dij - sj) * rgbmax1i))
                     sum(i) = sum(i) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = sj2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(i) = sum(i) - sj * tmpsd * dij2i * dumbo                
                  else if (dij .gt. ri + sj ) then 
		     temp =  0.5 *( sj/(r2-sj2) +
     $                    0.5 * dij1i * DLOG((dij-sj)/(dij+sj)))
                     sum(i) = sum(i) - temp
                  else if (dij .gt. DABS(ri-sj) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + ri*ri -sj2)
                     uij = 1.0 / (dij+sj)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uij + dij1i*DLOG(ri*uij))
                     sum(i) = sum(i) - temp
                  else if ( ri .lt. sj ) then
		     temp = 0.5 * (sj / (r2 - sj2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((sj - dij)/(sj + dij)))
                     sum(i) = sum(i) - temp

                  end if 
                  if ( dij .gt. rgbmax+si ) go to 100
                  if ( dij .gt. rgbmax-si ) then 
                     uji = 1.0 / ( dij - si )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uji+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - si2)+
     &                      2.0 * DLOG((dij - si) * rgbmax1i))
                     sum(j) = sum(j) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = si2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(j) = sum(j) - si * tmpsd * dij2i * dumbo                
                  else if (dij .gt. rj + si ) then 
		     temp =  0.5 *( si/(r2-si2) +
     $                    0.5 * dij1i * DLOG((dij-si)/(dij+si)))
                     sum(j) = sum(j) - temp
                  else if (dij .gt. DABS(rj-si) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + rj*rj -si2)
                     uji = 1.0 / (dij+si)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uji + dij1i*DLOG(rj*uji))
                     sum(j) = sum(j) - temp
                  else if ( rj .lt. si ) then
		     temp = 0.5 * (si / (r2 - si2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((si - dij)/(si + dij)))
                     sum(j) = sum(j) - temp

                  end if 
              end if

              if( gbnpbool ) then 
                  if( (P0(ptsaid(i)) + P0(ptsaid(j))) .gt. dij ) then
                      if ( (P0(ptsaid(i)) .gt. 2.5 ) .and. 
     &                     (P0(ptsaid(j)) .gt. 2.5 ) ) then
                           ineighbor(count) = j
                           count = count + 1
		           if ( j .gt. i ) then
			       overlappair(counta) = j
			       counta = counta + 1
                           end if
                      end if
                  end if
              end if 

101       continue  
          end do 
	  jbeg3 = jend3 + 1

C RE: For compatibility with old code the usually ecluded interactions are 
C also included in the GB calculation, but I feel uncomfortable about it
C
         jbegex = exc1(i-1)+1
         jendex = exc1(i)

         do k = jbegex,jendex
            j = exc2(k)
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then
c@
c		 if (lesid(i).ne.0) write(*,*) i,j,lesid(i),lesid(j)
		  xj  = coor(1, j)
		  yj  = coor(2, j)
		  zj  = coor(3, j)
		  rj  = rborn(j) - BOFFSET
                  xij = xi - xj
                  yij = yi - yj
                  zij = zi - zj
                  r2 = xij * xij + yij * yij + zij * zij
                  if ( r2 .gt. rgbmaxpsmax2 ) go to 100
                  dij1i = 1.0 / DSQRT(r2)
                  dij = r2 * dij1i
                  sj = fs(j) * rj
                  sj2 = sj * sj
 
c YS: from NAB code: 
c     followings are from the Appendix of Schaefer and Froemmel
c     JMB 216:1045-1066, 1990; Taylor series expansion for d >> s
c     is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
c     Andreas Svrcek-Seiler and Alexy Onufriev  

c RE
c the text below is now running twice, once for the derivative of the Born
c radius of i once for j
c
                  if ( dij .gt. rgbmax+sj ) go to 102
                  if ( dij .gt. rgbmax-sj ) then 
                     uij = 1.0 / ( dij - sj )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uij+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - sj2)+
     &                      2.0 * DLOG((dij - sj) * rgbmax1i))
                     sum(i) = sum(i) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = sj2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(i) = sum(i) - sj * tmpsd * dij2i * dumbo                
                  else if (dij .gt. ri + sj ) then 
		     temp =  0.5 *( sj/(r2-sj2) +
     $                    0.5 * dij1i * DLOG((dij-sj)/(dij+sj)))
                     sum(i) = sum(i) - temp
                  else if (dij .gt. DABS(ri-sj) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + ri*ri -sj2)
                     uij = 1.0 / (dij+sj)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uij + dij1i*DLOG(ri*uij))
                     sum(i) = sum(i) - temp
                  else if ( ri .lt. sj ) then
		     temp = 0.5 * (sj / (r2 - sj2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((sj - dij)/(sj + dij)))
                     sum(i) = sum(i) - temp

                  end if 
                  if ( dij .gt. rgbmax+si ) go to 102
                  if ( dij .gt. rgbmax-si ) then 
                     uji = 1.0 / ( dij - si )
		     temp = 0.125 * dij1i * (1.0 +
     &                      2.0 * dij * uji+ rgbmax2i * (r2 -
     &                      4.0 * rgbmax * dij - si2)+
     &                      2.0 * DLOG((dij - si) * rgbmax1i))
                     sum(j) = sum(j) - temp
                  else if ( dij .gt. 4.0*sj) then 
                     dij2i = dij1i * dij1i
                     tmpsd = si2*dij2i
                     dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * 
     &                     (TD + tmpsd * TDD)))                    
                     sum(j) = sum(j) - si * tmpsd * dij2i * dumbo                
                  else if (dij .gt. rj + si ) then 
		     temp =  0.5 *( si/(r2-si2) +
     $                    0.5 * dij1i * DLOG((dij-si)/(dij+si)))
                     sum(j) = sum(j) - temp
                  else if (dij .gt. DABS(rj-si) ) then 
                     theta = 0.5 * ri1i*dij1i*(r2 + rj*rj -si2)
                     uji = 1.0 / (dij+si)
		     temp = 0.25 * ( ri1i*(2.0D0-theta) -
     $                    uji + dij1i*DLOG(rj*uji))
                     sum(j) = sum(j) - temp
                  else if ( rj .lt. si ) then
		     temp = 0.5 * (si / (r2 - si2) + 2.0 * ri1i +  
     $                    0.5 * dij1i * DLOG((si - dij)/(si + dij)))
                     sum(j) = sum(j) - temp

                  end if 
              end if

              if( gbnpbool ) then 
                  if( (P0(ptsaid(i)) + P0(ptsaid(j))) .gt. dij ) then
                      if ( (P0(ptsaid(i)) .gt. 2.5 ) .and. 
     &                     (P0(ptsaid(j)) .gt. 2.5 ) ) then
                           ineighbor(count) = j
                           count = count + 1
		           if ( j .gt. i ) then
			       overlappair(counta) = j
			       counta = counta + 1
                           end if
                      end if
                  end if
              end if 

102       continue  
	 end do
	end do
C ** end of double i & j loopp. Below i loop only

	do i=1,natom
          if(gbobcbool) then 
C YS  Onufriev effective radii
C     note sum(i) is negative here
               psi(i) = -ri * sum(i)
               reff(i) = 1.0 / (ri1i-DTANH(( gbalpha - gbbeta*psi(i)+
     &                  gbgamma*psi(i)*psi(i))*psi(i))/rborn(i))

          else 
C YS "standard" HCT effective radii
               reff(i) = 1.0 /(ri1i+ sum(i))
	       if(reff(i) .lt. 0.0 ) reff(i) = 30.0
          end if 

          if ( gbnpbool ) then 
                ineighbor(count) = 0
                count = count + 1
		ineighbor_ptrs(i+1) = count
		overlappair_ptrs(i+1) = counta
          end if
c@
c          write(*,*) "BORN RADII> i=",i, " reff=", reff(i)
      end do 

C** end of Algorithm Step 1

c  skip polarized term
c      go to  1776
    
c -----------------------------------------------------------
c Algorithm Step 2
c Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275
c -----------------------------------------------------------

c start calculating the polarized GB term (off-diagonal, 
c interaction energy ) 

      do i = 1, natom
          sumdeijda(i) = 0.0d0
      end do


      atompick = 0
      epol = 0
      iexcl = 0
      jbeg1 = 1
      jbeg3 = 1

C RE: Why not add this directly to step 1?? For the moment I left as is...
C
      do i = 1,natom-1

c atom j is selected from list 

         qi = q(i)
         ri = reff(i)
         jend1 = point1(i)
         do k = jbeg1,jend1
            j = list1(k)
	    xij = coor(1,i) - coor(1,j)
	    yij = coor(2,i) - coor(2,j)
	    zij = coor(3,i) - coor(3,j)
	    r2 = xij*xij + yij*yij + zij*zij
            if (.not.(q(i) .eq.0D0) .AND. .not.( q(j) .eq. 0D0)
     &		.AND. .not.(r2.gt.cutele2)) then 
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then

c@
c		if (lesid(i).ne.0) then
c@	 write(*,*)' in list 1: i,j,lesid(i),lesid(j) ',i,j,lesid(i),lesid(j)
c		end if
			if (lesid(i).eq.lesid(j)) then
				wqi = qi/ptwei(i)
			else
				wqi = qi
			end if
		if (ifsndbool.eq.1) then
                  call egb_calc_pair_f(coor(1,i), coor(1,j),
     $              wqi, q(j), ri, reff(j), epol,
     $              dpot(1,i),dpot(1,j),
     $              moilvectorout(3*(i-1)+1),
     $              moilvectorout(3*(j-1)+1), 
     $              moilvector(3*(i-1)+1),
     $              moilvector(3*(j-1)+1),
     $              massfac(i), massfac(j),
     $              diel_ext, ifsndbool,i,j)
		else 
C@ RE
c		 write(*,*) ' before calc_pair  i j wqi q(j) ', i,j,wqi,q(j)
c		 write(*,*) ' ri reff(j) ',ri,reff(j)
		  call egb_calc_pair_f(coor(1,i),coor(1,j),
     $              wqi, q(j), ri, reff(j), epol,
     $              dpot(1,i),dpot(1,j),
     $              0.d0,
     $              0.d0, 
     $              0.d0,
     $              0.d0,
     $              massfac(i), massfac(j),
     $              diel_ext, ifsndbool,i,j)
c		 write(*,*) ' epol = ',epol
		end if
               end if
            end if 
         end do 

         jbeg1 = jend1+1
         jend3 = point3(i)
         do k = jbeg3,jend3
            j = list3(k)
	    xij = coor(1,i) - coor(1,j)
	    yij = coor(2,i) - coor(2,j)
	    zij = coor(3,i) - coor(3,j)
	    r2 = xij*xij + yij*yij + zij*zij
            if (.not.(q(i) .eq.0D0) .AND. .not.( q(j) .eq. 0D0)
     &		.AND. .not.(r2.gt.cutele2)) then 
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then
			if (lesid(i).eq.lesid(j)) then
				wqi = qi/ptwei(i)
			else
				wqi = qi
			end if
                   CALL egb_calc_pair_f(coor(1,i), coor(1,j),
     $              wqi, q(j), ri, reff(j), epol,
     $              dpot(1,i),dpot(1,j), 
     $              moilvectorout(3*(i-1)+1),
     $              moilvectorout(3*(j-1)+1),
     $              moilvector(3*(i-1)+1),
     $              moilvector(3*(j-1)+1),
     $              massfac(i), massfac(j),
     $              diel_ext, ifsndbool,i,j)
c@
		dpot_tmp = dpot(1,i)*dpot(1,i) + dpot(2,i)*dpot(2,i)
     $		+ dpot(3,i)*dpot(3,i)
		dpot_tmp = dsqrt(dpot_tmp)
		if (dpot_tmp.gt.50) then
		 write(*,*)' i, j dpot ' ,i,j,dpot_tmp
		end if
            end if 
	   end if
         end do 

         jbeg3 = jend3+1
         jbegex = exc1(i-1)+1
         jendex = exc1(i)

         do k = jbegex,jendex
            j = exc2(k)
            if (.not.(q(i) .eq. 0.0D0) .AND. .not.(q(j).eq.0.0D0)) then 
              if (lesid(i).eq.0 .or. lesid(j).eq.0 .or.
     &		(lesid(i).eq.lesid(j) .and. cplbl(i).eq.cplbl(j))) then
			if (lesid(i).eq.lesid(j)) then
				wqi = qi/ptwei(i)
			else
				wqi = qi
			end if
c@
c		if (lesid(i).ne.0) then
c	 write(*,*)' in x-list : i,j,lesid(i),lesid(j) ',i,j,lesid(i),lesid(j)
c		end if
c@
c		if (lesid(i).ne.0 .and. lesid(i).eq.lesid(j)) then
c		write(*,*)' i j qi qj wqi ',i,j,qi,q(j),wqi
c		end if
                   call egb_calc_pair_f(coor(1,i), coor(1,j),
     $                wqi, q(j), ri, reff(j), epol,
     $                dpot(1,i),dpot(1,j),
     $                moilvectorout(3*(i-1)+1),
     $                moilvectorout(3*(j-1)+1), 
     $                moilvector(3*(i-1)+1), 
     $                moilvector(3*(j-1)+1), 
     $                massfac(i), massfac(j),
     $                diel_ext, ifsndbool,i,j)

c                  write (6,66) 'loop 3 ',stepno,i-1,j-1,epol
c                  write (6,*) 'dpoti ',dpot(1,i),dpot(2,i),dpot(3,i)
c                  write (6,*) 'dpotj ',dpot(1,j),dpot(2,j),dpot(3,j)
               end if
            end if  
         end do 
      end do 


c --------------------------------------------------------------------- 
c Algorithm Step 3
c Algorithm see:  Tsui and Case, Biopolymers, 2001, 56, 275
c --------------------------------------------------------------------- 

c Compute diagonal term (self-energy)
c and derivatives of effective radius w.r.t Cartesian coord

      call egb_deriv_rad(natom,qmoil, diel_ext, epol)

1776  continue

c ---------------------------------------------------------------------
c Compute non-polarized solvation energy (surface area based)
c --------------------------------------------------------------------- 

      if (gbnpbool) then
	  call egb_nonpol2(natom, e_gbnonpol)
          e_gbsa = epol + e_gbnonpol
      else 
          e_gbsa = epol
      end if

      return 
      END 






