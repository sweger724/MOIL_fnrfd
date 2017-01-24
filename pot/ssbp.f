C----------------------------------------------------------
C Spherical Solvent Boundary Potential (SSBP)
C
C Authors:  Dmitrii Beglovd and Benoit Roux  (1994)
C           Department of Chemistry, University of Montreal
C

C Default is include all energy terms for the 
C    solvent boundary potential

c     ----------------------------
      SUBROUTINE FACTO(NFACT)
c     ----------------------------
c Calculate FACT=NFACT!
c
      include 'COMMON/LENGTH.BLOCK' 
      include 'COMMON/DEBUG.BLOCK' 
      include 'COMMON/SSBP.BLOCK'
C      IMPLICIT NONE
c Output variable
C      double precision FACT(*)
c Input variable
      INTEGER NFACT
c Local variables
      INTEGER I,J
c
      FACT(1)=1.D0
      IF(NFACT.GT.0)THEN
      DO 10 I=1,NFACT
      J=I+1
      FACT(J)=FACT(I)*I
   10 CONTINUE
      ELSE
      ENDIF
      RETURN
      END
c
c     -------------------------------------------------------
      SUBROUTINE SSBP1
c     -------------------------------------------------------
c
c Total Spherical Solvent Boundary Potential = EKIRK + EHSR + EANGU + 
c

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK' 
     
      ENSSBP=0.d0
      ENKIR=0.d0
      ENCAV=0.d0
      ENHSR=0.d0
      ENANGU=0.d0
      ENHARM=0.d0
c
c Get the radius from the center of the sphere for all particles
      CALL DEFRAD

c Kirkwood's multipolar expansion for reaction field
      IF (QKIRK) THEN
          CALL SMEXP
      ENDIF

c Cavity potential of mean force calculated from RISM-HNC
      IF (QCAVI) THEN
          CALL SSBCAVP
      ENDIF

c Hard-Sphere-Restriction contribution = P*V+sigma*S 
      IF (QHSR) THEN
          CALL SSBHSRP
      ENDIF

c Angular correction potential for isotropic distribution at boundary
      IF (QANGU) THEN
          CALL SSBANGP
      ENDIF 
      
c Half harmonic restraint when fixed boundary simulation
      IF (QFIX.and.QHARM) THEN
          CALL HALF_HARM()
      ENDIF

      ENSSBP=ENKIR+ENCAV+ENHSR+ENANGU+ENHARM 
      RETURN 
      END 
c     -------------------------------- 
      SUBROUTINE DEFRAD
c     --------------------------------

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONNECT.BLOCK'

      integer i,j,namel,level
      character*8 name

c Calculate radial distance (RADIST) from the center
c and RMAX by choice (QFIX/.NOT.QFIX)
      name='DEFRAD'
      namel=6
      do 10 i=1,npt
      rsxy(i)=coor(1,i)*coor(1,i)+coor(2,i)*coor(2,i)
      ratom2(i)=rsxy(i)+coor(3,i)*coor(3,i)
      radist(i)=dsqrt(ratom2(i))
      rrdist(i)=1.d0/radist(i)
 10   continue
      IPA=0
      NHARM = 0

      IF(.NOT.QFIX)THEN
         RMAXM=0.d0
      ELSE
         RMAXM=FXCAS
      ENDIF

      DO 11 J=1,NTSSBP

        I=LSTSSBP(J)
       IF(RADIST(I).GT.RMAXM)THEN

        IF(.NOT.QFIX)THEN
          RMAXM=RADIST(I)
          IPA=I
        ELSE
         if (QHARM) then
           NHARM=NHARM+1
           HARM(NHARM)=I
         else
           if (qkirk) then
           level = 1
           call alert(name,namel,'Atom is outside the fixed 
     1               rmax',17,level)
           endif
         endif
        ENDIF

       ENDIF

   11 CONTINUE

      RETURN
      END  

c     -------------------------------------------
      SUBROUTINE SMEXP
c     -------------------------------------------
c Reaction field (Kirkwood's) part

c      IMPLICIT NONE

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'


c     local variables
	integer istart,iend
	logical first
      INTEGER I,J,K,L,M,N,ILP,ICOD,NATOM
      REAL*8 CHTOT,DRM1,DDRM1,RMAX1,RRSXY,RRATOM2,TVARI,RRATOM3
Cdeb      REAL*8 DRM1,DDRM1,RMAX1,RSXY,RATOM2,TVARI,RATOM3
      REAL*8 COEF0,RMAXF,RMAXLU,RMAXS,COEF,COL,COLX,COLY,COLZ
      REAL*8 QLMS
     
      integer II,IIM1
      COMPLEX*16 COM0,COM
      double precision YRL(maxpt),XR2(maxpt),
     1                 YR2(maxpt),ZR2(maxpt)
      double precision COMR(nalm),COMI(nalm)
      double precision RCH(nalm),TVAR(nalm),DTVAR(nalm)
C-------nalm = maxpt * (nmult+1)
      double precision QLMSX(maxpt),QLMSY(maxpt),QLMSZ(maxpt)
      double precision YLMR1(maxpt),YLMR2(maxpt),
     1                 YLMR4(maxpt),YLMR3(maxpt)
      double precision YLMI1(maxpt),YLMI2(maxpt),
     1                 YLMI4(maxpt),YLMI3(maxpt)
      double precision RLR3(maxpt),ZYR3(maxpt),
     1                 XRL(maxpt),ZXR3(maxpt)
      double precision TVAR12(maxpt)
	data first/.true./
	save chtot,first,coef0

      natom = npt
	if (prll_on_off) then
c natom is now the number of particles for which the ssbp forces
c will be computed: on each processor from istart to iend
c
		call load_balance(natom,my_pe,num_pes,istart,iend)
	end if

c     Definition of deltaRmaxDiel (DRM1) from the total charge CHTOT
c     and RMAX
	if (first) then
      		CHTOT=0.d0 
      		DO 17 I=1,npt
      			CHTOT = CHTOT + PTCHG(I)
   17 		CONTINUE
      		CHTOT=DABS(CHTOT)
      		COEF0=-0.5d0*kofdie*(DIECST-1.d0) 
		first = .false.
	end if
      DRM1=DRMAX1-CHTOT*1.6D0*DEXP(-0.5D0*RMAXM)
      DDRM1=(DRMAX1-DRM1)*0.5D0+1.d0
      RMAX1=RMAXM+DRM1
 
c     Definitions of some useful arrays for the fast calculation of
c     multipoles and forces

      DO 4 I=1,npt
Cdeb      rsxy = coor(1,i)*coor(1,i) + coor(2,i)*coor(2,i) 
Cdeb      RATOM2=RADIST(I)**2 
      rratom2 = 1.d0/ratom2(i)
      rratom3 = rrdist(i)*rratom2
      rrsxy = 1.d0/rsxy(i)
Cdeb      RATOM3=RATOM2*RADIST(I)
      tvari = coor(3,i)*rrdist(i)
      TVAR(I)=1.d0 
      DTVAR(I)=0.d0
      rch(i)=ptchg(i)
Cdeb      COM0=((0.D0,1.D0)*coor(2,i)+(1.D0,0.D0)*coor(1,i))/radist(i)
      COM0=((0.D0,1.D0)*coor(2,i)+(1.D0,0.D0)*coor(1,i))*rrdist(i)
      COM=(1.D0,0.D0)
      COMR(I)=1.D0
      COMI(I)=0.D0

      IF (NMULT.NE.0) THEN 

      DO 3 J=1,NMULT    
      II = NATOM*J+I
      IIM1=II-NATOM
      TVAR(II)=TVAR(IIM1)*TVARI
      RCH(II)=RCH(IIM1)*RADIST(I) 
Cdeb      TVAR(NATOM*J+I)=TVAR(NATOM*(J-1)+I)*TVARI 
Cdeb      RCH(NATOM*J+I)=RCH(NATOM*(J-1)+I)*RADIST(I) 
      COM=COM*COM0
Cdeb      COMR(NATOM*J+I)=DREAL(COM)
Cdeb      COMI(NATOM*J+I)=DIMAG(COM)
Cyef may cause problems on f90 due to bug in the compiler
      COMR(II)=DREAL(COM)
c      COMR(II)=REAL(COM)
      COMI(II)=DIMAG(COM)
    3 continue

      END IF

Cdeb  TVAR12(I)=TVARI/(1.D0-TVARI**2)
      RLR3(I)=RSXY(I)*RRATOM3  

      ZYR3(I)=coor(3,i)*coor(2,i)*RRATOM3 
      XRL(I)=coor(1,i)*RRSXY
      ZXR3(I)=coor(3,i)*coor(1,i)*RRATOM3 
      YRL(I)=coor(2,i)*RRSXY
      XR2(I)=coor(1,i)*RRATOM2 
      YR2(I)=coor(2,i)*RRATOM2 
      ZR2(I)=coor(3,i)*RRATOM2 

      TVAR12(I)=ZR2(I)/RLR3(I)
    4 CONTINUE

      IF(.NOT.QFIX)THEN
Cdeb      RMAXF=DDRM1/RMAX1/RMAXM
      	RMAXF=DDRM1/(RMAX1*RMAXM)
      ENDIF
      RMAXLU=RMAX1
Cdeb      RMAXS=RMAX1**2
      RMAXS=RMAX1*RMAX1

c     Main loop K=l-index in the multipole Qlm
      DO 10 I=1,NMULT+1
      RMAXLU=RMAXLU/RMAXS
      K=I-1
      N=K

c     Operations with Legendre's polynomes
c     Definition
      CALL POLE(N,APOL,FACT)
c     Duplication
      CALL EQPOL(L,APOL1,N,APOL)
c     Differentiation
      CALL DPOL(L,APOL1)

Cdeb      COEF=COEF0/(DIECST+DREAL(K/I)) 
      COEF=COEF0/(DIECST+DFLOAT(K/I)) 
      COEF=COEF/(2*K+1)
      COEF=COEF*RMAXLU
      IF(.NOT.QFIX)THEN
      COL=(2*K+1)*COEF*RMAXF
      COLX=COL*coor(1,IPA)
      COLY=COL*coor(2,IPA)
      COLZ=COL*coor(3,IPA)
      ENDIF
c     Enclosed Main loop M=m-index in the multipole Qlm
      DO 10 J=1,I
      ICOD=1
      M=J-1

      IF(M.GT.0)THEN
      CALL EQPOL(N,APOL,L,APOL1)
      CALL DPOL(L,APOL1)
      ICOD=2
      ENDIF

      CALL FQLM(K,M,NATOM,N,L,APOL,APOL1,FACT,
     &     QLMS,QLMSX,QLMSY,QLMSZ,TVAR12,RLR3,ZYR3,XRL,ZXR3,
     &     YRL,XR2,YR2,ZR2,COMR,COMI,RCH,TVAR,
     &     YLMR1,YLMR2,YLMR3,YLMR4,YLMI1,YLMI2,YLMI3,YLMI4)

c     Subroutine calculating Qlm**2=QLMS
c     and derivatives QLMSX,QLMSY,QLMSZ
 
      QLMS=ICOD*QLMS
      ENKIR=ENKIR+COEF*QLMS
      DO 16 ILP=1,NATOM

Cdeb      DX(ILP)=DX(ILP)+ICOD*QLMSX(ILP)*COEF
Cdeb      DY(ILP)=DY(ILP)+ICOD*QLMSY(ILP)*COEF
Cdeb      DZ(ILP)=DZ(ILP)+ICOD*QLMSZ(ILP)*COEF

      dpot(1,ILP)=dpot(1,ILP)+ICOD*QLMSX(ILP)*COEF
      dpot(2,ILP)=dpot(2,ILP)+ICOD*QLMSY(ILP)*COEF
      dpot(3,ILP)=dpot(3,ILP)+ICOD*QLMSZ(ILP)*COEF

   16 CONTINUE

      IF(.NOT.QFIX)THEN

Cdeb      DX(IPA)=DX(IPA)-COLX*QLMS
Cdeb      DY(IPA)=DY(IPA)-COLY*QLMS
Cdeb      DZ(IPA)=DZ(IPA)-COLZ*QLMS 

      dpot(1,IPA)=dpot(1,IPA)-COLX*QLMS
      dpot(2,IPA)=dpot(2,IPA)-COLY*QLMS
      dpot(3,IPA)=dpot(3,IPA)-COLZ*QLMS 

Cdeb      ELSE
      ENDIF
   10 CONTINUE  
      RETURN
      END 

c     ----------------------------------------------------
      SUBROUTINE FQLM(LIND,MIND,NATOM,NPOL,LPOL,
     &           APOL,APOL1,FACT,QLMS,QLMSX,QLMSY,QLMSZ,
     &           TVAR12,RLR3,ZYR3,XRL,ZXR3,YRL,XR2,YR2,ZR2,
     &           COMR,COMI,RCH,TVAR,
     &           YLMR1,YLMR2,YLMR3,YLMR4,
     &           YLMI1,YLMI2,YLMI3,YLMI4)
c     ----------------------------------------------------
c     Calculation of the multipole moment Qlm with
c     l=lind,m=mind
c      IMPLICIT NONE
      INTEGER LIND,MIND,NATOM,NPOL,LPOL 

      REAL*8 QLMS,APOL(*),APOL1(*),FACT(*)
      REAL*8 QLMSX(*),QLMSY(*),QLMSZ(*)
      REAL*8 TVAR12(*),RLR3(*),ZYR3(*),XRL(*),ZXR3(*)
      REAL*8 YRL(*),XR2(*),YR2(*),ZR2(*)
      REAL*8 COMI(*),COMR(*)
      REAL*8 RCH(*),TVAR(*)
      REAL*8 YLMR1(*),YLMR2(*),YLMR4(*),YLMR3(*)
      REAL*8 YLMI1(*),YLMI2(*),YLMI4(*),YLMI3(*)

c     Local variables
      INTEGER J,I
      REAL*8 BVAR,QLMR,QLMI,RC,DVAR,DVAR1,YLMIT,YLMRT
      REAL*8 YLMIC,YLMRC,DVAR2,BVAR1,BVARR1,BVARI1
 
      BVAR=(2*LIND+1)*FACT(LIND-MIND+1)/FACT(LIND+MIND+1)
      QLMR=0.d0
      QLMI=0.d0
      DO 10 I=1,NATOM
      RC=RCH(I+LIND*NATOM)
      DVAR=0.d0 
      DVAR1=0.d0
      DO 15 J=1,NPOL+1
   15 DVAR=DVAR+APOL(J)*TVAR(I+(NPOL+1-J)*NATOM) 
      DO 16 J=1,LPOL+1
   16 DVAR1=DVAR1+APOL1(J)*TVAR(I+(LPOL+1-J)*NATOM)
      YLMIT=COMI(I+NATOM*MIND) 
      YLMRT=COMR(I+MIND*NATOM)
      YLMIC=YLMIT*DVAR
      YLMRC=YLMRT*DVAR
      DVAR2=DVAR1-DVAR*TVAR12(I)*MIND
      BVAR1=DVAR2*RLR3(I)
      YLMR3(I)=YLMRT*BVAR1 
      YLMI3(I)=YLMIT*BVAR1
      BVARR1=-DVAR2*ZYR3(I)
      BVARI1=MIND*XRL(I) 
      YLMR2(I)=YLMRT*BVARR1-YLMIC*BVARI1
      YLMI2(I)=YLMIT*BVARR1+YLMRC*BVARI1
      BVARR1=-DVAR2*ZXR3(I)
      BVARI1=-MIND*YRL(I)
      YLMR1(I)=YLMRT*BVARR1-YLMIC*BVARI1
      YLMI1(I)=YLMIT*BVARR1+YLMRC*BVARI1
      YLMI4(I)=YLMIC
      YLMR4(I)=YLMRC
      QLMI=QLMI+RC*YLMI4(I) 
      QLMR=QLMR+RC*YLMR4(I)
   10 CONTINUE
      QLMS=QLMR*QLMR+QLMI*QLMI
      QLMS=QLMS*BVAR
      DO11 I=1,NATOM
      RC=2*BVAR*RCH(I+LIND*NATOM)
      BVAR1=LIND*XR2(I)
      BVARR1=BVAR1*YLMR4(I)+YLMR1(I)
      BVARI1=BVAR1*YLMI4(I)+YLMI1(I)
      QLMSX(I)=RC*(BVARR1*QLMR+BVARI1*QLMI) 
      BVAR1=LIND*YR2(I)
      BVARR1=BVAR1*YLMR4(I)+YLMR2(I)
      BVARI1=BVAR1*YLMI4(I)+YLMI2(I)
      QLMSY(I)=RC*(BVARR1*QLMR+BVARI1*QLMI)
      BVAR1=LIND*ZR2(I)
      BVARR1=BVAR1*YLMR4(I)+YLMR3(I)
      BVARI1=BVAR1*YLMI4(I)+YLMI3(I)
      QLMSZ(I)=RC*(BVARR1*QLMR+BVARI1*QLMI)
   11 CONTINUE
      RETURN
      END  

c     --------------------
      SUBROUTINE DPOL(N,A)
c     --------------------
c Differentiation of Legendre's polynome of N degree
c with coefficients A.
c      IMPLICIT NONE
      double precision a(*)
      INTEGER N
clocal
      integer I
 
      IF(N.GT.0)THEN
      DO 10 I=1,N
   10 A(I)=A(I)*(N-I+1) 
      N=N-1
      ELSE
      N=0
      A(1)=0.D0
      ENDIF
      RETURN
      END

c     -----------------------------
      SUBROUTINE EQPOL(N1,A1,N2,A2)
c     -----------------------------
c Duplicate a Legendre's polynome of N2 degree with coefficients A2
c      IMPLICIT NONE
      double precision A1(*),A2(*)
      INTEGER N1,N2
clocal
      integer I

      DO 10 I=1,N2+1 
      A1(I)=A2(I)
   10 CONTINUE
      N1=N2
      RETURN
      END

c     -------------------------------
      SUBROUTINE POLE(NPOL,APOL,FACT)
c     -------------------------------
c Calculation of the coefficients (APOL) of the Legendre's
c polynome of the degree NPOL
c      IMPLICIT NONE
      double precision APOL(*),FACT(*)
      integer NPOL
clocal
      double precision BARG 
      INTEGER IBARG,NHAL
      INTEGER I,J,K,L,M  
 
Cdeb     BARG=DREAL(NPOL)/2.D0 
Cdeb      IBARG=INT(BARG)
Cdeb      IF(DABS(BARG-DREAL(IBARG)).GT.0.1D0)THEN
      IF ((NPOL/2*2).NE.NPOL) THEN
      NHAL=(NPOL-1)/2
      ELSE
      NHAL=NPOL/2
      ENDIF
      DO 5 I=1,NPOL+1
    5 APOL(I)=0.D0
      DO 10 J=1,NHAL+1
      M=J-1
      K=NPOL-M
      L=2*K+1 
      BARG=FACT(L)/FACT(J)/FACT(K+1)/FACT(L-NPOL)
      BARG=BARG/(2.D0)**NPOL
   10 APOL(1+2*M)=BARG*(-1)**M
      RETURN
      END

c     -----------------------------------------------------
      SUBROUTINE SSBHSRP
c     -----------------------------------------------------

c Calculation of Pressure*Volume+Surface_tension*Surface term
c Pressure=PRESI, Surface_tension=STENS

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
C     local
	logical first
      double precision aval,bval,val,pi
	save pi,first
	data first/.true./

	if (first) then
      		pi = (4.d0*datan(1.d0))
		first = .false.
	end if
      AVAL=ATMOSP*PRESI*RMAXM**3*PI*4.D0/3.D0
      BVAL=STENS*4.d0*PI*RMAXM*RMAXM
      ENHSR=AVAL+BVAL
      IF(.NOT.QFIX)THEN
      VAL=(3.D0*AVAL+2.D0*BVAL)/RMAXM
Cdeb      DX(IPA)=DX(IPA)+VAL*X(IPA)/RADIST(IPA)
Cdeb      DY(IPA)=DY(IPA)+VAL*Y(IPA)/RADIST(IPA)
Cdeb      DZ(IPA)=DZ(IPA)+VAL*Z(IPA)/RADIST(IPA) 
      dpot(1,IPA)=dpot(1,IPA)+VAL*coor(1,IPA)*RRDIST(IPA)
      dpot(2,IPA)=dpot(2,IPA)+VAL*coor(2,IPA)*RRDIST(IPA)
      dpot(3,IPA)=dpot(3,IPA)+VAL*coor(3,IPA)*RRDIST(IPA) 
Cdeb      ELSE
      ENDIF

      RETURN
      END

c     -----------------------------------------------------
      SUBROUTINE SSBCAVP
c     ------------------------------------------------------

c Calculation of the Cavity part.
c This constribution was calculated from the cavity potential of mean force
c between a LJ TIP3P-like oxygen atom and a large hard-sphere of radius RMAX.
c The polynomial approximation has been fitted using mathematica and yields
c a boundary potential that is very similar to the SBOUND method of C.L. Brooks.

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

Clocal     
      integer i,j 
      double precision rmax2,shi,dshi,xarg,pot2,dpot2,dut 

Cdeb      RMAX2=RMAXM+DRCAV-2.6d0
      RMAX2=RMAXM+DRMAX2-2.6d0
      CALL RARA(RMAX2,SHI,DSHI,ACAV) 
Cdeb      ESSBP2=ZERO
      DUT=0.d0
      DO 10 I=1,NTSSBP
      J=LSTSSBP(I)
      XARG=-RMAX2+RADIST(J)
      CALL COMPOT(XARG,POT2,DPOT2,BCAV)
      ENCAV=ENCAV+POT2

Cdeb      DX(J)=DX(J)+DPOT2*X(J)/RADIST(J)
Cdeb      DY(J)=DY(J)+DPOT2*Y(J)/RADIST(J)
Cdeb      DZ(J)=DZ(J)+DPOT2*Z(J)/RADIST(J)

      dpot(1,J)=dpot(1,J)+DPOT2*coor(1,J)*RRDIST(J)
      dpot(2,J)=dpot(2,J)+DPOT2*coor(2,J)*RRDIST(J)
      dpot(3,J)=dpot(3,J)+DPOT2*coor(3,J)*RRDIST(J)

      DUT=DUT+DPOT2 
   10 CONTINUE
      IF(.NOT.QFIX)THEN
      DSHI=DSHI*NTSSBP-DUT

Cdeb      DX(IPA)=DX(IPA)+DSHI*X(IPA)/RADIST(IPA)
Cdeb      DY(IPA)=DY(IPA)+DSHI*Y(IPA)/RADIST(IPA)
Cdeb      DZ(IPA)=DZ(IPA)+DSHI*Z(IPA)/RADIST(IPA) 

      dpot(1,IPA)=dpot(1,IPA)+DSHI*coor(1,IPA)*RRDIST(IPA)
      dpot(2,IPA)=dpot(2,IPA)+DSHI*coor(2,IPA)*RRDIST(IPA)
      dpot(3,IPA)=dpot(3,IPA)+DSHI*coor(3,IPA)*RRDIST(IPA) 
Cdeb      ELSE
      ENDIF
      ENCAV=ENCAV+SHI*NTSSBP
      RETURN
      END

c     --------------------------------------
      SUBROUTINE RARA(XARGA,CAPA,DCAPA,ACAV)
c     --------------------------------------
c Calculation of the Function-axis shift of the approximated cavity potential
c      IMPLICIT NONE
      double precision XARGA,CAPA,DCAPA,ACAV(*)
Clocal
      double precision xarga2,xarga3,rama 
ca0   ACAV(5)=-1.6649500d0
ca1   ACAV(1)=0.56198800d0
ca2   ACAV(2)=-0.072798148d0
ca3   ACAV(3)=0.00426122036d0
ca4   ACAV(4)=-0.0000925233817d0
cac   ACAV(6)=0.0840d0
cXc   ACAV(7)=15.39333D0
      RAMA=XARGA+2.6d0
      IF(RAMA.LT.5.D0)THEN
      WRITE (*,'(A)')' SSBP CAVITY: RMAX+DRCAV<MIN  APPROXIMATED'
      WRITE (*,'(A,F8.4)')
     &       '         MIN = 5.0 A, RMAX+DRCAV = ', RAMA
      ENDIF
      XARGA2=XARGA**2 
      XARGA3=XARGA2*XARGA 
      IF(XARGA.GT.ACAV(7))THEN
      CAPA=ACAV(6)  
      DCAPA=0.d0
      ELSE
      CAPA=ACAV(5)+ACAV(1)*XARGA+ACAV(2)*XARGA2+ACAV(3)*XARGA3+
     & ACAV(4)*XARGA*XARGA3
      DCAPA=ACAV(1)+2*XARGA*ACAV(2)+3*ACAV(3)*XARGA2+4*ACAV(4)*XARGA3
      ENDIF
      CAPA=CAPA+8.5d0
      RETURN
      END

c     ----------------------------------------
      SUBROUTINE COMPOT(XARGB,CAPB,DCAPB,BCAV)
c     ----------------------------------------
c Approximated calculation of the cavity potential
c      IMPLICIT NONE
      double precision XARGB,XARGB2,CAPB,DCAPB,BCAV(*)

c     Coefficients:
c     BCAV(1)= 1.319978287d0
c     BCAV(2)=-0.840953501d0
c     BCAV(3)=-0.001602388122d0
c     BCAV(4)=-8.392886499d0
c     BCAV(5)=BCAV(2)+BCAV(4)
c     BCAV(6)= 1.6d0
c     BCAV(7)=-8.4751210228d0
      IF(XARGB.LT.-5.D0)THEN
      CAPB=BCAV(7) 
      DCAPB=0.D0
      ELSE
      XARGB2=XARGB**2
      IF(XARGB.GT.0.D0)THEN
      CAPB=BCAV(5)+BCAV(6)*XARGB2 
      DCAPB=2.D0*XARGB*BCAV(6) 
      ELSE
      CAPB=BCAV(2)/(1.D0+XARGB2/BCAV(1))+BCAV(3)*XARGB2+BCAV(4) 
      DCAPB= 2.D0*XARGB*
     & (BCAV(3)-BCAV(2)/BCAV(1)/(1.D0+XARGB2/BCAV(1))**2)
      ENDIF
      ENDIF
      RETURN
      END

c     -------------------------------------
      SUBROUTINE SSBANGP
c     ----------------------------

c Calculation of the Angular Potential 

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

c     Local variables:
      INTEGER I,J 
      double precision RIPA,DUTIPA,RADLOC,VARDP,VANG,FRCNST,DPOT5,POT5     
      double precision VALO(3),DVALO(3),DVALH(2,3),VALH(2,3)


Cdeb  RIPA=RMAX5 
      RIPA=RMAXM
      IF(.NOT.QFIX)RIPA=RADIST(IPA)
Cdeb      SSBPE5=ZERO
      DUTIPA=0.d0
      DO 10 I=1,NTSSBP
      J=LSTSSBP(I)
      RADLOC=RADIST(J)
      VARDP=RADLOC-RIPA-DRHA

Cdeb      VALO(1)=X(J)
Cdeb      VALO(2)=Y(J)
Cdeb      VALO(3)=Z(J)
Cdeb      VALH(1,1)=X(J+1)
Cdeb      VALH(1,2)=Y(J+1)
Cdeb      VALH(1,3)=Z(J+1)
Cdeb      VALH(2,1)=X(J+2)
Cdeb      VALH(2,2)=Y(J+2)
Cdeb      VALH(2,3)=Z(J+2)

      VALO(1)=coor(1,J)
      VALO(2)=coor(2,J)
      VALO(3)=coor(3,J)
      VALH(1,1)=coor(1,J+1)
      VALH(1,2)=coor(2,J+1)
      VALH(1,3)=coor(3,J+1)
      VALH(2,1)=coor(1,J+2)
      VALH(2,2)=coor(2,J+2)
      VALH(2,3)=coor(3,J+2)

      IF(VARDP.GT.0.D0)THEN
      CALL ANGLAIP(rrdist(i),VALO,VALH,DVALO,DVALH,VANG,CANG)
      DPOT5=VARDP*VANG
      DUTIPA=DUTIPA+DPOT5
      POT5=0.5D0*DPOT5*VARDP  
      FRCNST=0.5D0*VARDP**2
      ENANGU=ENANGU+POT5 

Cdeb      DX(J)=DX(J)+DPOT5*VALO(1)/RADLOC+FRCNST*DVALO(1)
Cdeb      DY(J)=DY(J)+DPOT5*VALO(2)/RADLOC+FRCNST*DVALO(2)
Cdeb      DZ(J)=DZ(J)+DPOT5*VALO(3)/RADLOC+FRCNST*DVALO(3)
Cdeb      DX(J+1)=DX(J+1)+FRCNST*DVALH(1,1)
Cdeb      DY(J+1)=DY(J+1)+FRCNST*DVALH(1,2)
Cdeb      DZ(J+1)=DZ(J+1)+FRCNST*DVALH(1,3)
Cdeb      DX(J+2)=DX(J+2)+FRCNST*DVALH(2,1)
Cdeb      DY(J+2)=DY(J+2)+FRCNST*DVALH(2,2)
Cdeb      DZ(J+2)=DZ(J+2)+FRCNST*DVALH(2,3)

      dpot(1,J)=dpot(1,J)+DPOT5*VALO(1)*rrdist(j)+FRCNST*DVALO(1)
      dpot(2,J)=dpot(2,J)+DPOT5*VALO(2)*rrdist(j)+FRCNST*DVALO(2)
      dpot(3,J)=dpot(3,J)+DPOT5*VALO(3)*rrdist(j)+FRCNST*DVALO(3)
      dpot(1,J+1)=dpot(1,J+1)+FRCNST*DVALH(1,1)
      dpot(2,J+1)=dpot(2,J+1)+FRCNST*DVALH(1,2)
      dpot(3,J+1)=dpot(3,J+1)+FRCNST*DVALH(1,3)
      dpot(1,J+2)=dpot(1,J+2)+FRCNST*DVALH(2,1)
      dpot(2,J+2)=dpot(2,J+2)+FRCNST*DVALH(2,2)
      dpot(3,J+2)=dpot(3,J+2)+FRCNST*DVALH(2,3)
      ENDIF
   10 CONTINUE
      IF(.NOT.QFIX)THEN

Cdeb      DX(IPA)=DX(IPA)-DUTIPA*X(IPA)/RADIST(IPA)
Cdeb      DY(IPA)=DY(IPA)-DUTIPA*Y(IPA)/RADIST(IPA)
Cdeb      DZ(IPA)=DZ(IPA)-DUTIPA*Z(IPA)/RADIST(IPA)

      dpot(1,IPA)=dpot(1,IPA)-DUTIPA*coor(1,IPA)*RRDIST(IPA)
      dpot(2,IPA)=dpot(2,IPA)-DUTIPA*coor(2,IPA)*RRDIST(IPA)
      dpot(3,IPA)=dpot(3,IPA)-DUTIPA*coor(3,IPA)*RRDIST(IPA)
Cdeb      ELSE
      ENDIF
      RETURN
      END

c     ---------------------------------------------------------- 
      SUBROUTINE ANGLAIP(RRDIST,VALO,VALH,DVALO,DVALH,VANG,CANG)
c     ----------------------------------------------------------
c Calculation of the angular part of the angular potential
c This contribution was developed empirically to make more isotropic the 
c orientational distribution function of the waters near the outer shell.
c Warning:  It works only for a 3 site water models such as TIP3P.
c The order of the atoms in the psf MUST BE (oxygen, hydrogen,hydrogen)
      INTEGER I,J
      REAL*8 OHVEC(2,3),RRDIST,VALO(3),DVALO(3),DVALH(2,3),VALH(2,3)
      REAL*8 EXAN(2),REBMUN1(2),REBMUN2(2),REBMUN3(2),CO(2),DEXAN(2) 
      REAL*8 CO1,CO2,CO3,CO4,VANG,CANG(*),SCALMUL,OHABS 
c Coeffitients:
c CANG(5)=0.840661d0
c CANG(4)=- 1.20064327666d0
c CANG(3)=- 3.067139576698d0
c CANG(2)=1.766839115555d0
c CANG(1)= 2.4085919002d0
      DO 30 J=1,2
      SCALMUL=0.d0
      DO 10 I=1,3
      OHVEC(J,I)=VALH(J,I)-VALO(I)
      SCALMUL=SCALMUL+VALO(I)*OHVEC(J,I)
   10 CONTINUE
      OHABS=OHVEC(J,1)**2+OHVEC(J,2)**2+OHVEC(J,3)**2
      REBMUN1(J)=RRDIST/DSQRT(OHABS)  
      REBMUN2(J)=SCALMUL/OHABS
      REBMUN3(J)=SCALMUL*RRDIST*RRDIST
      CO(J)=SCALMUL*REBMUN1(J) 
      CO1=CO(J)
      CO2=CO1*CO(J)
      CO3=CO2*CO(J)
      CO4=CO3*CO(J)
      DEXAN(J)=2*
     & (4.D0*CANG(1)*CO3+3.D0*CANG(2)*CO2+2.D0*CANG(3)*CO1+CANG(4))
      EXAN(J)=2*
     & (CANG(1)*CO4+CANG(2)*CO3+CANG(3)*CO2+CANG(4)*CO1+CANG(5))
   30 CONTINUE
c     
      DO 15 I=1,3
      DVALH(1,I)=DEXAN(1)*REBMUN1(1)*(-REBMUN2(1)*OHVEC(1,I)+VALO(I))
      DVALH(2,I)=DEXAN(2)*REBMUN1(2)*(-REBMUN2(2)*OHVEC(2,I)+VALO(I)) 
      DVALO(I)=DEXAN(1)*REBMUN1(1)*
     & (-VALO(I)*(REBMUN3(1)+1.D0)+OHVEC(1,I)*(1.D0+REBMUN2(1)))
      DVALO(I)=DVALO(I)+DEXAN(2)*REBMUN1(2)*
     & (-VALO(I)*(REBMUN3(2)+1.D0)+OHVEC(2,I)*(1.D0+REBMUN2(2)))
   15 CONTINUE
c    
      VANG=EXAN(1)+EXAN(2) 
      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012
