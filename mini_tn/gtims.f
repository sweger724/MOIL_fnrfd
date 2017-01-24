      SUBROUTINE GTIMS(V,GV,N,X,G,W,LW,SFUN,FIRST,DELTA
     * ,ACCRCY,XNORM,EFCALL)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER EFCALL,N,LW
      INTEGER LGV,LZ1,LZK,LYK,LDIAGB,LSR,LYR,LHYR,LHG
      INTEGER LV,LSK,LHYK,LPK,LEMAT,LWTEST,IHG,I
      DOUBLE PRECISION V(N),GV(N),DINV,DELTA,G(N)
      DOUBLE PRECISION F,X(N),W(LW),ACCRCY,DSQRT,XNORM
      LOGICAL FIRST
      EXTERNAL SFUN
      COMMON/SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
C
C THIS ROUTINE COMPUTES THE PRODUCT OF THE MATRIX G TIMES THE VECTOR
C V AND STORES THE RESULT IN THE VECTOR GV (FINITE-DIFFERENCE VERSION)
C
      IF (.NOT. FIRST) GO TO 20
      DELTA = DSQRT(ACCRCY)*(1.D0+XNORM)
      FIRST = .FALSE.
20    CONTINUE
      DINV = 1.D0/DELTA
      IHG = LHG
      DO 30 I = 1,N
         W(IHG) = X(I) + DELTA*V(I)
         IHG = IHG + 1
30    CONTINUE
      CALL SFUN(N,W(LHG),F,GV,EFCALL)
      DO 40 I = 1,N
         GV(I) = (GV(I) - G(I))*DINV
40    CONTINUE
      RETURN
      END