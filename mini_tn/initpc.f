      SUBROUTINE INITPC(DIAGB,EMAT,N,W,LW,MODET,
     *     UPD1,YKSK,GSK,YRSR,LRESET)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,LW,MODET,LZ1,LZK,LV,LSK,LYK,LDIAGB,LGV
      INTEGER LSR,LYR,LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
      DOUBLE PRECISION DIAGB(N),EMAT(N),W(LW)
      DOUBLE PRECISION YKSK,GSK,YRSR
      LOGICAL LRESET,UPD1
      COMMON/SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
      CALL INITP3(DIAGB,EMAT,N,LRESET,YKSK,YRSR,W(LHYK),
     *     W(LSK),W(LYK),W(LSR),W(LYR),MODET,UPD1)
      RETURN
      END
