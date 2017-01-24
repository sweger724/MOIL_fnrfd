      SUBROUTINE SSBFGS(N,GAMMA,SJ,YJ,HJV,HJYJ,YJSJ,YJHYJ,
     *     VSJ,VHYJ,HJP1V)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION GAMMA,YJSJ,YJHYJ,VSJ,VHYJ
      DOUBLE PRECISION SJ(N),YJ(N),HJV(N),HJYJ(N),HJP1V(N)
C
C SELF-SCALED BFGS
C
      INTEGER I
      DOUBLE PRECISION BETA,DELTA
      DELTA = (1.D0 + GAMMA*YJHYJ/YJSJ)*VSJ/YJSJ
     *     - GAMMA*VHYJ/YJSJ
      BETA = -GAMMA*VSJ/YJSJ
      DO 10 I = 1,N
         HJP1V(I) = GAMMA*HJV(I) + DELTA*SJ(I) + BETA*HJYJ(I)
10    CONTINUE
      RETURN
      END
