      SUBROUTINE ZTIME(N,X,IPIVOT)
c      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      integer n,i
      DOUBLE PRECISION X(N)
      INTEGER          IPIVOT(N)
C
C THIS ROUTINE MULTIPLIES THE VECTOR X BY THE CONSTRAINT MATRIX Z
C
      DO 10 I = 1,N
         IF (IPIVOT(I) .NE. 0) X(I) = 0.D0
10    CONTINUE
      RETURN
      END
