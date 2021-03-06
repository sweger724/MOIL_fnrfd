      SUBROUTINE TN (IERROR, N, X, F, G, W, LW,
     1      SFUN, EFCALL, MAXFUN, XTOL, STEPMX)
c      IMPLICIT          DOUBLE PRECISION (A-H,O-Z)
      INTEGER           IERROR, N, LW, EFCALL
      integer maxfun,maxit,msglvl,nmax,i
      DOUBLE PRECISION  X(N), G(N), F, W(LW)
C
C THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM
C
C            MINIMIZE F(X)
C               X
C
C WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS
C A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA
C THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984),
C PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES
C NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A
C GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW.
C IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS
C ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE.
C
C SUBROUTINE PARAMETERS:
C
C IERROR - (INTEGER) ERROR CODE
C          ( 0 => NORMAL RETURN)
C          ( 2 => MORE THAN MAXFUN EVALUATIONS)
C          ( 3 => LINE SEARCH FAILED TO FIND
C          (          LOWER POINT (MAY NOT BE SERIOUS)
C          (-1 => ERROR IN INPUT PARAMETERS)
C N      - (INTEGER) NUMBER OF VARIABLES
C X      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL
C          ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION.
C G      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL
C          VALUE OF THE GRADIENT
C F      - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE
C          OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE
C          OF THE OBJECTIVE FUNCTION AT THE SOLUTION
C W      - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N
C LW     - (INTEGER) THE DECLARED DIMENSION OF W
C SFUN   - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION
C          AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE
C          THE CALLING SEQUENCE
C             SUBROUTINE SFUN (N, X, F, G)
C             INTEGER           N
C             DOUBLE PRECISION  X(N), G(N), F
C
C THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE
C LMQN.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE
C OF THIS ALGORITHM SHOULD CALL LMQN DIRECTLY.
C
C----------------------------------------------------------------------
C THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON
C ALGORITHM.  THE PARAMETERS ARE:
C
C ETA    - SEVERITY OF THE LINESEARCH
C MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
C XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
C STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
C ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
C MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
C          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
C MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP
C
      DOUBLE PRECISION ETA, ACCRCY, XTOL, STEPMX, MCHPR1
      EXTERNAL         SFUN
C
C SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE
C
      MAXIT = N/2
      IF (MAXIT .GT. 50) MAXIT = 50
      IF (MAXIT .LE. 0) MAXIT = 1
      MSGLVL = 1
C      MAXFUN = 150*N
      ETA = .25D0
C      STEPMX = 1.D1
      ACCRCY = 1.D2*MCHPR1()
C      XTOL = DSQRT(ACCRCY)
C
C MINIMIZE THE FUNCTION
C
      CALL LMQN (IERROR, N, X, F, G, W, LW, SFUN,
     *     MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL, EFCALL)
C
C PRINT THE RESULTS
C
      IF (IERROR .NE. 0) WRITE(*,800) IERROR
      WRITE(*,810) F
      IF (MSGLVL .LT. 1) RETURN
      WRITE(*,820)
      NMAX = 10
      IF (N .LT. NMAX) NMAX = N
      WRITE(*,830) (I,X(I),I=1,NMAX)
      RETURN
800   FORMAT(//,' ERROR CODE =', I3)
810   FORMAT(//,' OPTIMAL FUNCTION VALUE = ', 1PD22.15)
820   FORMAT(10X, 'CURRENT SOLUTION IS (AT MOST 10 COMPONENTS)', /,
     *       14X, 'I', 11X, 'X(I)')
830   FORMAT(10X, I5, 2X, 1PD22.15)
      END
