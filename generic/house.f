      SUBROUTINE HOUSE (Z,NZ,N,D,E,ERROR)
C       Diagonalization routine written by Ryszard Czerminski
c based on householder algorithm
c Input: Z a real double precision matrix
c        NZ=N  matrix dimension
c output: Z eigenvector matrix
c         Z(i=1,n;j) is the j-th eigenvector
c         D vecotor of length N with eigenvalues
c         error, integer hopefully zero on output
c
      INTEGER ERROR,N,NZ
      double precision MACHEP
      double precision D(N),E(N),Z(NZ,N)
      integer ii,i,j,k,l,m,ip1,jp1,mml,nm1
	double precision b,c,f,g,h,p,r,s,small,tol,HH
c       DATA SMALL/5.4D-79/      IBM 370
c       MACHEP=2.D0**(-56)       IBM 370
C       DATA SMALL/2.23D-308/    IBM PC real*8 ms-fortran 4.0
C       MACHEP=2.D0**(-52)       ibm pc real*8 ms-fortran 4.0
        DATA SMALL/2.225D-35/
C       DSQRT(X)=DDSQRT(X)
C       DSIGN(X,Y)=DDSIGN(X,Y)
C       DABS(X)=DDABS(X)
        MACHEP=2.D0**(-52)
C
      TOL=SMALL/MACHEP
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 2
         F = Z(I-1,I)
         G = 0.0D0
      IF (L .LT. 1) GO TO 140
      DO 120 K = 1, L
  120    G = G + Z(K,I) * Z(K,I)
  140    H = G + F * F
         IF (G .GT. TOL) GO TO 160
      E(I) = F
      H = 0.0D0
      GO TO 280
  160    L = L + 1
         G = -DSIGN(DSQRT(H),F)
         E(I) = G
         H = H - F * G
         Z(I-1,I) = F - G
         F = 0.0D0
         DO 240 J = 1, L
            Z(I,J) = Z(J,I)
            Z(J,I) = Z(J,I) / H
            G = 0.0D0
            DO 180 K = 1, J
  180       G = G + Z(K,J) * Z(I,K)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
  200      G = G + Z(J,K) * Z(K,I)
  220       E(J) = G / H
            F = F + G * Z(J,I)
  240      CONTINUE
            HH = F / (H + H)
            DO 260 J = 1, L
              F = Z(I,J)
              G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(K,J) = Z(K,J) - F * E(K) - G * Z(I,K)
  260      CONTINUE
  280      D(I) = H
  300 CONTINUE
  320 D(1) = 0.0D0
      E(1) = 0.0D0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
         DO 360 J = 1, L
            G = 0.0D0
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
  380     D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
  500 CONTINUE
      ERROR = 0
      IF (N .EQ. 1) GO TO 1500
      DO 1100 I = 2, N
 1100 E(I-1) = E(I)
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
      DO 1240 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H
         DO 1110 M = L, N
            IF (DABS(E(M)) .LE. B) GO TO 1120
 1110      CONTINUE
 1120      IF (M .EQ. L) GO TO 1220
 1130      IF (J .EQ. 30) GO TO 1400
         J = J + 1
          P = (D(L+1) - D(L)) / (2.0D0 * E(L))
          R = DSQRT(P * P + 1.0D0)
          H = D(L) - E(L) / (P + DSIGN(R,P))
          DO 1140 I = L, N
 1140    D(I) = D(I) - H
          F = F + H
          P = D(M)
          C = 1.0D0
          S = 0.0D0
          MML = M - L
          DO 1200 II = 1, MML
            I = M - II
            G = C * E(I)
      H = C * P
      IF (DABS(P) .LT. DABS(E(I))) GO TO 1150
            C = E(I) / P
      R = DSQRT(C * C + 1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 1160
 1150        C = P / E(I)
            R = DSQRT(C * C + 1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C / R
 1160        P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 1180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
 1180        CONTINUE
 1200    CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (DABS(E(L)) .GT. B) GO TO 1130
 1220    D(L) = D(L) + F
 1240 CONTINUE
      NM1 = N - 1
      DO 1300 I = 1, NM1
         K = I
         P = D(I)
            IP1 = I + 1
         DO 1260 J = IP1, N
            IF (D(J) .LE. P) GO TO 1260
            K = J
            P = D(J)
 1260    CONTINUE
          IF (K .EQ. I) GO TO 1300
          D(K) = D(I)
          D(I) = P
          DO 1280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
 1280    CONTINUE
 1300 CONTINUE
      GO TO 1500
 1400 ERROR = L
 1500 RETURN
      END
