        subroutine invmat(c,n,dtnrm,detm)
        include 'COMMON/LENGTH.BLOCK'
        integer l,n,k,m,j(mcoupc+50)
        real*8 pd,dd,detm,cc,s,dtnrm,c(mcoupc,mcoupc)
        PD=1.D00
        DO 124 L=1,N
                DD=0.D00
                DO 123 K=1,N
123             DD=DD+C(L,K)*C(L,K)
                DD=DSQRT(DD)
124             PD=PD*DD
        DETM=1.D00
        DO 125 L=1,N
125     J(L+20)=L
        DO 144 L=1,N
                CC=0.D00
                M=L
                DO 135 K=L,N
                        IF((DABS(CC)-DABS(C(L,K))).GE.0.D00) GO TO 135
126                      M=K
                        CC=C(L,K)
135             CONTINUE
127             IF (L.EQ.M) GO TO 138
128             K=J(M+20)
                J(M+20)=J(L+20)
                J(L+20)=K
                DO 137 K=1,N
                        S=C(K,L)
                        C(K,L)=C(K,M)
137             C(K,M)=S
138             C(L,L)=1.D00
                DETM=DETM*CC
                DO 139 M=1,N
139             C(L,M)=C(L,M)/CC
                DO 142 M=1,N
                        IF(L.EQ.M) GO TO 142
129                     CC=C(M,L)
                        IF (CC.EQ.0.D00) GO TO 142
130                     C(M,L)=0.D00
                        DO 141 K=1,N
141                     C(M,K)=C(M,K)-CC*C(L,K)
142             CONTINUE
144     CONTINUE
        DO 143 L=1,N
                IF (J(L+20).EQ.L) GO TO 143
131             M=L
132             M=M+1
                IF(J(M+20).EQ.L) GO TO 133
136             IF (N.GT.M) GO TO 132
133             J(M+20)=J(L+20)
                DO 163 K=1,N
                        CC=C(L,K)
                        C(L,K)=C(M,K)
163             C(M,K)=CC
                J(L+20)=L
143     CONTINUE
        DETM=DABS(DETM)
        DTNRM=DETM/PD
        RETURN
        END
