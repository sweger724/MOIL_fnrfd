       subroutine CORFUN(NCOR,NRUN,DT,EI,EJ,ECF,intCF,ffluc)

C subroutine to compute the Correlation Function
C (taken mostly from Allen-Tildeslay)
C input:
C 1) NCOR = number of time steps on which the correlation function is computed
C 2) NRUN = total number of steps of the simulation
C 3) DT   = time interval between two configurations used
C 4) EI   = first variable in CORFUN
C 5) EJ   = second variable in CORFUN (ie <EI(0)EJ(t)> is computed)
C 6) ECF  = correlation function
C 7) intCF= integral of the correlation function 
C 8) ffluc= if true, compute <EI(0)EJ(t)> - <EI><EJ>, otherwise <EI(0)EJ(t)>

       implicit none

       integer lstr,i,j,t,NRUN,NCOR,EMAX
       double precision DT,EIAVE,EJAVE,EI(NRUN),EJ(NRUN)
       double precision dE0,dEI(NRUN),dEJ(NRUN)
       double precision ECF(NCOR)
       double precision NCF(NCOR),EVAR_R
       double precision intCF
       logical ffluc
       character*80 str

c compute the averages
       EIAVE = 0.d0
       EJAVE = 0.d0
       if (ffluc) then
        DO i=1,NRUN
         EIAVE = EIAVE + EI(i)
         EJAVE = EJAVE + EJ(i)
        ENDDO
        EIAVE = EIAVE/dble(NRUN)
        EJAVE = EJAVE/dble(NRUN)
       endif
c obtain the dEI(i) and dEJ(i) 
       DO i=1,NRUN
        dEI(i) = EI(i) - EIAVE
        dEJ(i) = EJ(i) - EJAVE
       ENDDO
c compute the correlation function <dE(i)dE(0)>
       DO i=1,NCOR
        ECF(i) = 0.d0
        NCF(i) = 0.d0
       ENDDO
       DO i=1,NRUN
        dE0 = dEI(i)
        EMAX = MIN(NRUN,i+NCOR)
        DO j=i,EMAX
         t = j-i+1
         ECF(t) = ECF(t) + dE0*dEJ(j)
         NCF(t) = NCF(t) + 1.d0
        ENDDO
       ENDDO
       DO i=1,NCOR
        ECF(i) = ECF(i)/NCF(i)
       ENDDO
       intCF = 0.d0

c compute the integral 
       DO i=1,NCOR-1
        intCF = intCF + (ECF(i)+ECF(i+1))*DT*0.5d0
       ENDDO

       return
       end
