        subroutine eCG_backboneNB(T,r,A,B,C,E,dEdr,dEdA,dEdB,dEdC,deriv)
C       HB backbone potential according to
C      J. Phys. Chem. B 2004, 108, 9421 - 9438
C      Parametrixation of Backbone-Electrostatic and Multibody Contributions to the 
C      UNRES Force Field for Protein-Structure Prediction from Ab Initio Energy Surfaces 
C      of Model Systems;  Liwo, Oldziej, Czaolewski, Kozlowska, Scheraga 

        implicit none

        double precision r,A,B,C,E,dEdr,dEdA,dEdB,dEdC
        integer T
        logical deriv

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision s,s6,s8,s12,dd

C        if (r.gt. rpp(T)) then 
            s = rpp(T)/(3.5*r-2.5*rpp(T))
            dd = 3.5d0
C        else
C            s = rpp(T)/(0.8*r+0.2*rpp(T))
C            dd = 0.8d0
C        end if

        
        s6 = s**6
        s12 = s6**2
        s8 = s**8

        E = App(T) * (A - 3.d0*B*C)/(r**3)
     &    - Bpp(T) * (4.d0 + (A-3.d0*B*C)**2 
     &               - 3.d0 * (B**2 + C**2))/(r**6)
     &    +  epp(T) * (s12 - 2.d0*s6)

        if(deriv) then
        
        dEdr = -3.d0 * App(T) * (A - 3.d0*B*C)/(r**4)
     &       +  6.d0 * Bpp(T) * (4.d0 + (A-3.d0*B*C)**2
     &                          - 3.d0 * (B**2 + C**2))/(r**7)
     &        - dd * 12.d0 * epp(T) * (s12 - s6)*s/rpp(T)
        

        dEdA = App(T)/(r**3) - Bpp(T)*(2.d0*(A-3.d0*B*C))/(r**6)

        dEdB = -3.d0*App(T)*C/(r**3) 
     &       + 6.d0*Bpp(T)*((A - 3.d0*B*C) * C + B)/(r**6)

        dEdC = -3.d0*App(T)*B/(r**3)
     &       + 6.d0*Bpp(T)*((A - 3.d0*B*C) * B + C)/(r**6)


        endif

        
C       if(E .gt. E_CG_max/wCG_HB) then
C           write(6,*)"E_HB: r,E,dEdr,dEdA,dEdB,dEdC: ",
C     &                r,E*wCG_HB,dEdr,dEdA,dEdB,dEdC
C       endif
        return
        end
