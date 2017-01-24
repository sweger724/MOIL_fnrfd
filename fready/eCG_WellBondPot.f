      subroutine eCG_WellBondPot(x,ibond,E,dE)
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c calculate bond energy and forces for i-th bond.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision x,db,dE,E, dG1, dG2
      integer ibond
      double precision Mix, dE1, dE2, dE3, dE12, E1, E2, E3, E12

c  DATASTRUCTURE CGBond
c  column   1       2      3      4      5        6        7    8        9      10      11   12
c  data    type   r_eq1  e_min1   k_1   r_eq2   e_min_2   k_2  beta1   r_eq3  e_min3   k_3   beta2
c
c  type encodes how many wells are there (1,2,3)
c
c   for type==1     E = 1/2 k_1 (r - r_eq1)^2
c
c   for type==2     E_1 = 1/2 k_1 (r - r_eq1)^2  + e_min1
c                   E_2 = 1/2 k_2 (r - r_eq2)^2  + e_min2
c                   E = 1/2 * ( E_1 + E_2 - sqrt((E_1-E_2)^2 + beta_1) )
c
c   for type==3    E_12 = E (from type==2)
c                  E_3 = 1/2 k_3 (r - r_eq3)^2  + e_min3
c                  E = 1/2 * ( E_12 + E_3 - sqrt((E_12-E_3)^2 + beta_2) )


        if (CGBond(ibond,1).eq.1) then
          db = x - CGBond(ibond,2)
          dE = CGBond(ibond,4) * db       
          E = 0.5d0 * dE * db
        else 
           db = x - CGBond(ibond,2)
           dE1 = CGBond(ibond,4) * db 
           E1 = 0.5d0 * dE1 * db + CGBond(ibond,3)

           db = x - CGBond(ibond,5)
           dE2 = CGBond(ibond,7) * db
           E2 = 0.5d0 * dE2 * db + CGBond(ibond,6)
          
           E = Mix(E1,E2,dG1,dG2,CGBond(ibond,8))
           dE = dG1 * dE1 + dG2 * dE2
          if ( CGBond(ibond,1) .eq. 3) then
             E12 = E
             dE12 = dE

             db =  x - CGBond(ibond,9)
             dE3 = CGBond(ibond,11) * db
             E3 = 0.5d0 * dE3 * db + CGBond(ibond,10)

             E = Mix(E12,E3,dG1,dG2,CGBond(ibond,12))
             dE = dG1 * dE12 + dG2 * dE3
          endif 
        endif
C         write(6,*)"Bond: r, Ene, dE, req1, req2, req3, type:",
C     &               x, E, dE, CGBond(ibond,2),CGBond(ibond,5),
C     & CGBond(ibond,9), ibond
        if (E .gt. E_CG_max) then
        write(6,*)"Bond: E, dE, r, req1, req2, req3:"
        write(6,*)E, dE,x,CGBond(ibond,2),CGBond(ibond,5)
     &            ,CGBond(ibond,9)
        write(6,*)"Type: ",ibond
        endif
        
      return
      end

      function Mix(U1,U2,dM_dU1,dM_dU2,beta)
        implicit none

        double precision Mix, U1,U2, dM_dU1, dM_dU2, beta
        double precision s, s3

          s = 1.d0 / dsqrt( (U1-U2)**2 + beta)
          s3 = s**3

          dM_dU1 = 0.5d0 * (1.d0 - (U1 - U2)*s)
          dM_dU2 = 0.5d0 * (1.d0 + (U1 - U2)*s)
          Mix = 0.5d0 * (U1 + U2 - dsqrt( (U1-U2)**2 + beta ) )
C         write(6,*)"Mix: U1, U2, beta, U",U1, U2, beta, Mix
          return 
      end
