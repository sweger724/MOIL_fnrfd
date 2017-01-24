      subroutine eCG_WellAnglePot(x,iang,E,dE)
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c calculate bond energy and forces for i-th bond.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision x,db,dE,E, dG1, dG2
      integer iang
      double precision Mix, dE1, dE2, dE3, dE12, E1, E2, E3, E12

c  DATASTRUCTURE CGAngle
c  column   1       2      3      4      5        6        7    8        9      10      11   12
c  data    type   r_eq1  e_min1   k_1   r_eq2   e_min_2   k_2  beta1   r_eq3  e_min3   k_3   beta2
c
c  type encodes how many wells are there (1,2,3)
c
c   for type==1     E = 1/2 k_1 (r - r_eq1)^2
c
c   for type==2     E_1 = 1/2 k_1 (r - r_eq1)^2
c                   E_2 = 1/2 k_2 (r - r_eq2)^2
c                   E = 1/2 * ( E_1 + E_2 - sqrt((E_1-E_2)^2 + beta_1) )
c
c   for type==3    E_12 = E (from type==2)
c                  E_3 = 1/2 k_3 (r - r_eq3)^2
c                  E = 1/2 * ( E_12 + E_3 - sqrt((E_12-E_3)^2 + beta_2) )


        if (CGAngle(iang,1).eq.1) then
          db = x - CGAngle(iang,2)
          dE = CGAngle(iang,4) * db       
          E = 0.5d0 * dE * db
        else
C         dE = 0.d0
C         E =0.d0
C         return
          db = x - CGAngle(iang,2)
          dE1 = CGAngle(iang,4) * db 
           E1 = 0.5d0 * dE1 * db + CGAngle(iang,3)

           db = x - CGAngle(iang,5)
          dE2 = CGAngle(iang,7) * db
           E2 = 0.5d0 * dE2 * db + CGAngle(iang,6)
          
           E = Mix(E1,E2,dG1,dG2,CGAngle(iang,8))
          dE = dG1 * dE1 + dG2 * dE2
          if ( CGAngle(iang,1) .eq. 3) then
             E12 = E
             dE12 = dE

             db =  x - CGAngle(iang,9)
             dE3 = CGAngle(iang,11) * db
             E3 = 0.5d0 * dE3 * db + CGAngle(iang,10)

             E = Mix(E12,E3,dG1,dG2,CGAngle(iang,12))
             dE = dG1 * dE12 + dG2 * dE3
          endif 
        endif

        if (E .gt. E_CG_max) then
        write(6,*)"Angle: r, Ene, dE, req1, req2, req3, type:",
     &               x, E, dE, CGAngle(iang,2),
     &               CGAngle(iang,5),CGAngle(iang,9), iang
        endif
        
      return
      end
