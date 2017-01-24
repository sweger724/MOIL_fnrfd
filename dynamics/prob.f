          subroutine probability(n)
          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/CONNECT.BLOCK'
          include 'COMMON/VELOC.BLOCK'
          include 'COMMON/SPECL.BLOCK'
          include 'COMMON/ENERGY.BLOCK'
          include 'COMMON/CONVERT.BLOCK'
          include 'COMMON/EWALD.BLOCK'

          double precision  delta2,gama
          double precision  tempvx,tempvy,tempvz
          double precision  tempv,v,h
          double precision  tempe,temp,P
          real   radom 
          integer i,k,namel,level,n
          character*11 name
          data name/'probability'/
          data namel/11/

c F1  -  force from the repulsion function at crossing point
c F2  -  force from the morse potential function at crossing point
c Force   -  Force = DABS(F2-F1)

c delta - splitting between the rpulsion and morse potential curves
c delta2 - square of delta   

c  v    - velocity  ( A/s )
c          unit of tempv is  A/time   here   time*tconv = ps
c  h    -  Planck constant (kcal*s/mol) 

         if(Force.le.0. ) then
            level=1
            call alert(name,namel,'Force in Landau-zener
     $ formula "Forc" is 0',16,level)
          end if
          if (delt .le. 0. ) then 
            level = 1
            call alert(name,namel,'delt is 0.',16,level)
          end if

          delta2 = delt * delt
           i = imb1(n)
           k = imb2(n)
           tempvx = velo(1,k) - velo(1,i)
           tempvy = velo(2,k) - velo(2,i)
           tempvz = velo(3,k) - velo(3,i)
           tempv  = dsqrt( tempvx*tempvx + tempvy*tempvy 
     $             + tempvz*tempvz )
           v  = (tempv * 1.d12) / tconv
           
           h  = 6.626176d-37 * 6.022d23 /4.18d0

           temp = h * v * Force
           gama = delta2 * 2.d0 * 3.1415d0 / temp
           
           tempe = DEXP(-3.1415d0 * gama / 2d0)
           P = 1 - tempe
          
           call RANLUX(radom,1)

           if ( radom .le. P ) then
            write(*,*)
            write(*,*)'--------------------'
            write(*,*)' A switch is made '
            if (repyes(n)) then
             repyes(n) = .false.
             emyes(n)  = .true.
             write(*,1000)n,P,radom
1000         format(1x,' n = ',i7,' P = ',f7.5,' random = ',f7.5)
             write(*,*)' Morse is used '
             else
              repyes(n) = .true.
              emyes(n)  = .false.
             write(*,1000)n,P,radom
             write(*,*)' Repulsion is used '
            end if
            write(*,*)'--------------------'
            write(*,*)
            call eforce()
           end if
          return
          end
