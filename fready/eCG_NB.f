        Subroutine eCG_NB(r,T,E,df)

C       Lennard Jones CG energy terms
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision rx,ry,rz,r2,s,s2,s3,s4,s6,s10,s12,a,b,E,df,r
        double precision xi,yi,zi,dxi,dyi,dzi, x(0:9)
        integer k,i,j,jbeg,jend, T, loop

        E = 0.d0
        df = 0.d0

        if (r .gt. 13.5d0 ) return

        r2 = r**2
        s2 = 1.0d0/r2
        s3 = s2/r
        s4 = s2**2
        s6 = s4*s2
        s12 = s6**2
        s10 = s6*s4

C                 write(6,*) "r: ",r, " Type: ", T, " cutt: ",LJa(T,10)
                  if (r .gt. LJa(T,10)) then
                    x(0) = 1.d0
                    x(1) = r - LJa(T,11)
                    do loop = 2,9
                      x(loop) = x(loop-1)*x(1)
                    end do
                    E = LJa(T,0)
                    df = 0.d0
                    do loop = 1,9
                      E = E + x(loop) * LJa(T,loop)
                      df = df + x(loop-1) * LJa(T,loop) * loop
                    end do
                    df = df/r
C                   write(6,*)"Type r, E1, df: ",T,r,E,df
C                   write(6,*)"LJa:",LJa(T,1),LJa(T,8),LJa(T,9)
                  else

C  6_2 potential
                    E = 1.d0*(LJr(T,1) * s6  + LJr(T,2) * s2 + LJr(T,3))
                    df =1.d0*(-6.d0*LJr(T,1)*s6*s2 -2.d0*LJr(T,2)*s4)

C                   write(6,*)"LJr:",LJr(T,1),LJr(T,2),LJr(T,3)
C                   write(6,*)"Type, r, E, df: ",T,r,E,df
                  endif
C Add default force
             if (T .le. Namino**2) then
                  E = E  + LJr(T,4) * s12  + LJr(T,5) * s6
                 df = df - 12.d0*LJr(T,4)*s12*s2 -6.d0*LJr(T,5)*s6*s2
             else 
                  E =  E + LJr(T,4) * s12
                  df = df  - 12.d0*LJr(T,4)*s12*s2
             endif


            if (smooth_hardcore .and. E .gt. 0.6) then 
               df = df / (10*(E+0.4))
               E = (0.6 + log(E + 0.4)/10)
            end if

C                if (E .gt. E_CG_max) then
C                   write(6,*)"Type, r, E, df: ",T,r,E,df
C                endif

        return 
        end






