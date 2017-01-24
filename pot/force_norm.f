        subroutine force_norm(gradf)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        double precision gradf
c local
        integer i
        double precision tmp

        gradf = 0.d0

C beginning of the old loop

C        do 1 j=1,npt
C        i = prtc_pointer(j)
        
        do 1 i=1,npt
         gradf = gradf + dpot(1,i)*dpot(1,i) + dpot(2,i)*dpot(2,i)
     1           + dpot(3,i)*dpot(3,i)
1       continue


        gradf = dsqrt(gradf/(3*npt))

        if (gradf.gt.1.d6) then
         do 2 i=1,npt
         tmp =  dpot(1,i)*dpot(1,i) + dpot(2,i)*dpot(2,i)
     1           + dpot(3,i)*dpot(3,i)
          if (tmp.gt.1.d3) then
                write(6,*)
                write(6,*)' ---------------------'
                write(6,*)' norm force  = ',tmp, ' for  atom ',i
                write(6,*)' ---------------------'
                write(6,*)
          end if
2       continue
        end if

        return 
        end

