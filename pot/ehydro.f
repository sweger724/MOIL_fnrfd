        subroutine ehydro()

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'


        integer l,i,j,k
        double precision les_scale,avg2
        double precision rx,ry,rz,r2,r,df,e,h,hydro_th2

        e_hyd   = 0.d0
        hydro_th2 = hydro_th * hydro_th
        avg2      = 2.d0*avg_hydro
        do 100 l=1,nbeta - 1
           i = betap(l)
           do 200 k = l + 1,nbeta
                les_scale = ptwei(i)
                j = betap(k)
                if ((lesid(i).ne.0) .and.
     *              (lesid(i) .eq.lesid(j)))  then
                    if (cplbl(i) .ne. cplbl(j)) go to 200
                end if
                if(lesid(i).ne.lesid(j))
     *                  les_scale = les_scale*ptwei(j)
                h = (cbeta(l)+cbeta(k)+avg2)*les_scale
                rx = coor(1,i) - coor(1,j)
                ry = coor(2,i) - coor(2,j)
                rz = coor(3,i) - coor(3,j)
                r2=rx*rx+ry*ry+rz*rz
                if (r2 .gt. hydro_th2) go to 200
                r=dsqrt(r2)
                e = hydro_scale*h*r
                df = hydro_scale*h/r
c debug
c               write(stdo,*) 'debug',i,j,r,e
                      rx = df*rx
                      ry = df*ry
                      rz = df*rz
                      dpot(1,i) = dpot(1,i) + rx
                      dpot(2,i) = dpot(2,i) + ry
                      dpot(3,i) = dpot(3,i) + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz
                e_hyd = e_hyd + e
200       continue
100     continue
        return
        end
