        subroutine ener14()


        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
c mauro -muta
        include 'COMMON/MUTA.BLOCK'


        double precision rx,ry,rz,r2,s2,a,b,e1,e2,df,df1,df2
        double precision s,s6
        double precision pick, epstmp

        integer i,j,k


        e_vdw14 = 0.0d0
        e_el14  = 0.0d0

        if (eelyes) then
                epstmp = 1.0d0
        else
                epstmp =  0.0d0
        end if

        if (evdyes) then
                pick = 1.0d0
        else
                pick = 0.0d0
        end if


                 do 100 i=1,totspe
                        j = spec1(i)
                        k = spec2(i)
                        rx = coor(1,j) - coor(1,k)
                        ry = coor(2,j) - coor(2,k)
                        rz = coor(3,j) - coor(3,k)
                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2

                        if (p14(1,i).gt.1.d-12) then
                         s6=s2*s2*s2
                         a = p14(1,i)*s6*s6*pick
                         b = p14(2,i)*s6*pick
                         e1 = (a - b)
                         df1 = -6.0d0*s2*(a+e1)
                        else
                         e1  = 0.d0
                         df1 = 0.d0
                        end if

                        if (dabs(p14(3,i)).gt.1.d-12) then
                         s = dsqrt(s2)
                         e2 = p14(3,i)*s*epstmp
                         df2 = -e2*s2
                        else
                         e2  = 0.d0
                         df2 = 0.d0
                        end if
                         
                        df = df1 + df2
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        dpot(1,j) = dpot(1,j) + rx
                        dpot(2,j) = dpot(2,j) + ry
                        dpot(3,j) = dpot(3,j) + rz
                        dpot(1,k) = dpot(1,k) - rx
                        dpot(2,k) = dpot(2,k) - ry
                        dpot(3,k) = dpot(3,k) - rz
                        e_vdw14 = e_vdw14 + e1 
                        e_el14 = e_el14 + e2
c mauro - muta <U2-U1>l
                        if (mutaid(j).eq.1 .or. mutaid(k).eq.1) then
                         tmp_e_lambda=tmp_e_lambda+e1/lambda
                         tmp_e_lambda=tmp_e_lambda+e2/lambda
                   else if (mutaid(j).eq.2 .or. mutaid(k).eq.2) then
                         tmp_e_lambda=tmp_e_lambda-e1/(1.0d0-lambda)
                         tmp_e_lambda=tmp_e_lambda-e2/(1.0d0-lambda)
                        endif

               !write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~'
               !write(stdo,*)'i and j are ',j,k,dsqrt(r2)
               !write(stdo,'(a,f18.10)')'e1--vdW is ',e1
               !write(stdo,'(a,f8.4)')'e2--elec is ',e2
               !write(stdo,*)'epsgm12 of i and j are',epsgm12(j),epsgm12(k)
               !write(stdo,*)'epsgm6 of i and j are',epsgm6(j),epsgm6(k)
               !write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~~' 
100              continue


        return 
        end
