        subroutine nbfinit()

C       cdie 

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'

        double precision rx,ry,rz,r2,e1,e0
        double precision xi,yi,zi,dxi,dyi,dzi

        integer i,j,k,jbeg1,jend1,jbeg2,jend2

        save e0

        e_vdw = 0.0d0
        e_el=0.d0


        jbeg1=1
        jbeg2=1

        do 400 i=1,npt-1


                jend1 = point1(i)
                xi  = coor(1,i)
                yi  = coor(2,i)
                zi  = coor(3,i)
                dxi = 0.d0
                dyi = 0.d0
                dzi = 0.d0
                e0  = 1.d6*ptwei(i)

                if (jbeg1.lt.jend1) then

                do 200 k=jbeg1,jend1
                        j=list1(k)
                        rx = xi - coor(1,j)
                        ry = yi - coor(2,j)
                        rz = zi - coor(3,j)
                        r2=rx*rx+ry*ry+rz*rz
                        e1 = e0*exp(-40.d0*r2)
                        if (lesid(i).ne.lesid(j)) e1 = e1*ptwei(j)

                        rx = -2.d0*e1*rx*40.d0
                        ry = -2.d0*e1*ry*40.d0
                        rz = -2.d0*e1*rz*40.d0
                        dxi = dxi + rx
                        dyi = dyi + ry
                        dzi = dzi + rz
                        dpot(1,j) = dpot(1,j) - rx
                        dpot(2,j) = dpot(2,j) - ry
                        dpot(3,j) = dpot(3,j) - rz
                        e_vdw = e_vdw + e1 
200             continue
                end if
                jbeg1 = jend1 + 1
  
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c


              jend2 = point2(i)


              if (jbeg2.le.jend2) then
              do 205 k=jbeg2,jend2
                      j=list2(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                        e1 = e0*exp(-r2)
                        if (lesid(i).ne.lesid(j)) e1 = e1/ptwei(j)

                        rx = -2.d0*e1*rx
                        ry = -2.d0*e1*ry
                        rz = -2.d0*e1*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                        dpot(1,j) = dpot(1,j) - rx
                        dpot(2,j) = dpot(2,j) - ry
                        dpot(3,j) = dpot(3,j) - rz

                      e_vdw = e_vdw + e1 
205           continue
              end if
              jbeg2 = jend2 + 1



                dpot(1,i) = dpot(1,i) + dxi
                dpot(2,i) = dpot(2,i) + dyi
                dpot(3,i) = dpot(3,i) + dzi

400     continue

        return 
        end
