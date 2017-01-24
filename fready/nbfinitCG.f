        ! this subroutine adopt finite repulsion between all
        ! particles in order to avoid sterical clashes... 
        
        ! used typicaly for initial correcting of "bad" structures
        subroutine nbfinitCG()

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision rx,ry,rz,r2,e1,e0
        double precision xi,yi,zi,dxi,dyi,dzi

        integer i,j,k

        save e0

        e_vdw = 0.0d0
        e0  = 0.5d3

        do i=1,npt-1

                xi  = coor(1,i)
                yi  = coor(2,i)
                zi  = coor(3,i)
                dxi = 0.d0
                dyi = 0.d0
                dzi = 0.d0

                do k=LJ_list1(i-1)+1,LJ_list1(i)
                        j=LJ_list2(k)
                        rx = xi - coor(1,j)
                        ry = yi - coor(2,j)
                        rz = zi - coor(3,j)
                        
                        r2=rx*rx+ry*ry+rz*rz
                        e1 = e0*exp(-0.3d0*r2)

                        rx = -2.d0*e1*rx*0.3d0
                        ry = -2.d0*e1*ry*0.3d0
                        rz = -2.d0*e1*rz*0.3d0
                        dxi = dxi + rx
                        dyi = dyi + ry
                        dzi = dzi + rz
                        dpot(1,j) = dpot(1,j) - rx
                        dpot(2,j) = dpot(2,j) - ry
                        dpot(3,j) = dpot(3,j) - rz
                        e_vdw = e_vdw + e1 
C                        if(r2.lt.16) write(6,*)"repuls:",sqrt(r2),e1
                end do
  
                dpot(1,i) = dpot(1,i) + dxi
                dpot(2,i) = dpot(2,i) + dyi
                dpot(3,i) = dpot(3,i) + dzi

        end do

        return 
        end
