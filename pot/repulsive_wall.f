        subroutine repulsive_wall()
c
c calculating repulsive "wall" as flat exponential repulsion from +/- B/
c along the wnorm direction
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/REPWALL.BLOCK'

        integer i,j
        double precision norm,d12,d6,d1
c
c Disperssion replusion part
c
        e_rwall = 0.d0
        do 2 i=1,npt
        do j=1,nwalls
                 d1   = coor(normw(j),i)-w0(j)
C <------------0------------>
C  |   0->                 |
C positive force d>0
C  |                 <-0   |
C negative force d<0 

                 d2 = d1*d1
                 d6 = d2*d2*d2
                 d12 = d6*d6
                e_rwall    = e_rwall + weps/d12
           dpot(normw(j),i) = dpot(normw(j),i) - 12.d0*weps/(d12*d1)
        enddo
2       continue

        return
        end

