        subroutine stdc(too_small,fsave)

c       A subroutine that preforme several steepest decent iterations   

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        logical too_small
        double precision fsave(3,maxpt)

c local
        double precision rsave(3,maxpt)
        integer i,j,l
        integer nsteps

        double precision e_ref
        double precision step_size,step_size0,very_small
        data step_size,step_size0,very_small/0.001d0,0.1d0,1.d-20/
        data nsteps/20/

        do 100 i = 1,nsteps
                call eforce()
                do 1 j = 1,inofrz
                        l = nofreez(j)
                        rsave(1,l) = coor(1,l)
                        rsave(2,l) = coor(2,l)
                        rsave(3,l) = coor(3,l)
                        fsave(1,l) = dpot(1,l)
                        fsave(2,l) = dpot(2,l)
                        fsave(3,l) = dpot(3,l)
1               continue
                e_ref = e_total
                step_size = step_size * 1.5d0
                if (step_size .gt. step_size0) 
     *                          step_size = step_size0
2               continue
                do 3 j = 1,inofrz
                        l = nofreez(j)
                        coor(1,l) = rsave(1,l) - fsave(1,l) * step_size
                        coor(2,l) = rsave(2,l) - fsave(2,l) * step_size
                        coor(3,l) = rsave(3,l) - fsave(3,l) * step_size
3               continue
                call eforce()
                if (e_ref .gt. e_total) then
                   step_size = step_size*0.5d0
                   if (step_size .lt. very_small) then
                    write(stdo,*)'**************************'
                    write(stdo,*)'    STEP TOO SMALL'
                    write(stdo,*)'**************************'
                    too_small = .true.
                    return
                   end if
                   go to 2
                end if
100     continue
        too_small = .false.
        return
        end

