      subroutine interface_constraint()

      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        character*8 name
        integer namel,i,j,l,k,k1,k2

        double precision rms1, rms2, K_u, dd, ff, K_2
        double precision coor3(3,maxpt), dF(3,maxpt)
        double precision upper,lower,neigh
        logical first

        data first/.true./

        save first, lower, upper, neigh, K_u, K_2

        name = 'umbrella'
        namel= 8

        if (first) then
            lower = dx - dx_umbrella
            upper = dx
            neigh = dx
            first = .false.
            K_u = 0.3d4 * allmass/12.0
            K_2 = K_u*2.d0
            write(6,*) "KKK: ", K_u,K_2
            stop
        end if

        k1 = (myCell-1)*npt
        k2 = (myCell2-1)*npt

        rms1 = 0.d0
        rms2 = 0.d0

        e_Kumb = 0.d0
        !write(6,*) "KKK", K_u,K_2,ddx
        
        do i = 1,npt
          do l=1,3
            coor3(l,i) = coor(l,i)
            dF(l,i) = 0.d0
          end do
        end do


        ! align the coordinate to the centers frame
        do i = 1,Ncells+2
          l=(i-1)*npt+1
          call rmsd_weight(npt,coor3(1,1),centers(1,l),rms,.false.,mmm)
        end do

        if (my_pe .ne. num_pes -1) return

        do j = 1,npt
          do l =1,3
            rms1 = rms1 + mmm(j)*(coor3(l,j) - centers(l,k1+j))**2
            rms2 = rms2 + mmm(j)*(coor3(l,j) - centers(l,k2+j))**2
          end do
        end do

        rms1 = rms1/allmass
        rms2 = rms2/allmass
        !write(6,*)"rrr:", rms1, rms2,allmass

        do i = 1, Ncells
          if (cell(i) .ne. myCell .and. cell(i).ne. myCell2) then
            k = (cell(i)-1)*npt
            rms = 0.d0
            do j = 1,npt
              do l =1,3
                rms = rms + mmm(j)*(coor3(l,j) - centers(l,k+j))**2
              end do
            end do

            rms = rms/allmass
C            if (i .eq. 14) 
C             write(6,*)"rrr:",rms, rms1, rms2

            if (rms1 -rms .gt. neigh) then
              e_Kumb = e_Kumb + 0.5d0 * K_2 * (rms + neigh - rms1)**2
              
              ff = 2.d0 * K_2/allmass * (rms + neigh - rms1)
              
              do j = 1,npt
                do l = 1,3
                  dd = ff * (centers(l,k1+j) - centers(l,k+j)) *mmm(j)
                  dF(l,j) = dF(l,j) + dd
                end do
              end do
            end if
            
          end if
        end do  

        if (rms1 - rms2 .gt. upper) then
          e_Kumb = e_Kumb + 0.5d0*K_u * (rms1 - rms2 -upper)**2
          !write(6,*)"xxx",e_Kumb,rms1 + rms2, upper

          ff = 2.d0 * K_u/allmass * (rms1-rms2 - upper)
        !write(6,*) "RRR:",rms1,rms2,myCell,myCell2,rms1-rms2-upper


          do j = 1, npt
            do l = 1,3
              dd = ff * ( centers(l,k2+j) - centers(l,k1+j)) * mmm(j)
              dF(l,j) = dF(l,j) + dd
            end do 
          end do
        end if

        if (rms1-rms2 .lt. lower) then
          e_Kumb = e_Kumb + 0.5d0 * K_u * (rms1-rms2-lower)**2
          !write(6,*)"yyy",e_Kumb, rms1 + rms2,lower

          ff = 2.d0 * K_u/allmass * (rms1-rms2-lower)
        !write(6,*) "RRR:",rms1,rms2,myCell,myCell2,rms1-rms2-lower

          do j = 1, npt
            do l = 1,3
              dd = ff * ( centers(l,k2+j) - centers(l,k1+j) )* mmm(j)
              dF(l,j) = dF(l,j) + dd
            end do
          end do
        end if 

        do i=1,npt
             dpot(1,i) = dpot(1,i) + dF(1,i)
             dpot(2,i) = dpot(2,i) + dF(2,i)
             dpot(3,i) = dpot(3,i) + dF(3,i)
        end do
        !write(6,*)"EEE",e_Kumb

        return
      end
