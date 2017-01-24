      subroutine angle_constraint()

      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'

        character*8 name
        integer namel,i,j,l,k,k1,k2,i1,i2,i3,i4, rr

        double precision D0, D1, D2, K_u, dd, ff, K_2
        double precision dF(3,maxpt)
        double precision upper,lower,neigh, myT, torsion
        double precision pi2,aa,aa_min,aa1,aa2,aa3, DistAngle
        double precision bb,bb_min,bb1,bb2,bb3
        logical first
        data first/.true./

        save first,lower,upper,neigh, K_u, K_2

        name = 'umbrella'
        namel= 8

        pi2 = 3.14159265 * 2.d0

        if (first) then
          lower = dx - dx_umbrella
          upper = dx
          neigh = dx 
          first = .false.
         !K_u = 1000.d0 ! 100  ! for 18 cells
         !K_u = 100.d0 ! for wh5
         !K_2 = K_u
          K_u = umbrella_K1
          K_2 = umbrella_K2
        end if

        k1 = (myCell-1)*npt
        k2 = (myCell2-1)*npt

        e_Kumb = 0.d0
        do i = 1,npt
          do l=1,3
            dF(l,i) = 0.d0
          end do
        end do
        call GetReducedDerivs()

        D1 = DistAngle( ReducedCoor(1,myCell), myReduced(1), Nreduced)
        D2 = DistAngle( ReducedCoor(1,myCell2), myReduced(1),Nreduced)
        
        do i = 1, Ncells
          if (cell(i) .ne. myCell .and. cell(i).ne. myCell2) then
            
            D0 = DistAngle(ReducedCoor(1,cell(i)),myReduced(1),Nreduced)
            
            if ( D1 - D0 .gt. neigh) then
              e_Kumb = e_Kumb + 0.5d0 * K_2 * (D1 - D0 - neigh)**2
              
              ff = 2.d0 * K_2 * (D1 -D0 - neigh)

              do rr = 1, Nreduced
                aa1 = ReducedCoor(rr,cell(i)) - myReduced(rr)
                aa2 = ReducedCoor(rr,cell(i))-myReduced(rr)-pi2
                aa3 = ReducedCoor(rr,cell(i))-myReduced(rr)+pi2
                aa_min = min(abs(aa1),abs(aa2),abs(aa3)) 
                if (abs(aa1).eq.aa_min) aa = aa1
                if (abs(aa2).eq.aa_min) aa = aa2
                if (abs(aa3).eq.aa_min) aa = aa3
                aa = aa * tor_weight2(rr)


                bb1 = ReducedCoor(rr,myCell) - myReduced(rr)
                bb2 = ReducedCoor(rr,myCell)-myReduced(rr)-pi2
                bb3 = ReducedCoor(rr,myCell)-myReduced(rr)+pi2
                bb_min = min(abs(bb1),abs(bb2),abs(bb3))
                if (abs(bb1).eq.bb_min) bb = bb1
                if (abs(bb2).eq.bb_min) bb = bb2
                if (abs(bb3).eq.bb_min) bb = bb3
                bb = bb * tor_weight2(rr)
              
                do j = 1,4
                  do l = 1,3
                    dd = ff * myDReduced(l,j,rr) * (aa - bb)
                    dF(l,torsions(rr,j)) = dF(l,torsions(rr,j)) + dd
                  end do
                end do
              end do
            end if
            
          end if
        end do  

        if (D1 - D2 .gt. upper) then
          e_Kumb = e_Kumb + 0.5d0*K_u * (D1 - D2 -upper)**2
     
          ff = 2.d0 * K_u * (D1-D2 - upper)

          do rr = 1, Nreduced
            aa1 = ReducedCoor(rr,myCell2) - myReduced(rr)
            aa2 = ReducedCoor(rr,myCell2) - myReduced(rr) - pi2
            aa3 = ReducedCoor(rr,myCell2) - myReduced(rr) + pi2
            aa_min = min(abs(aa1),abs(aa2),abs(aa3))
            if (abs(aa1).eq.aa_min) aa = aa1
            if (abs(aa2).eq.aa_min) aa = aa2
            if (abs(aa3).eq.aa_min) aa = aa3
            aa = aa * tor_weight2(rr)
 
            bb1 = ReducedCoor(rr,myCell) - myReduced(rr)
            bb2 = ReducedCoor(rr,myCell) - myReduced(rr) - pi2
            bb3 = ReducedCoor(rr,myCell) - myReduced(rr) + pi2
            bb_min = min(abs(bb1),abs(bb2),abs(bb3))
            if (abs(bb1).eq.bb_min) bb = bb1
            if (abs(bb2).eq.bb_min) bb = bb2
            if (abs(bb3).eq.bb_min) bb = bb3
            bb = bb * tor_weight2(rr)
 
            do j = 1, 4
              do l = 1,3
                dd = ff * myDReduced(l,j,rr) * (aa - bb)
                dF(l,torsions(rr,j)) = dF(l,torsions(rr,j)) + dd
              end do 
            end do
            
          end do ! rr
        end if

        if (D1-D2 .lt. lower) then
          e_Kumb = e_Kumb + 0.5d0 * K_u * (D1-D2-lower)**2

          ff = 2.d0 * K_u * (D1-D2-lower)

          do rr = 1, Nreduced
            aa1 = ReducedCoor(rr,myCell2) - myReduced(rr)
            aa2 = ReducedCoor(rr,myCell2) - myReduced(rr) - pi2
            aa3 = ReducedCoor(rr,myCell2) - myReduced(rr) + pi2
            aa_min = min(abs(aa1),abs(aa2),abs(aa3))
            if (abs(aa1).eq.aa_min) aa = aa1
            if (abs(aa2).eq.aa_min) aa = aa2
            if (abs(aa3).eq.aa_min) aa = aa3
            aa = aa * tor_weight2(rr)

            bb1 = ReducedCoor(rr,myCell) - myReduced(rr)
            bb2 = ReducedCoor(rr,myCell) - myReduced(rr) - pi2
            bb3 = ReducedCoor(rr,myCell) - myReduced(rr) + pi2
            bb_min = min(abs(bb1),abs(bb2),abs(bb3))
            if (abs(bb1).eq.bb_min) bb = bb1
            if (abs(bb2).eq.bb_min) bb = bb2
            if (abs(bb3).eq.bb_min) bb = bb3
            bb = bb * tor_weight2(rr)

c
c myDReduced is the derivative of atom j of  torsion rr with corresponding coordinate
c j (j=x,y,z) 
c
            do j = 1, 4
              do l = 1,3
                dd = ff * myDReduced(l,j,rr) * (aa - bb)
                dF(l,torsions(rr,j)) = dF(l,torsions(rr,j)) + dd
              end do
            end do            
          end do ! rr
        end if 

        do i=1,npt
             dpot(1,i) = dpot(1,i) + dF(1,i)
             dpot(2,i) = dpot(2,i) + dF(2,i)
             dpot(3,i) = dpot(3,i) + dF(3,i)
        end do

        return
      end
