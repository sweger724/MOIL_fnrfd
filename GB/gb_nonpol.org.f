c gb_nonpol.f
c    computing non-polar solvation energy in GB model
c
c    YS, Mar 10, 2005
c
c refs: Weiser et al JCC 2,217(1999)

 
      SUBROUTINE init_gb_nonpol
c     Initialization of parameters for nonpolar solvation
c     energy calculation.

      IMPLICIT none
      INCLUDE 'COMMON/LENGTH.BLOCK'
      INCLUDE 'COMMON/SGB.BLOCK'
      INCLUDE 'COMMON/CONNECT.BLOCK'
      INCLUDE 'COMMON/PROPERT.BLOCK'
    
      logical SASA_NO_NLR

c begin
      write(*,*) "YS GB_NONPOL> initial gb_nonpol  ..." 
      SASA_NO_NLR = .TRUE. 

      if ( SASA_NO_NLR ) then
c --------------------------------------------------------
c 22 particle types for LCPO surface: 
c refs: Weiser et al JCC 2,217(1999), Appendix B


c hydrogens  said = 1
          P0(1) = 1.4  
          P1(1) = 0.0
          P2(1) = 0.0
          P3(1) = 0.0
          P4(1) = 0.0

c sp3     aliphatic C, # bonded neighbor = 1, said = 2
          P0(2) =  1.70 + 1.40
          P1(2) =  0.7887 
          P2(2) = -0.28063 
          P3(2) = -0.0012968 
          P4(2) =  0.00039328  

c sp3     aliphatic C, # bonded neighbor = 2, said = 3
          P0(3) =  1.70 + 1.40
          P1(3) =  0.56482 
          P2(3) = -0.19608 
          P3(3) = -0.0010219 
          P4(3) =  0.0002658  

c sp3     aliphatic C, # bonded neighbor = 3, said = 4
          P0(4) =  1.70 + 1.40
          P1(4) =  0.23348 
          P2(4) = -0.072627 
          P3(4) = -0.00020079 
          P4(4) =  0.00007967  

c sp3     aliphatic C, # bonded neighbor = 4, said = 5
          P0(5) =  1.70 + 1.40
          P1(5) =  0.0 
          P2(5) =  0.0 
          P3(5) =  0.0 
          P4(5) =  0.0  

c sp2      C, # bonded neighbor = 2, said = 6
          P0(6) =  1.70 + 1.40
          P1(6) =  0.51245 
          P2(6) = -0.15966  
          P3(6) = -0.00019781 
          P4(6) =  0.00016392  

c sp2      C, # bonded neighbor = 3, said = 7
          P0(7) =  1.70 + 1.40
          P1(7) =  0.070344 
          P2(7) = -0.019015  
          P3(7) = -0.000022009 
          P4(7) =  0.000016875  

c oxyg    en in hydroxyl,sp3, bonded neighbor 1, said = 8
          P0(8) =  1.60 + 1.40
          P1(8) =  0.77914 
          P2(8) = -0.25262  
          P3(8) = -0.0016056 
          P4(8) =  0.00035071  

c oxyg    en ,sp3, bonded neighbor 2, said = 9 
          P0(9) =  1.60 + 1.40
          P1(9) =  0.49392 
          P2(9) = -0.16038  
          P3(9) = -0.00015512 
          P4(9) =  0.00016453  

c carb    onyl group oxygen sp2 bonded neighbour 1, said = 10 
          P0(10) =  1.60 + 1.40
          P1(10) =  0.68563 
          P2(10) = -0.1868  
          P3(10) = -0.00135573 
          P4(10) =  0.00023743  

c carboxyl and phosphate group oxygen, said = 11
          P0(11) =  1.60 + 1.40
          P1(11) =  0.88857 
          P2(11) = -0.33421  
          P3(11) = -0.0018683 
          P4(11) =  0.00049372  

c sp3     N, bonded neighbor 1, said = 12
          P0(12) =  1.65 + 1.40
          P1(12) =  0.78602 
          P2(12) = -0.29198  
          P3(12) = -0.0006537 
          P4(12) =  0.00036247  

c sp3     N, bonded neighbor 2, said = 13
          P0(13) =  1.65 + 1.40
          P1(13) =  0.22599 
          P2(13) = -0.036648  
          P3(13) = -0.0012297 
          P4(13) =  0.000080038  

c sp3     N, bonded neighbor 3, said = 14
          P0(14) =  1.65 + 1.4
          P1(14) =  0.051481
          P2(14) = -0.012603
          P3(14) = -0.00032006
          P4(14) =  0.000024774

c sp2     N, bonded neighbor, 1, said = 15
          P0(15) =  1.65 + 1.4
          P1(15) =  0.73511
          P2(15) = -0.22116
          P3(15) = -0.00089148
          P4(15) =  0.0002523

c sp2     N, bonded neighbor, 2, said = 16
          P0(16) =  1.65 + 1.4
          P1(16) =  0.41102
          P2(16) = -0.12254
          P3(16) = -0.000075448
          P4(16) =  0.00011804

c sp2     N, bonded neighbor 3, said = 17
          P0(17) =  1.65 + 1.4
          P1(17) =  0.062577
          P2(17) = -0.017874
          P3(17) = -0.00008312
          P4(17) =  0.000019849

c sulp    hur SH, bonded neighbor 1, said = 18
          P0(18) =  1.9 + 1.4
          P1(18) =  0.7722
          P2(18) = -0.26393
          P3(18) = -0.0010629
          P4(18) =  0.0002179

c sulp    hur S, bonded neighbor 2, said = 19
          P0(19) =  1.9 + 1.4
          P1(19) =  0.54581
          P2(19) = -0.19477
          P3(19) = -0.0012873
          P4(19) =  0.00029247

c P, b    onded neighbor 3, said = 20
          P0(20) =  1.9 + 1.4
          P1(20) =  0.3865
          P2(20) = -0.18249
          P3(20) = -0.0036598
          P4(20) =  0.0004264

c P, b    onded neighbor 4, said = 21
          P0(21) =  1.9 + 1.4
          P1(21) =  0.03873
          P2(21) = -0.0089339
          P3(21) =  0.0000083582
          P4(21) =  0.0000030381

c Cl,     said = 22
          P0(22) =  1.9 + 1.4
          P1(22) =  0.98318
          P2(22) = -0.40437
          P3(22) =  0.00011249
          P4(22) =  0.00040001
c -------------------------------------------------
c  using NLR to computing the SASA
C --------------------------------------------------
      else
c hydrogens  said = 1
          P0(1) = 1.4  
          P1(1) = 0.0
          P2(1) = 0.0
          P3(1) = 0.0
          P4(1) = 0.0

c sp3     aliphatic C, # bonded neighbor = 1, said = 2
          P0(2) =  1.70 + 1.40
          P1(2) =  0.86840 
          P2(2) = -0.41776 
          P3(2) = -0.0008575 
          P4(2) =  0.00078104  

c sp3     aliphatic C, # bonded neighbor = 2, said = 3
          P0(3) =  1.70 + 1.40
          P1(3) =  0.62286 
          P2(3) = -0.28190 
          P3(3) = -0.0024698 
          P4(3) =  0.00053606  

c sp3     aliphatic C, # bonded neighbor = 3, said = 4
          P0(4) =  1.70 + 1.40
          P1(4) =  0.28368 
          P2(4) = -0.12982 
          P3(4) = -0.0015757 
          P4(4) =  0.00024514  

c sp3     aliphatic C, # bonded neighbor = 4, said = 5
          P0(5) =  1.70 + 1.40
          P1(5) =  0.0 
          P2(5) =  0.0 
          P3(5) =  0.0 
          P4(5) =  0.0  

c sp2      C, # bonded neighbor = 2, said = 6
          P0(6) =  1.70 + 1.40
          P1(6) =  0.61006 
          P2(6) = -0.24859  
          P3(6) = -0.0019453 
          P4(6) =  0.00027405  

c sp2      C, # bonded neighbor = 3, said = 7
          P0(7) =  1.70 + 1.40
          P1(7) =  0.089938 
          P2(7) = -0.036938  
          P3(7) = -0.000022009 
          P4(7) =  0.000047525  

c stop here
c oxyg    en in hydroxyl,sp3, bonded neighbor 1, said = 8
          P0(8) =  1.60 + 1.40
          P1(8) =  0.77914 
          P2(8) = -0.25262  
          P3(8) = -0.0016056 
          P4(8) =  0.00035071  

c oxyg    en ,sp3, bonded neighbor 2, said = 9 
          P0(9) =  1.60 + 1.40
          P1(9) =  0.49392 
          P2(9) = -0.16038  
          P3(9) = -0.00015512 
          P4(9) =  0.00016453  

c carb    onyl group oxygen sp2 bonded neighbour 1, said = 10 
          P0(10) =  1.60 + 1.40
          P1(10) =  0.68563 
          P2(10) = -0.1868  
          P3(10) = -0.00135573 
          P4(10) =  0.00023743  

c carb    oxyl and phosphate group oxygen, said = 11
          P0(11) =  1.60 + 1.40
          P1(11) =  0.88857 
          P2(11) = -0.33421  
          P3(11) = -0.0018683 
          P4(11) =  0.00049372  

c sp3     N, bonded neighbor 1, said = 12
          P0(12) =  1.65 + 1.40
          P1(12) =  0.78602 
          P2(12) = -0.29198  
          P3(12) = -0.0006537 
          P4(12) =  0.00036247  

c sp3     N, bonded neighbor 2, said = 13
          P0(13) =  1.65 + 1.40
          P1(13) =  0.22599 
          P2(13) = -0.036648  
          P3(13) = -0.0012297 
          P4(13) =  0.000080038  

c sp3     N, bonded neighbor 3, said = 14
          P0(14) =  1.65 + 1.4
          P1(14) =  0.051481
          P2(14) = -0.012603
          P3(14) = -0.00032006
          P4(14) =  0.000024774

c sp2     N, bonded neighbor, 1, said = 15
          P0(15) =  1.65 + 1.4
          P1(15) =  0.73511
          P2(15) = -0.22116
          P3(15) = -0.00089148
          P4(15) =  0.0002523

c sp2     N, bonded neighbor, 2, said = 16
          P0(16) =  1.65 + 1.4
          P1(16) =  0.41102
          P2(16) = -0.12254
          P3(16) = -0.000075448
          P4(16) =  0.00011804

c sp2     N, bonded neighbor 3, said = 17
          P0(17) =  1.65 + 1.4
          P1(17) =  0.062577
          P2(17) = -0.017874
          P3(17) = -0.00008312
          P4(17) =  0.000019849

c sulp    hur SH, bonded neighbor 1, said = 18
          P0(18) =  1.9 + 1.4
          P1(18) =  0.7722
          P2(18) = -0.26393
          P3(18) = -0.0010629
          P4(18) =  0.0002179

c sulp    hur S, bonded neighbor 2, said = 19
          P0(19) =  1.9 + 1.4
          P1(19) =  0.54581
          P2(19) = -0.19477
          P3(19) = -0.0012873
          P4(19) =  0.00029247

c P, b    onded neighbor 3, said = 20
          P0(20) =  1.9 + 1.4
          P1(20) =  0.3865
          P2(20) = -0.18249
          P3(20) = -0.0036598
          P4(20) =  0.0004264

c P, b    onded neighbor 4, said = 21
          P0(21) =  1.9 + 1.4
          P1(21) =  0.03873
          P2(21) = -0.0089339
          P3(21) =  0.0000083582
          P4(21) =  0.0000030381

c Cl,     said = 22
          P0(22) =  1.9 + 1.4
          P1(22) =  0.98318
          P2(22) = -0.40437
          P3(22) =  0.00011249
          P4(22) =  0.00040001
      end if
c ------------------------------------------------------------------
      return
      END 

      SUBROUTINE calc_nonpol_energy(natom, e_gbnp)
c     Enonpol = surften* SASA
c     SASA based on LCPO model:  Weiser et al JCC, 20, 217(1999)

      IMPLICIT none
      INCLUDE 'COMMON/LENGTH.BLOCK'
      INCLUDE 'COMMON/SGB.BLOCK'
      INCLUDE 'COMMON/CONNECT.BLOCK'
      INCLUDE 'COMMON/COORD.BLOCK'
      INCLUDE 'COMMON/ENERGY.BLOCK'

c args type
      INTEGER natom
      DOUBLE PRECISION e_gbnp

c local variables
      DOUBLE PRECISION totsasa, r2, dij1i
      DOUBLE PRECISION xi, yi, zi, xj, yj, zj, xk, yk, zk
      INTEGER count, count2, icount, i,j,k, ip, jp, kp
      DOUBLE PRECISION si, sumAij, sumAjk, sumAijAjk, sumdAijddijdxi
      DOUBLE PRECISION sumdAijddijdyi,sumdAijddijdzi,sumdAijddijdxiAjk
      DOUBLE PRECISION sumdAijddijdyiAjk,sumdAijddijdziAjk, rij,
     &                 tmpaij, Aij, dAijddij
      DOUBLE PRECISION dAijddijdxj, dAijddijdyj, dAijddijdzj
      DOUBLE PRECISION sumdAjkddjkdxj,sumdAjkddjkdyj,sumdAjkddjkdzj,
     &                 p3p4Aij
      DOUBLE PRECISION rjk2, djk1i, rjk, vdw2dif, tmpajk,
     &                 Ajk, sumAjk2, dAjkddjk 
      DOUBLE PRECISION dAjkddjkdxj, dAjkddjkdyj, dAjkddjkdzj, 
     &                 lastxj, lastyj, lastzj
      DOUBLE PRECISION dAidxj, dAidyj, dAidzj, Ai, Aidxi, Aidyi,
     &                 Aidzi, dAidxi, dAidyi, dAidzi

      DOUBLE PRECISION PIx4, PIx2, PIx1
      PARAMETER (PIx4 = 12.5663706143591724639918,
     &           PIx2 =  6.2831853071795862319959,
     &           PIx1 =  3.1415926535897931159979)

c begin
      totsasa = 0.0D0
      write(*,*) "GB_NONPOL> computing SASA ..." 
      write(*,*) "GB_NONPOL> surften = ", surften 



c print out the neigbour list for debugging
c      count = 1 
c      do i = 1, natom
c          write(*,*) "NB_LIST> -----------------------------"
c          write(*,*) "NB_LIST> Neighours of atom ", i, " are:"
c501       if ( ineighbor(count) .ne. 0 ) then 
c             write (*, *) "NB_LIST> atom:", ineighbor(count)
c             count = count + 1
c          go to 501 
c          end if
c          count = count + 1
c      end do


      count = 1 
      do i = 1, natom
          ip = ptsaid(i)
c          write(*,*) "GB_NONPOL> i=", i, " ptsaid =", ptsaid(i) 
          if ( ineighbor(count) .eq. 0 ) then
              count = count + 1
          else
c obtaining Ai
              si=PIx4 * P0(ip)*P0(ip)
              sumAij = 0.0
              sumAjk = 0.0
              sumAjk2 = 0.0
              sumAijAjk = 0.0
              sumdAijddijdxi = 0.0
              sumdAijddijdyi = 0.0
              sumdAijddijdzi = 0.0
              sumdAijddijdxiAjk = 0.0
              sumdAijddijdyiAjk = 0.0
              sumdAijddijdziAjk = 0.0
              
              icount=count
770           j=ineighbor(count) 
              jp=ptsaid(j)
              xi = coor(1,i)
              yi = coor(2,i)
              zi = coor(3,i)
              xj = coor(1,j)
              yj = coor(2,j)
              zj = coor(3,j)
              r2=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj)
              dij1i = 1.0D0 / DSQRT(r2)
              rij = r2 * dij1i
              vdw2dif = P0(ip)*P0(ip) - 
     &                  P0(jp)*P0(jp)
              tmpaij = 2.0 * P0(ip) - rij - vdw2dif * dij1i
              Aij=PIx1 * P0(ip) * tmpaij
              dAijddij = PIx1*P0(ip)*(dij1i*dij1i*vdw2dif -1.0)         
              dAijddijdxj = dAijddij*(xj-xi)*dij1i 
              dAijddijdyj = dAijddij*(yj-yi)*dij1i 
              dAijddijdzj = dAijddij*(zj-zi)*dij1i 

              sumAij = sumAij + Aij

              count2 = icount
              sumAjk2 = 0.0
              sumdAjkddjkdxj = 0.0 
              sumdAjkddjkdyj = 0.0 
              sumdAjkddjkdzj = 0.0 

              p3p4Aij=-surften * (P3(ip) + P4(ip)*Aij)

780           k = ineighbor(count2) 
              kp = ptsaid(k)
              if ( j .eq. k ) go to 785
           
              xk = coor(1,k)
              yk = coor(2,k)
              zk = coor(3,k)
              rjk2 = (xj - xk) * (xj - xk) +
     &               (yj - yk) * (yj - yk) +
     &               (zj - zk) * (zj - zk) 
              djk1i = 1.0 / DSQRT(rjk2)
              rjk = rjk2 * djk1i
              if ( ( P0(jp) + P0(kp)) .gt. rjk ) then
                  vdw2dif = P0(jp) * P0(jp) - P0(kp) * P0(kp)
                  tmpajk = 2.0 * P0(jp) - rjk -vdw2dif * djk1i
                  Ajk = PIx1 * P0(jp) * tmpajk
                  sumAjk = sumAjk + Ajk
                  sumAjk2 = sumAjk2 + Ajk
                  dAjkddjk = PIx1*P0(jp)*(djk1i*djk1i*vdw2dif-1.0)
                  dAjkddjkdxj = dAjkddjk * (xj - xk ) * djk1i
                  dAjkddjkdyj = dAjkddjk * (yj - yk ) * djk1i
                  dAjkddjkdzj = dAjkddjk * (zj - zk ) * djk1i

c  add code here to update forces, ignore for now ....
                  dpot(1,k) = dpot(1,k) + dAjkddjkdxj*p3p4Aij
                  dpot(2,k) = dpot(2,k) + dAjkddjkdyj*p3p4Aij
                  dpot(3,k) = dpot(3,k) + dAjkddjkdzj*p3p4Aij

                  sumdAjkddjkdxj = sumdAjkddjkdxj + dAjkddjkdxj
                  sumdAjkddjkdyj = sumdAjkddjkdyj + dAjkddjkdyj
                  sumdAjkddjkdzj = sumdAjkddjkdzj + dAjkddjkdzj
     
              end if

785           count2 = count2 + 1
              if ( ineighbor(count2) .ne. 0 ) then
                  go to 780
              else
                  count2 = icount
              end if

              sumAijAjk = sumAijAjk + Aij * sumAjk2
              sumdAijddijdxi = sumdAijddijdxi - dAijddijdxj
              sumdAijddijdyi = sumdAijddijdyi - dAijddijdyj
              sumdAijddijdzi = sumdAijddijdzi - dAijddijdzj

              sumdAijddijdxiAjk = 
     &               sumdAijddijdxiAjk - dAijddijdxj * sumAjk2 
              sumdAijddijdyiAjk = 
     &               sumdAijddijdyiAjk - dAijddijdyj * sumAjk2 
              sumdAijddijdziAjk = 
     &               sumdAijddijdziAjk - dAijddijdzj * sumAjk2 

              lastxj = dAijddijdxj * sumAjk2 + Aij * sumdAjkddjkdxj
              lastyj = dAijddijdyj * sumAjk2 + Aij * sumdAjkddjkdyj
              lastzj = dAijddijdzj * sumAjk2 + Aij * sumdAjkddjkdzj
       
              dAidxj = surften * (P2(ip) * dAijddijdxj + 
     &                 P3(ip) * sumdAjkddjkdxj + P4(ip) * lastxj )
              dAidyj = surften * (P2(ip) * dAijddijdyj + 
     &                 P3(ip) * sumdAjkddjkdyj + P4(ip) * lastyj )
              dAidzj = surften * (P2(ip) * dAijddijdzj + 
     &                 P3(ip) * sumdAjkddjkdzj + P4(ip) * lastzj )

c here updating force array
              dpot(1,j) = dpot(1,j) + dAidxj
              dpot(2,j) = dpot(2,j) + dAidyj
              dpot(3,j) = dpot(3,j) + dAidzj

              count = count + 1
              if (ineighbor(count) .ne. 0 ) then
                  go to 770
              else
                  count = count + 1
              end if


              Ai = P1(ip) * si + P2(ip) * sumAij + P3(ip) * sumAjk +
     &             P4(ip) * sumAijAjk

c              write(*, *) "SASA> -------------------------"              
c              write(*, *) "SASA AREA> sasa of Atom", i , " is ", Ai 
c              write(*, *) "SASA> ip of Atom", i , " is ", ip 
c              write(*, *) "SASA> si=", si, 
c     &                    "    P1(ip)*si=",P1(ip)*si 
c              write(*, *) "SASA> sumAij=",sumAij,
c     &                    "    P2(ip)*sumAij=", P2(ip)*sumAij
c              write(*, *) "SASA> sumAjk=",sumAjk,
c     &                    "    P3(ip)*sumAjk=", P3(ip)*sumAjk
c              write(*, *) "SASA> sumAjkAjk=",sumAijAjk,
c     &                    "    P4(ip)*sumAijAjk=", P4(ip)*sumAijAjk
c            

              dAidxi = surften * (P2(ip) * sumdAijddijdxi +
     &                 P4(ip) * sumdAijddijdxiAjk ) 
              dAidyi = surften * (P2(ip) * sumdAijddijdyi +
     &                 P4(ip) * sumdAijddijdyiAjk ) 
              dAidzi = surften * (P2(ip) * sumdAijddijdzi +
     &                 P4(ip) * sumdAijddijdziAjk ) 

c here updating the force array, ignore for now
              dpot(1,i) = dpot(1,i) + dAidxi
              dpot(2,i) = dpot(2,i) + dAidyi
              dpot(3,i) = dpot(3,i) + dAidzi


c summing Ai
              totsasa = totsasa + Ai
c              write(*, *) "SASA TOT> totsasa up to ", i ,
c     &                    " is ", totsasa
          end if 
      end do

      write(*,*) "SASA> TOTSASA =  ", totsasa 
      e_gbnp = surften*totsasa
      write(*,*) "SASA> surften =  ", surften 
      write(*,*) "SASA> e_gbnp =  ", e_gbnp 
      write(*,*) "GB_NONPOL> nonpol energy calculated ..." 
      return
      end


c debug code to check derivatives

      SUBROUTINE check_derivative(natom)

      IMPLICIT none
      INCLUDE 'COMMON/LENGTH.BLOCK'
      INCLUDE 'COMMON/SGB.BLOCK'
      INCLUDE 'COMMON/CONNECT.BLOCK'
      INCLUDE 'COMMON/COORD.BLOCK'
      INCLUDE 'COMMON/ENERGY.BLOCK'

c args type
      INTEGER natom, i

c local variables
      DOUBLE PRECISION deltad,temp,deri, U1,U2, diff,sumdiff
      PARAMETER (deltad= 0.0000001D0)
c begin
      write (*,*) "natom=", natom
      sumdiff = 0D0
      do i = 1, natom
          write(*,*) "CHK_DERI> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"  
          write(*,*) "CHK_DERI> Atom:", i  

c calculate numerical dX
          temp = coor(1,i)
          coor(1,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(1,i) = temp - deltad
          call eforce()
          U2 = e_total
          coor(1,i) = temp
          deri =  ((U1-U2)/(2*deltad))
          diff = DABS((deri - dpot(1,i))/deri)
          write(*,*) "CHK_DERI> Numerical  dX: ", deri  
          write(*,*) "CHK_DERI> Analytical dX: ", dpot(1,i) 
          write(*,*) "CHK_DERI> diff         : ", diff 
          sumdiff = sumdiff + diff
      
c calculate numerical dY
          temp = coor(2,i)
          coor(2,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(2,i) = temp - deltad
          call eforce()
          U2 = e_total
          coor(2,i) = temp
          deri =  ((U1-U2)/(2*deltad))
          diff = DABS((deri - dpot(2,i))/deri)
          write(*,*) "CHK_DERI> Numerical  dY: ", deri  
          write(*,*) "CHK_DERI> Analytical dY: ", dpot(2,i) 
          write(*,*) "CHK_DERI> diff         : ", diff 
          sumdiff = sumdiff + diff

c calculate numerical dZ
          temp = coor(3,i)
          coor(3,i) = temp + deltad
          call eforce()
          U1 = e_total
          coor(3,i) = temp - deltad
          call eforce()
          U2 = e_total
          coor(3,i) = temp
          deri = ((U1-U2)/(2*deltad))
          diff = DABS((deri - dpot(3,i))/deri)
          write(*,*) "CHK_DERI> Numerical  dZ: ", deri  
          write(*,*) "CHK_DERI> Analytical dZ: ", dpot(3,i) 
          write(*,*) "CHK_DERI> diff         : ", diff 
          sumdiff = sumdiff + diff

      end do
      write(*,*) "CHK_DERI> ***********************************" 
      write(*,*) "CHK_DERI> Average diff : ", sumdiff/(3*natom) 
      
      return
      end
