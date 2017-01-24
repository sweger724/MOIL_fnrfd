      SUBROUTINE egb_nonpol2(natom, e_gbnp)
c Calculate the non-polarized solvation energy and 
c also the first derivatives.
c Enonpol = surften* SASA
c SASA based on LCPO model:  Weiser et al JCC, 20, 217(1999)

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'

c args type
      integer natom
      double precision e_gbnp

c local variables
      double precision totsasa, r2, dij1i
      double precision xi, yi, zi, xj, yj, zj, xk, yk, zk
      integer count, count2, icount, i,j,k, ip, jp, kp
      double precision si, sumAij, sumAjk, sumAijAjk, sumdAijddijdxi
      double precision sumdAijddijdyi,sumdAijddijdzi,sumdAijddijdxiAjk
      double precision sumdAijddijdyiAjk,sumdAijddijdziAjk, rij,
     &                 tmpaij, Aij, dAijddij
      double precision dAijddijdxj, dAijddijdyj, dAijddijdzj
      double precision sumdAjkddjkdxj,sumdAjkddjkdyj,sumdAjkddjkdzj,
     &                 p3p4Aij
      double precision rjk2, djk1i, rjk, vdw2dif, tmpajk,
     &                 Ajk, sumAjk2, dAjkddjk 
      double precision dAjkddjkdxj, dAjkddjkdyj, dAjkddjkdzj, 
     &                 lastxj, lastyj, lastzj
      double precision dAidxj, dAidyj, dAidzj, Ai, Aidxi, Aidyi,
     &                 Aidzi, dAidxi, dAidyi, dAidzi

      double precision PIx4, PIx2, PIx1
      PARAMETER (PIx4 = 12.5663706143591724639918,
     &           PIx2 =  6.2831853071795862319959,
     &           PIx1 =  3.1415926535897931159979)

 
c --------------------------------------------------------------------- 
c begin
      totsasa = 0.0D0
c     write(*,*) "GB_NONPOL2> egb_nonpol2 called  ..." 

c      write(*,*) "GB_NONPOL> at the begining of  gb_nonpol"
c      do  i=1,npt
c         write(*,*)'GB_NONPOL> i dpot ',i,dpot(1,i),dpot(2,i),dpot(3,i)
c      end do 



c print out the neigbour list for debugging 
c (!!! do not delete, comment out in the production code )
c      count = 1 
c      do i = 1, natom
c         
c          write(*,*) "NB_LIST> -----------------------------"
c          write(*,*) "NB_LIST> Neighours of atom ", i, " are:"
c          do while ( ineighbor(count) .ne. 0 )  
c             write (*, *) "NB_LIST> atom:", ineighbor(count)
c             count = count + 1
c          end do 
c          count = count + 1
c          write(*,*) "NB_LIST2> 1st Neighours of atom ", i, " are:"
c	  write(*,*) "NB_LIST2> ", ineighbor(ineighbor_ptrs(i))
c	  write(*,*) "NB_LIST2> number of neigbor atoms of atom: ", i
c	  write(*,*) "NB_LIST2> ",ineighbor_ptrs(i+1) - 
c     &	                          ineighbor_ptrs(i) - 1
c      end do


c --------------------------------------------------------------------- 
c *** make sure first derivatives array start from zeros for debugging
c     i.e. dpot stores only the 1st derivatives of gb non-polarized term
c     w.r.t. coordinates.
c     should comment out in the production code
c      do i=1, natom
c          dpot(1,i) = 0.0d0
c          dpot(2,i) = 0.0d0
c          dpot(3,i) = 0.0d0
c      end do
c --------------------------------------------------------------------- 


      count = 1 
      do i = 1, natom
          Ai = 0 
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
              vdw2dif = P0(ip)*P0(ip) - P0(jp)*P0(jp)
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

c               write(*, *) "SASA> -------------------------"              
c               write(*, *) "SASA AREA> sasa of Atom", i , " is ", Ai 

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

c here updating the force array
              dpot(1,i) = dpot(1,i) + dAidxi
              dpot(2,i) = dpot(2,i) + dAidyi
              dpot(3,i) = dpot(3,i) + dAidzi


c summing Ai
              totsasa = totsasa + Ai
c              write(*, *) "SASA TOT> totsasa up to ", i ,
c     &                    " is ", totsasa

c YS:  this "end if " corresponding to " if ( ineighbor(count) .eq. 0 ) then"
          end if 
c          write(*, *) "SASA> -------------------------"              
c          write(*, *) "SASA AREA> sasa of Atom", i , " is ", Ai 
      end do


c      write(*,*) "SASA> TOTSASA =  ", totsasa 
       e_gbnp = surften*totsasa
c      write(*,*) "SASA> e_gbnp =  ", e_gbnp 

c      write(*,*) "GB_NONPOL> at the end of gb_nonpol"
c      do  i=1,npt
c         write(*,*)'GB_NONPOL> i dpot ',i,dpot(1,i),dpot(2,i),dpot(3,i)
c      end do 

      return
      end
