      SUBROUTINE egb_np2nd(natom)
c Calculate the second derivatives of
c non-polarized solvation energy 
c

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SGB2DRV.BLOCK'
c args
      integer natom

c     local variables used to calculate 2nd derivatives declared here:
      double precision  r2, dij1i
      double precision xi, yi, zi, xj, yj, zj, xk, yk, zk
      integer count,  i,j,k, ip, jp, kp
      double precision sumdAijddijdyiAjk,sumdAijddijdziAjk, rij,
     &                 tmpaij, Aij, dAijddij
      double precision dAijddijdxj, dAijddijdyj, dAijddijdzj
      double precision sumdAjkddjkdxj,sumdAjkddjkdyj,sumdAjkddjkdzj,
     &                 p3p4Aij
      double precision rjk2, djk1i, rjk, vdw2dif, tmpajk,
     &                 Ajk, sumAjk, sumAjk2, dAjkddjk 

      double precision PIx4, PIx2, PIx1
      PARAMETER (PIx4 = 12.5663706143591724639918,
     &           PIx2 =  6.2831853071795862319959,
     &           PIx1 =  3.1415926535897931159979)

c.....  for case 1
      double precision 
     & d2Aidxdx,d2Aijdxdx,sumd2Aijdxdx, c1term1_xx, 
     & d2Aidxdy,d2Aijdxdy,sumd2Aijdxdy, c1term1_xy, 
     & d2Aidxdz,d2Aijdxdz,sumd2Aijdxdz, c1term1_xz, 
     & d2Aidydy,d2Aijdydy,sumd2Aijdydy, c1term1_yy, 
     & d2Aidydz,d2Aijdydz,sumd2Aijdydz, c1term1_yz, 
     & d2Aidzdz,d2Aijdzdz,sumd2Aijdzdz, c1term1_zz 
      double precision
     &  d2dijdxdx,sumc1term1_xx,
     &  d2dijdydy,sumc1term1_xy,
     &  d2dijdzdz,sumc1term1_xz,
     &  d2dijdxdy,sumc1term1_yy,
     &  d2dijdydz,sumc1term1_yz,
     &  d2dijdxdz,sumc1term1_zz


      double precision d2Aijddijddij, 
     &            ddijdxi, ddijdxj,
     &            ddijdyi, ddijdyj,
     &            ddijdzi, ddijdzj
      integer jj, jstart, jend
      integer kk, kstart, kend, indx, jndx, kndx

c..... for case 2
      double precision ddjkdxj, ddjkdyj,ddjkdzj,
     &                 ddjkdxk, ddjkdyk,ddjkdzk
      double precision d2Ajkddjkddjk, 
     &  sumkdAjkdxj, dAjkdxj,dAjkdxk,dAijdxj,dAijdxi,  
     &  sumkdAjkdyj, dAjkdyj,dAjkdyk,dAijdyj,dAijdyi,  
     &  sumkdAjkdzj, dAjkdzj,dAjkdzk,dAijdzj,dAijdzi

      double precision
     & d2Aidxjdxj, sumd2Aijdxjdxj,d2Aijdxjdxj,
     & d2Aidxjdyj, sumd2Aijdyjdyj,d2Aijdyjdyj,
     & d2Aidxjdzj, sumd2Aijdzjdzj,d2Aijdzjdzj,
     & d2Aidyjdyj, sumd2Aijdxjdyj,d2Aijdxjdyj,
     & d2Aidyjdzj, sumd2Aijdxjdzj,d2Aijdxjdzj,
     & d2Aidzjdzj, sumd2Aijdyjdzj,d2Aijdyjdzj
      double precision
     & sumkd2Ajkdxjdxj, d2Ajkdxjdxj,d2djkdxjdxj,
     & sumkd2Ajkdxjdyj, d2Ajkdxjdyj,d2djkdyjdyj,
     & sumkd2Ajkdxjdzj, d2Ajkdxjdzj,d2djkdzjdzj,
     & sumkd2Ajkdyjdyj, d2Ajkdyjdyj,d2djkdxjdyj,
     & sumkd2Ajkdyjdzj, d2Ajkdyjdzj,d2djkdxjdzj,
     & sumkd2Ajkdzjdzj, d2Ajkdzjdzj,d2djkdyjdzj

       double precision 
     & c2term1_xx, c2term2_xx, c2term3_xx, sumc2term123_xx, 
     & c2term1_xy, c2term2_xy, c2term3_xy, sumc2term123_xy, 
     & c2term1_xz, c2term2_xz, c2term3_xz, sumc2term123_xz, 
     & c2term1_yy, c2term2_yy, c2term3_yy, sumc2term123_yy, 
     & c2term1_yz, c2term2_yz, c2term3_yz, sumc2term123_yz, 
     & c2term1_zz, c2term2_zz, c2term3_zz, sumc2term123_zz 

c....... for case (3)
       double precision  
     & d2Ajkdxkdxk,d2Aidxkdxk,
     & d2Ajkdxkdyk,d2Aidxkdyk,
     & d2Ajkdxkdzk,d2Aidxkdzk,
     & d2Ajkdykdyk,d2Aidykdyk,
     & d2Ajkdykdzk,d2Aidykdzk,
     & d2Ajkdzkdzk,d2Aidzkdzk
    
       double precision rik,rik2, dAikddik, d2Aikddikddik,dik1i, 
     & d2dikdxkdxk, d2Aikdxkdxk 
              
       double precision 
     & ddikdxk  
c ---------------------------------------------------------------------      
c...... off diagonal terms 
c...... for case (4) off-diagonal xi, xj 
       integer pairndx, overlappair_idx 
       double precision
     &  d2Aidxidxj, d2Aijdxidxj, d2dijdxidxj,
     &  d2Aidxidyj, d2Aijdxidyj, d2dijdxidyj,
     &  d2Aidxidzj, d2Aijdxidzj, d2dijdxidzj,
     &  d2Aidyidyj, d2Aijdyidyj, d2dijdyidyj,
     &  d2Aidyidzj, d2Aijdyidzj, d2dijdyidzj,
     &  d2Aidzidzj, d2Aijdzidzj, d2dijdzidzj
       double precision  
     &   c4term1_xx,  c4term2_xx,  
     &   c4term1_xy,  c4term2_xy,  
     &   c4term1_xz,  c4term2_xz,  
     &   c4term1_yy,  c4term2_yy,  
     &   c4term1_yz,  c4term2_yz,  
     &   c4term1_zz,  c4term2_zz  

c...... for case (5) off-diag  xi, xk
       double precision  
     &  d2Aidxidxk, c5term1_xx,
     &  d2Aidxidyk, c5term1_xy,
     &  d2Aidxidzk, c5term1_xz,
     &  d2Aidyidyk, c5term1_yy,
     &  d2Aidyidzk, c5term1_yz,
     &  d2Aidzidzk, c5term1_zz

c...... for case (6) off-daig  xj, xk
       double precision 
     &  d2Aidxjdxk,c6term1_xx,c6term2_xx,d2djkdxjdxk, d2Ajkdxjdxk,
     &  d2Aidxjdyk,c6term1_xy,c6term2_xy,d2djkdxjdyk, d2Ajkdxjdyk,
     &  d2Aidxjdzk,c6term1_xz,c6term2_xz,d2djkdxjdzk, d2Ajkdxjdzk,
     &  d2Aidyjdyk,c6term1_yy,c6term2_yy,d2djkdyjdyk, d2Ajkdyjdyk,
     &  d2Aidyjdzk,c6term1_yz,c6term2_yz,d2djkdyjdzk, d2Ajkdyjdzk,
     &  d2Aidzjdzk,c6term1_zz,c6term2_zz,d2djkdzjdzk, d2Ajkdzjdzk

c begin           

c     write(*,*) "GB_NONPOL2> calcate 2nd derivative:"
              
              
c (1) start with the current atom i, start i loop 

c --------------------------------------------------------------------- 
c make sure diag is clear for debugging 
c comment out in production code  
c      do i = 1, 6 * natom
c          diag(i) = 0.0d0
c      end do
c.....initialize the offdiaggbnp 
      do i=1, 180 * natom
          offdiaggbnp(i) = 0.0d0
      end do

c --------------------------------------------------------------------- 
      do i = 1,natom
         ip = ptsaid(i)
c........case 1: diagnol, wrt. to atom i 
         sumd2Aijdxdx = 0.0d0
         sumd2Aijdxdy = 0.0d0
         sumd2Aijdxdz = 0.0d0
         sumd2Aijdydy = 0.0d0
         sumd2Aijdydz = 0.0d0
         sumd2Aijdzdz = 0.0d0

         sumd2Aijdxjdxj = 0.0d0
         sumd2Aijdxjdyj = 0.0d0
         sumd2Aijdxjdzj = 0.0d0
         sumd2Aijdyjdyj = 0.0d0
         sumd2Aijdyjdzj = 0.0d0
         sumd2Aijdzjdzj = 0.0d0

         sumc1term1_xx = 0.0d0 
         sumc1term1_xy = 0.0d0 
         sumc1term1_xz = 0.0d0 
         sumc1term1_yy = 0.0d0 
         sumc1term1_yz = 0.0d0 
         sumc1term1_zz = 0.0d0 

c........find the next neighbour atom of i (atom j) 
         jstart = ineighbor_ptrs(i)
         jend = ineighbor_ptrs(i+1) - 2

c (2)... start j loop 

         do jj = jstart, jend
            j = ineighbor(jj)
            jp = ptsaid(j)
c..... 	    write(*,*) "GB_NONPOL2> next neigbhor of atom:", i, "is", j
            
            xi = coor(1,i)
            yi = coor(2,i)
            zi = coor(3,i)
            xj = coor(1,j)
            yj = coor(2,j)
            zj = coor(3,j)
            r2 = (xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj)
            dij1i = 1.0d0 / DSQRT(r2) 
	    rij = r2 * dij1i
            vdw2dif = P0(ip) * P0(ip) - P0(jp) * P0(jp)          
	    tmpaij = 2.0d0 * P0(ip) - rij - vdw2dif * dij1i
   	    Aij = PIx1 * P0(ip) * tmpaij
c (3)....    calculate d2Aijdxdx etc 6 terms
c (3.1) ....  calculate d2dijdxdx etc for case 1 and 2 

            d2dijdxdx = -(xj-xi)*(xj-xi)*dij1i*dij1i*dij1i+dij1i
            d2dijdxdy = -(xj-xi)*(yj-yi)*dij1i*dij1i*dij1i
            d2dijdxdz = -(xj-xi)*(zj-zi)*dij1i*dij1i*dij1i
            d2dijdydy = -(yj-yi)*(yj-yi)*dij1i*dij1i*dij1i+dij1i
            d2dijdydz = -(yj-yi)*(zj-zi)*dij1i*dij1i*dij1i
            d2dijdzdz = -(zj-zi)*(zj-zi)*dij1i*dij1i*dij1i+dij1i

c........... for off-diagonal terms
            d2dijdxidxj=(xj-xi)*(xj-xi)*dij1i*dij1i*dij1i-dij1i 
            d2dijdxidyj=(xj-xi)*(yj-yi)*dij1i*dij1i*dij1i 
            d2dijdxidzj=(xj-xi)*(zj-zi)*dij1i*dij1i*dij1i 
            d2dijdyidyj=(yj-yi)*(yj-yi)*dij1i*dij1i*dij1i-dij1i 
            d2dijdyidzj=(yj-yi)*(zj-zi)*dij1i*dij1i*dij1i 
            d2dijdzidzj=(zj-zi)*(zj-zi)*dij1i*dij1i*dij1i-dij1i 

c (3.2) ..  calculate dAijddij
            dAijddij = PIx1 * P0(ip) * (dij1i*dij1i*vdw2dif - 1.0d0)

c (3.3)     calculate ddijdxi etc
            ddijdxi = -(xj-xi) * dij1i
            ddijdxj = -ddijdxi
	    ddijdyi = -(yj-yi) * dij1i
            ddijdyj = -ddijdyi
	    ddijdzi = -(zj-zi) * dij1i
            ddijdzj = -ddijdzi

c.... for case (2)  
            dAijdxj = dAijddij * ddijdxj
            dAijdyj = dAijddij * ddijdyj
            dAijdzj = dAijddij * ddijdzj

	    dAijdxi = dAijddij * ddijdxi 
	    dAijdyi = dAijddij * ddijdyi 
	    dAijdzi = dAijddij * ddijdzi 
            
c (3.4)     calculate d2Aijdijdij etc
            d2Aijddijddij= -PIx2 * P0(ip) * vdw2dif *
     &	         dij1i * dij1i * dij1i  
c (3.5) based on eqn. A.7
            d2Aijdxdx = dAijddij * d2dijdxdx + ddijdxi * ddijdxi * 
     &                 d2Aijddijddij 
            d2Aijdxdy = dAijddij * d2dijdxdy + ddijdxi * ddijdyi * 
     &                 d2Aijddijddij 
            d2Aijdxdz = dAijddij * d2dijdxdz + ddijdxi * ddijdzi * 
     &                 d2Aijddijddij 
            d2Aijdydy = dAijddij * d2dijdydy + ddijdyi * ddijdyi * 
     &                 d2Aijddijddij 
            d2Aijdydz = dAijddij * d2dijdydz + ddijdyi * ddijdzi * 
     &                 d2Aijddijddij 
            d2Aijdzdz = dAijddij * d2dijdzdz + ddijdzi * ddijdzi * 
     &                 d2Aijddijddij 

c.......... for off-diagonal terms
            d2Aijdxidxj = dAijddij * d2dijdxidxj +  
     &                    ddijdxi * ddijdxj * d2Aijddijddij       
            d2Aijdxidyj = dAijddij * d2dijdxidyj +  
     &                    ddijdxi * ddijdyj * d2Aijddijddij       
            d2Aijdxidzj = dAijddij * d2dijdxidzj +  
     &                    ddijdxi * ddijdzj * d2Aijddijddij       
            d2Aijdyidyj = dAijddij * d2dijdyidyj +  
     &                    ddijdyi * ddijdyj * d2Aijddijddij       
            d2Aijdzidzj = dAijddij * d2dijdzidzj +  
     &                    ddijdzi * ddijdzj * d2Aijddijddij       
            d2Aijdyidzj = dAijddij * d2dijdyidzj +  
     &                    ddijdyi * ddijdzj * d2Aijddijddij       

c ** for case (2)
            d2Aijdxjdxj = dAijddij * d2dijdxdx + ddijdxj * ddijdxj * 
     &                 d2Aijddijddij 
            d2Aijdxjdyj = dAijddij * d2dijdxdy + ddijdxj * ddijdyj * 
     &                 d2Aijddijddij 
            d2Aijdxjdzj = dAijddij * d2dijdxdz + ddijdxj * ddijdzj * 
     &                 d2Aijddijddij 
            d2Aijdyjdyj = dAijddij * d2dijdydy + ddijdyj * ddijdyj * 
     &                 d2Aijddijddij 
            d2Aijdyjdzj = dAijddij * d2dijdydz + ddijdyj * ddijdzj * 
     &                 d2Aijddijddij 
            d2Aijdzjdzj = dAijddij * d2dijdzdz + ddijdzj * ddijdzj * 
     &                 d2Aijddijddij 
c (4)
            sumd2Aijdxdx = sumd2Aijdxdx + d2Aijdxdx
            sumd2Aijdxdy = sumd2Aijdxdy + d2Aijdxdy
            sumd2Aijdxdz = sumd2Aijdxdz + d2Aijdxdz
            sumd2Aijdydy = sumd2Aijdydy + d2Aijdydy
            sumd2Aijdydz = sumd2Aijdydz + d2Aijdydz
            sumd2Aijdzdz = sumd2Aijdzdz + d2Aijdzdz

            sumd2Aijdxjdxj = sumd2Aijdxjdxj + d2Aijdxjdxj
            sumd2Aijdxjdyj = sumd2Aijdxjdyj + d2Aijdxjdyj
            sumd2Aijdxjdzj = sumd2Aijdxjdzj + d2Aijdxjdzj
            sumd2Aijdyjdyj = sumd2Aijdyjdyj + d2Aijdyjdyj
            sumd2Aijdyjdzj = sumd2Aijdyjdzj + d2Aijdyjdzj
            sumd2Aijdzjdzj = sumd2Aijdzjdzj + d2Aijdzjdzj

c (5) k loop
            sumAjk = 0.0d0 
            sumkdAjkdxj = 0.0d0
            sumkdAjkdyj = 0.0d0
            sumkdAjkdzj = 0.0d0

            sumkd2Ajkdxjdxj = 0.0d0
            sumkd2Ajkdxjdyj = 0.0d0
            sumkd2Ajkdxjdzj = 0.0d0
            sumkd2Ajkdyjdyj = 0.0d0
            sumkd2Ajkdyjdzj = 0.0d0
            sumkd2Ajkdzjdzj = 0.0d0

            kstart = ineighbor_ptrs(i)
            kend = ineighbor_ptrs(i+1)-2
            do 100 kk = kstart, kend
	        k = ineighbor(kk)
	        if ( k .eq. j ) goto 100  
		kp = ptsaid(k)
		xk = coor(1,k)
		yk = coor(2,k)
		zk = coor(3,k)
		rjk2 = 
     &               (xj - xk) * (xj - xk) +
     &               (yj - yk) * (yj - yk) +
     &               (zj - zk) * (zj - zk) 
                djk1i  = 1.0d0/DSQRT(rjk2)
                rjk = rjk2 * djk1i 
                if ( (P0(jp) + P0(kp)) .gt. rjk ) then
c (5.1)            calculate sumAjk
                   vdw2dif = P0(jp) * P0(jp) - P0(kp) * P0(kp) 
		   tmpajk = 2.0d0 * P0(jp) - rjk - vdw2dif * djk1i
		   Ajk = PIx1 * P0(jp) * tmpajk
		   sumAjk = sumAjk + Ajk


c (5.2) .......... calculate dAjkdxj etc 3 term
                   dAjkddjk = PIx1 * P0(jp) * 
     &		              (djk1i*djk1i*vdw2dif - 1.0d0)
                   ddjkdxj = -(xk - xj)* djk1i
                   ddjkdxk = - ddjkdxj
                   ddjkdyj = -(yk - yj)* djk1i
                   ddjkdyk = - ddjkdyj
                   ddjkdzj = -(zk - zj)* djk1i
                   ddjkdzk = - ddjkdzj

                   dAjkdxj = dAjkddjk * ddjkdxj
                   dAjkdxk = dAjkddjk * ddjkdxk
		   dAjkdyj = dAjkddjk * ddjkdyj
		   dAjkdyk = dAjkddjk * ddjkdyk
		   dAjkdzj = dAjkddjk * ddjkdzj
		   dAjkdzk = dAjkddjk * ddjkdzk

		   sumkdAjkdxj = sumkdAjkdxj + dAjkdxj
		   sumkdAjkdyj = sumkdAjkdyj + dAjkdyj
		   sumkdAjkdzj = sumkdAjkdzj + dAjkdzj


c (5.3) .......... calculate d2Ajkdxjdxj etc 6 term                   
c (5.3.1) ......... d2djkdxjdxj ect 6 terms
                   d2djkdxjdxj= -(xk-xj)*(xk-xj)*djk1i*djk1i*djk1i+djk1i
		   d2djkdxjdyj= -(xk-xj)*(yk-yj)*djk1i*djk1i*djk1i
		   d2djkdxjdzj= -(xk-xj)*(zk-zj)*djk1i*djk1i*djk1i
                   d2djkdyjdyj= -(yk-yj)*(yk-yj)*djk1i*djk1i*djk1i+djk1i
		   d2djkdyjdzj= -(yk-yj)*(zk-zj)*djk1i*djk1i*djk1i
                   d2djkdzjdzj= -(zk-zj)*(zk-zj)*djk1i*djk1i*djk1i+djk1i
c (5.3.2)
                   d2Ajkddjkddjk = -PIx2 * P0(jp)* vdw2dif*
     &		                    djk1i*djk1i*djk1i
c six terms (A.7)
                   d2Ajkdxjdxj = dAjkddjk * d2djkdxjdxj +
     &                           ddjkdxj * ddjkdxj * d2Ajkddjkddjk 
                   d2Ajkdxjdyj = dAjkddjk * d2djkdxjdyj +
     &                           ddjkdxj * ddjkdyj * d2Ajkddjkddjk 
                   d2Ajkdxjdzj = dAjkddjk * d2djkdxjdzj +
     &                           ddjkdxj * ddjkdzj * d2Ajkddjkddjk 
                   d2Ajkdyjdyj = dAjkddjk * d2djkdyjdyj +
     &                           ddjkdyj * ddjkdyj * d2Ajkddjkddjk 
                   d2Ajkdyjdzj = dAjkddjk * d2djkdyjdzj +
     &                           ddjkdyj * ddjkdzj * d2Ajkddjkddjk 
                   d2Ajkdzjdzj = dAjkddjk * d2djkdzjdzj +
     &                           ddjkdzj * ddjkdzj * d2Ajkddjkddjk 
                   
c (5.3.3) ........... add above six terms to the sumkd2Ajkdxjdxj etc
                   sumkd2Ajkdxjdxj = sumkd2Ajkdxjdxj + d2Ajkdxjdxj
                   sumkd2Ajkdxjdyj = sumkd2Ajkdxjdyj + d2Ajkdxjdyj
                   sumkd2Ajkdxjdzj = sumkd2Ajkdxjdzj + d2Ajkdxjdzj
                   sumkd2Ajkdyjdyj = sumkd2Ajkdyjdyj + d2Ajkdyjdyj
                   sumkd2Ajkdyjdzj = sumkd2Ajkdyjdzj + d2Ajkdyjdzj
                   sumkd2Ajkdzjdzj = sumkd2Ajkdzjdzj + d2Ajkdzjdzj

c (5.4)            for case (3) d2Ajkdxkdxk etj
c                  note: d2djkdxjdxj = d2djkdxkdxk
                   d2Ajkdxkdxk = dAjkddjk * d2djkdxjdxj +
     &                           ddjkdxk * ddjkdxk * d2Ajkddjkddjk
                   d2Ajkdxkdyk = dAjkddjk * d2djkdxjdyj +
     &                           ddjkdxk * ddjkdyk * d2Ajkddjkddjk
                   d2Ajkdxkdzk = dAjkddjk * d2djkdxjdzj +
     &                           ddjkdxk * ddjkdzk * d2Ajkddjkddjk
                   d2Ajkdykdyk = dAjkddjk * d2djkdyjdyj +
     &                           ddjkdyk * ddjkdyk * d2Ajkddjkddjk
                   d2Ajkdykdzk = dAjkddjk * d2djkdyjdzj +
     &                           ddjkdyk * ddjkdzk * d2Ajkddjkddjk
                   d2Ajkdzkdzk = dAjkddjk * d2djkdzjdzj +
     &                           ddjkdzk * ddjkdzk * d2Ajkddjkddjk

                   d2Aidxkdxk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdxkdxk
                   d2Aidxkdyk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdxkdyk
                   d2Aidxkdzk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdxkdzk
                   d2Aidykdyk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdykdyk
                   d2Aidykdzk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdykdzk
                   d2Aidzkdzk =  
     &                          (P3(ip) + P4(ip) * Aij) * d2Ajkdzkdzk

c                   write(*,*) "## i,j,k, d2Aidxkdxk", i,j,k,d2Aidxkdxk
c ......... .......update diag array
                   kndx = 6 * (k - 1) + 1 
                   diag(kndx) = diag(kndx) + d2Aidxkdxk * surften
		   kndx = kndx + 1
                   diag(kndx) = diag(kndx) + d2Aidxkdyk * surften
		   kndx = kndx + 1
                   diag(kndx) = diag(kndx) + d2Aidxkdzk * surften
		   kndx = kndx + 1
                   diag(kndx) = diag(kndx) + d2Aidykdyk * surften
		   kndx = kndx + 1
                   diag(kndx) = diag(kndx) + d2Aidykdzk * surften
		   kndx = kndx + 1
                   diag(kndx) = diag(kndx) + d2Aidzkdzk * surften

c --------------------------------------------------------------------- 
c.............     starts off-diagonal terms 
c ................for case (5) off-diagonal i, k pair
                  c5term1_xx = dAijdxi * dAjkdxk
                  c5term1_yy = dAijdyi * dAjkdyk
                  c5term1_zz = dAijdzi * dAjkdzk

		  if ( i .lt. k) then 
                      c5term1_xy = dAijdxi * dAjkdyk
                      c5term1_xz = dAijdxi * dAjkdzk
                      c5term1_yz = dAijdyi * dAjkdzk
                  else
                      c5term1_xy = dAijdyi * dAjkdxk
                      c5term1_xz = dAijdzi * dAjkdxk
                      c5term1_yz = dAijdzi * dAjkdyk
                  end if
		    

		  d2Aidxidxk = P4(ip) * c5term1_xx
		  d2Aidxidyk = P4(ip) * c5term1_xy
		  d2Aidxidzk = P4(ip) * c5term1_xz
		  d2Aidyidyk = P4(ip) * c5term1_yy
		  d2Aidyidzk = P4(ip) * c5term1_yz
		  d2Aidzidzk = P4(ip) * c5term1_zz
c		  write(*,*) "#5 i,k,d2Aidxidyk:",i,k,d2Aidxidyk

		  pairndx = overlappair_idx(i,k)
		  indx = 6 * (pairndx - 1) + 1
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidxidxk * surften
		  indx = indx + 1
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidxidyk * surften
		  indx = indx + 1
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidxidzk * surften

		  indx = indx + 1 
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidyidyk * surften
		  indx = indx + 1 
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidyidzk * surften

		  indx = indx + 1 
		  offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &		                      d2Aidzidzk * surften

c ................for case (6) off-daiagonal j,k pair
		  d2djkdxjdxk = (xk - xj) * (xk - xj) *
     &		                djk1i * djk1i * djk1i - djk1i
		  d2djkdxjdyk = (xk - xj) * (yk - yj) *
     &		                djk1i * djk1i * djk1i 
		  d2djkdxjdzk = (xk - xj) * (zk - zj) *
     &		                djk1i * djk1i * djk1i 
		  d2djkdyjdyk = (yk - yj) * (yk - yj) *
     &		                djk1i * djk1i * djk1i - djk1i
		  d2djkdyjdzk = (yk - yj) * (zk - zj) *
     &		                djk1i * djk1i * djk1i 
		  d2djkdzjdzk = (zk - zj) * (zk - zj) *
     &		                djk1i * djk1i * djk1i - djk1i

                  d2Ajkddjkddjk = -PIx2 * P0(jp) * 
     &		                   (vdw2dif * djk1i * djk1i * djk1i)
c ............... eqn A.7
                  d2Ajkdxjdxk = dAjkddjk * d2djkdxjdxk +
     &                          ddjkdxj * ddjkdxk * d2Ajkddjkddjk
                  d2Ajkdxjdyk = dAjkddjk * d2djkdxjdyk +
     &                          ddjkdxj * ddjkdyk * d2Ajkddjkddjk
                  d2Ajkdxjdzk = dAjkddjk * d2djkdxjdzk +
     &                          ddjkdxj * ddjkdzk * d2Ajkddjkddjk
                  d2Ajkdyjdyk = dAjkddjk * d2djkdyjdyk +
     &                          ddjkdyj * ddjkdyk * d2Ajkddjkddjk
                  d2Ajkdyjdzk = dAjkddjk * d2djkdyjdzk +
     &                          ddjkdyj * ddjkdzk * d2Ajkddjkddjk
                  d2Ajkdzjdzk = dAjkddjk * d2djkdzjdzk +
     &                          ddjkdzj * ddjkdzk * d2Ajkddjkddjk

		  c6term1_xx = dAijdxj * dAjkdxk
		  c6term1_yy = dAijdyj * dAjkdyk
		  c6term1_zz = dAijdzj * dAjkdzk
                  if ( j .lt. k) then 
		      c6term1_xy = dAijdxj * dAjkdyk
		      c6term1_xz = dAijdxj * dAjkdzk
		      c6term1_yz = dAijdyj * dAjkdzk
                  else
		      c6term1_xy = dAijdyj * dAjkdxk
		      c6term1_xz = dAijdzj * dAjkdxk
		      c6term1_yz = dAijdzj * dAjkdyk
		  end if   


		  c6term2_xx = Aij * d2Ajkdxjdxk
		  c6term2_xy = Aij * d2Ajkdxjdyk
		  c6term2_xz = Aij * d2Ajkdxjdzk
		  c6term2_yy = Aij * d2Ajkdyjdyk
		  c6term2_yz = Aij * d2Ajkdyjdzk
		  c6term2_zz = Aij * d2Ajkdzjdzk

		  d2Aidxjdxk = P3(ip) * d2Ajkdxjdxk + P4(ip) * 
     &                         ( c6term1_xx + c6term2_xx)
		  d2Aidxjdyk = P3(ip) * d2Ajkdxjdyk + P4(ip) * 
     &                         ( c6term1_xy + c6term2_xy)
		  d2Aidxjdzk = P3(ip) * d2Ajkdxjdzk + P4(ip) * 
     &                         ( c6term1_xz + c6term2_xz)
		  d2Aidyjdyk = P3(ip) * d2Ajkdyjdyk + P4(ip) * 
     &                         ( c6term1_yy + c6term2_yy)
		  d2Aidyjdzk = P3(ip) * d2Ajkdyjdzk + P4(ip) * 
     &                         ( c6term1_yz + c6term2_yz)
		  d2Aidzjdzk = P3(ip) * d2Ajkdzjdzk + P4(ip) * 
     &                         ( c6term1_zz + c6term2_zz)
c		  write(*,*) "#6 j,k,d2Aidxjdyk:",j,k,d2Aidxjdyk
      
                  pairndx = overlappair_idx(j,k)
		  indx = 6 * (pairndx - 1) +1
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidxjdxk * surften 
		  indx = indx + 1
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidxjdyk * surften 
		  indx = indx + 1 
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidxjdzk * surften 
		  indx = indx + 1
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidyjdyk * surften 
		  indx = indx + 1
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidyjdzk * surften 
		  indx = indx + 1 
		  offdiaggbnp(indx) =  offdiaggbnp(indx) + 
     &                   	       d2Aidzjdzk * surften 
                   
                end if
c   k loop ended
100	    continue 
c --------------------------------------------------------------------- 
c            write(*,*) "******************"

c (6) 
c (6.1)
            c1term1_xx = d2Aijdxdx * sumAjk
            c1term1_xy = d2Aijdxdy * sumAjk
            c1term1_xz = d2Aijdxdz * sumAjk
            c1term1_yy = d2Aijdydy * sumAjk
            c1term1_yz = d2Aijdydz * sumAjk
            c1term1_zz = d2Aijdzdz * sumAjk

            sumc1term1_xx = sumc1term1_xx + c1term1_xx
            sumc1term1_xy = sumc1term1_xy + c1term1_xy
            sumc1term1_xz = sumc1term1_xz + c1term1_xz
            sumc1term1_yy = sumc1term1_yy + c1term1_yy
            sumc1term1_yz = sumc1term1_yz + c1term1_yz
            sumc1term1_zz = sumc1term1_zz + c1term1_zz
c (6.2)
            c2term1_xx = Aij * sumkd2Ajkdxjdxj
            c2term1_xy = Aij * sumkd2Ajkdxjdyj
            c2term1_xz = Aij * sumkd2Ajkdxjdzj
            c2term1_yy = Aij * sumkd2Ajkdyjdyj
            c2term1_yz = Aij * sumkd2Ajkdyjdzj
            c2term1_zz = Aij * sumkd2Ajkdzjdzj

            c2term2_xx = dAijdxj * sumkdAjkdxj + dAijdxj * sumkdAjkdxj
            c2term2_xy = dAijdxj * sumkdAjkdyj + dAijdyj * sumkdAjkdxj
            c2term2_xz = dAijdxj * sumkdAjkdzj + dAijdzj * sumkdAjkdxj
            c2term2_yy = dAijdyj * sumkdAjkdyj + dAijdyj * sumkdAjkdyj
            c2term2_yz = dAijdyj * sumkdAjkdzj + dAijdzj * sumkdAjkdyj
            c2term2_zz = dAijdzj * sumkdAjkdzj + dAijdzj * sumkdAjkdzj

            c2term3_xx = d2Aijdxjdxj * sumAjk
            c2term3_xy = d2Aijdxjdyj * sumAjk
            c2term3_xz = d2Aijdxjdzj * sumAjk
            c2term3_yy = d2Aijdyjdyj * sumAjk
            c2term3_yz = d2Aijdyjdzj * sumAjk
            c2term3_zz = d2Aijdzjdzj * sumAjk

            sumc2term123_xx =  c2term1_xx +  c2term2_xx + c2term3_xx
            sumc2term123_xy =  c2term1_xy +  c2term2_xy + c2term3_xy
            sumc2term123_xz =  c2term1_xz +  c2term2_xz + c2term3_xz
            sumc2term123_yy =  c2term1_yy +  c2term2_yy + c2term3_yy
            sumc2term123_yz =  c2term1_yz +  c2term2_yz + c2term3_yz
            sumc2term123_zz =  c2term1_zz +  c2term2_zz + c2term3_zz

c (6.3) ... for case (2) wrt j atom, update diag array
            d2Aidxjdxj = P2(ip)* d2Aijdxjdxj + P3(ip)*sumkd2Ajkdxjdxj +
     &                   P4(ip)*sumc2term123_xx 
            d2Aidxjdyj = P2(ip)* d2Aijdxjdyj + P3(ip)*sumkd2Ajkdxjdyj +
     &                   P4(ip)*sumc2term123_xy 
            d2Aidxjdzj = P2(ip)* d2Aijdxjdzj + P3(ip)*sumkd2Ajkdxjdzj +
     &                   P4(ip)*sumc2term123_xz 
            d2Aidyjdyj = P2(ip)* d2Aijdyjdyj + P3(ip)*sumkd2Ajkdyjdyj +
     &                   P4(ip)*sumc2term123_yy 
            d2Aidyjdzj = P2(ip)* d2Aijdyjdzj + P3(ip)*sumkd2Ajkdyjdzj +
     &                   P4(ip)*sumc2term123_yz 
            d2Aidzjdzj = P2(ip)* d2Aijdzjdzj + P3(ip)*sumkd2Ajkdzjdzj +
     &                   P4(ip)*sumc2term123_zz 
c            write(*,*) "## i, j, d2Aidxjdxj", i,j, d2Aidxjdxj

c ......... update diag array
            jndx = 6 * (j - 1) + 1 
            diag(jndx) = diag(jndx) + d2Aidxjdxj * surften
	    jndx = jndx + 1
            diag(jndx) = diag(jndx) + d2Aidxjdyj * surften
	    jndx = jndx + 1
            diag(jndx) = diag(jndx) + d2Aidxjdzj * surften
	    jndx = jndx + 1
            diag(jndx) = diag(jndx) + d2Aidyjdyj * surften
	    jndx = jndx + 1
            diag(jndx) = diag(jndx) + d2Aidyjdzj * surften
	    jndx = jndx + 1
            diag(jndx) = diag(jndx) + d2Aidzjdzj * surften

c .......... starts off-diagonal terms
c ..........case (4) off diagonal i j pair   
	    c4term1_xx = dAijdxi * sumkdAjkdxj
	    c4term1_yy = dAijdyi * sumkdAjkdyj
	    c4term1_zz = dAijdzi * sumkdAjkdzj

c !!! very tricky here
            if (j .gt. i) then
		c4term1_xy = dAijdxi * sumkdAjkdyj
	        c4term1_xz = dAijdxi * sumkdAjkdzj
	        c4term1_yz = dAijdyi * sumkdAjkdzj
            else 
	        c4term1_xy = dAijdyi * sumkdAjkdxj
	        c4term1_xz = dAijdzi * sumkdAjkdxj
	        c4term1_yz = dAijdzi * sumkdAjkdyj
            end if


	    c4term2_xx =  d2Aijdxidxj * sumAjk
	    c4term2_xy =  d2Aijdxidyj * sumAjk
	    c4term2_xz =  d2Aijdxidzj * sumAjk
	    c4term2_yy =  d2Aijdyidyj * sumAjk
	    c4term2_yz =  d2Aijdyidzj * sumAjk
	    c4term2_zz =  d2Aijdzidzj * sumAjk

	    d2Aidxidxj = P2(ip) * d2Aijdxidxj + 
     &	                 P4(ip) * ( c4term1_xx + c4term2_xx)
	    d2Aidxidyj = P2(ip) * d2Aijdxidyj + 
     &	                 P4(ip) * ( c4term1_xy + c4term2_xy)
	    d2Aidxidzj = P2(ip) * d2Aijdxidzj + 
     &	                 P4(ip) * ( c4term1_xz + c4term2_xz)
	    d2Aidyidyj = P2(ip) * d2Aijdyidyj + 
     &	                 P4(ip) * ( c4term1_yy + c4term2_yy)
	    d2Aidyidzj = P2(ip) * d2Aijdyidzj + 
     &	                 P4(ip) * ( c4term1_yz + c4term2_yz)
	    d2Aidzidzj = P2(ip) * d2Aijdzidzj + 
     &	                 P4(ip) * ( c4term1_zz + c4term2_zz)
            
c            write(*,*) "## i,j,d2Aidxidyj",i,j,d2Aidxidyj
c            write(*,*) "## i,j,c4term1_xy",i,j,c4term1_xy
c            write(*,*) "## i,j,c4term2_xy",i,j,c4term2_xy

            pairndx = overlappair_idx(i,j)
	    indx = 6 * (pairndx - 1) +1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidxidxj * surften
	    indx = indx + 1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidxidyj * surften
	    indx = indx + 1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidxidzj * surften
	    indx = indx + 1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidyidyj * surften
	    indx = indx + 1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidyidzj * surften
	    indx = indx + 1 
	    offdiaggbnp(indx) = offdiaggbnp(indx) + 
     &	                        d2Aidzidzj * surften

         end do 
c (7)    end do jj.  j loop ended
c --------------------------------------------------------------------- 

c (8) update diag array
c (8.1) for atom i case (1)
         d2Aidxdx = P2(ip) * sumd2Aijdxdx + P4(ip) * sumc1term1_xx
         d2Aidxdy = P2(ip) * sumd2Aijdxdy + P4(ip) * sumc1term1_xy
         d2Aidxdz = P2(ip) * sumd2Aijdxdz + P4(ip) * sumc1term1_xz
         d2Aidydy = P2(ip) * sumd2Aijdydy + P4(ip) * sumc1term1_yy
         d2Aidydz = P2(ip) * sumd2Aijdydz + P4(ip) * sumc1term1_yz
         d2Aidzdz = P2(ip) * sumd2Aijdzdz + P4(ip) * sumc1term1_zz

c	 write(*,*) "## i, d2Aidxdx", i, d2Aidxdx

	 indx = 6 * (i - 1) + 1 
	 diag(indx) = diag(indx) + d2Aidxdx * surften
	 indx = indx + 1
	 diag(indx) = diag(indx) + d2Aidxdy * surften
	 indx = indx + 1
	 diag(indx) = diag(indx) + d2Aidxdz * surften
	 indx = indx + 1
	 diag(indx) = diag(indx) + d2Aidydy * surften
	 indx = indx + 1
	 diag(indx) = diag(indx) + d2Aidydz * surften
	 indx = indx + 1
	 diag(indx) = diag(indx) + d2Aidzdz * surften

c.....end of do i=1, natom
      end do
      return
      END
