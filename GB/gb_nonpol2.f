c gb_nonpol2.f
c    computing non-polar solvation energy in GB model
c    and also the 2nd derivatives
c
c    YS, Nov 17, 2005
c
c refs: Weiser et al JCC 2,217(1999)

 
      SUBROUTINE init_gb_nonpol
c     Initialization of parameters for nonpolarized solvation
c     energy calculation.

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PROPERT.BLOCK'
    
      logical SASA_NO_NLR

c begin
      write(*,*) "YS GB_NONPOL> initial gb_nonpol  ..." 
c 
c neigbouring list reduction (NLR) is not implemented 
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

c oxygen in hydroxyl,sp3, bonded neighbor 1, said = 8
          P0(8) =  1.60 + 1.40
          P1(8) =  0.77914 
          P2(8) = -0.25262  
          P3(8) = -0.0016056 
          P4(8) =  0.00035071  

c oxygen ,sp3, bonded neighbor 2, said = 9 
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
c  Not implemented yet - YS (May 16, 2005)
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
      return
      END 

c ------------------------------------------------------------------
      SUBROUTINE egb_np2nd(natom)
c Calculate the secondary derivatives of
c non-polarized solvation energy and 
c

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/SGB2DRV.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
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

c --------------------------------------------------------------------- 
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
      include 'COMMON/GETVEC.BLOCK'

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

c --------------------------------------------------------------------- 
      integer FUNCTION overlappair_idx(i,j) 
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
c args
      integer i, j
c locals
      integer tempi, tempj, temp, jstart, jend, jj
c begin
      tempi = i
      tempj = j
      if ( tempj .eq. tempi ) then
          overlappair_idx = 0
	  return
      end if

      if ( tempj .le. tempi ) then
          temp = tempi
	  tempi = tempj
	  tempj = temp
      end if

      jstart = overlappair_ptrs(tempi)
      jend = overlappair_ptrs(tempi+1) - 1
      if ( jend .lt. jstart ) then
          overlappair_idx = 0
	  return
      end if 

      do jj = jstart, jend 
          if ( overlappair(jj) .eq. tempj ) then
	     overlappair_idx = jj
	     return
          end if
      end do

      return
      END 
c --------------------------------------------------------------------- 
      SUBROUTINE testderi2nd(natom)

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

c local variables
      integer  i, j,indx, overlappair_idx,k
      double precision e_gbnp, e_gbnp1, e_gbnp2, e_gbnp3, e_gbnp4,
     &   diag_xx, diag_xy, diag_xz,
     &   diag_yy, diag_yz, diag_zz
      double precision tempx, tempy, tempz, dx, dy, dz 
      double precision anal_1st1, anal_1st2, num_2nd, d1iff
      double precision gbnp1, gbnp2, num_1st, anal_2nd
c begin
      dx = 0.00000001d0
      dy = 0.00000001d0
      dz = 0.00000001d0
C@
	write(*,*) "maxpt2d = ",maxpt2d

      do i=1, natom
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xx for atom", i  
         tempx = coor(1,i)
         coor(1,i) = tempx + dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(1,i)
         coor(1,i) = tempx - dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(1,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dx) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dx)

         coor(1,i) = tempx 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 1   
         diag_xx = diag(indx) 
         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,*) "xx> 2nd> num,anal", num_2nd,diag_xx  
	 write(*,*) "xx> 2nd> diff", num_2nd - diag_xx 
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xy for atom", i  
         tempy = coor(2,i)
         coor(2,i) = tempy + dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(1,i)
         coor(2,i) = tempy - dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(1,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dy) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dy)

         coor(2,i) = tempy 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 2   
         diag_xy = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,*) "xy> 2nd> num,anal", num_2nd,diag_xy  
	 write(*,*) "xy> 2nd> diff", num_2nd - diag_xy 
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xz for atom", i  
         tempx = coor(1,i)
         coor(1,i) = tempx + dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(3,i)
         coor(1,i) = tempx - dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(3,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dx) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dx)

         coor(1,i) = tempx 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 3   
         diag_xz = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,*) "xz> 2nd> num,anal", num_2nd,diag_xz  
	 write(*,*) "xz> 2nd> diff", num_2nd - diag_xz 
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_yy for atom", i  
         tempy = coor(2,i)

         coor(2,i) = tempy + dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(2,i)

         coor(2,i) = tempy - dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(2,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dy) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dy)

         coor(2,i) = tempy 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 4   
         diag_yy = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,*) "yy> 2nd> num,anal", num_2nd,diag_yy  
	 write(*,*) "yy> 2nd> diff", num_2nd - diag_yy 
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_yz for atom", i  
         tempz = coor(3,i)
         coor(3,i) = tempz + dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(2,i)
         coor(3,i) = tempz - dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(2,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dz) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dz)

         coor(3,i) = tempz 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 5   
         diag_yz = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,*) "yz> 2nd> num,anal", num_2nd,diag_yz  
	 write(*,*) "yz> 2nd> diff", num_2nd - diag_yz 
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_zz for atom", i  
         tempz = coor(3,i)

         coor(3,i) = tempz + dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(3,i)

         coor(3,i) = tempz - dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(3,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dz) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dz)

         coor(3,i) = tempz 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 6   
         diag_zz = diag(indx) 
         write(*,*) "zz> 1st> num, anal", num_1st, dpot(3,i)
         write(*,*) "zz> 1st> ", num_1st - dpot(3,i)
         write(*,*) "zz> 2nd> num,anal", num_2nd,diag_zz  
	 write(*,*) "zz> 2nd> diff", num_2nd - diag_zz 

      end do
c ---------------------------------------------------------------------
      write(*,*) "** testing off-diagonal term" 
      do i = 1, natom-1
         do j = i+1, natom
             indx = overlappair_idx(i,j)
             if ( indx .gt. 0 ) then
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xx term"
                tempx = coor(1, j)
                coor(1,j) = tempx + dx
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(1,j) = tempx - dx
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dx)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dx) 
                
                coor(1,j) = tempx
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+1)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "xxo 2nd> i,j", i, j
                write(*,*) "xxo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "xxo 2nd> diff:", num_2nd - anal_2nd 
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xy term"
                tempy = coor(2, j)
                coor(2,j) = tempy + dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(2,j) = tempy - dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(2,j) = tempy
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+2)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "xyo 2nd> i,j", i, j
                write(*,*) "xyo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "xyo 2nd> diff:", num_2nd - anal_2nd 
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dz)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dz) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+3)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "xzo 2nd> i,j", i, j
                write(*,*) "xzo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "xzo 2nd> diff:", num_2nd - anal_2nd 
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal yy term"
                tempy = coor(2, j)
                coor(2,j) = tempy + dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(2,i)
                 
                coor(2,j) = tempy - dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(2,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(2,j) = tempy
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+4)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "yyo 2nd> i,j", i, j
                write(*,*) "yyo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "yyo 2nd> diff:", num_2nd - anal_2nd 
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal yz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(2,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(2,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dz)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dz) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+5)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "yzo 2nd> i,j", i, j
                write(*,*) "yzo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "yzo 2nd> diff:", num_2nd - anal_2nd 
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal zz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(3,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(3,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+6)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
                write(*,*) "zzo 2nd> i,j", i, j
                write(*,*) "zzo 2nd> num,anal:", num_2nd, anal_2nd 
                write(*,*) "zzo 2nd> diff:", num_2nd - anal_2nd 
	     end if 
	 end do 
      end do
c --------------------------------------------------------------------- 

c      write(*,*) "DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
c      call egb_np2nd(natom)

999   return
      end

