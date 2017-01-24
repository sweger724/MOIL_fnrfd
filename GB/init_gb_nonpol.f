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
