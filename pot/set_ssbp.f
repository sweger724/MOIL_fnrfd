C23456789012345678901234567890123456789012345678901234567890123456789012
       SUBROUTINE SET_SSBP()
       
       IMPLICIT NONE
       
       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/SSBP.BLOCK'  

       integer stdo,stdi 

       stdi = 5
       stdo = 6

      LMA    = NMULT+1
      LMA2   = NMULT*2+1

      if (qkirk) then
      call facto(lma2)
      endif

      if (qcavi) then
      write(stdo,'(2a)')'Cavity potential applied to atom selection',
     &                  'only (based on RISM-HNC with water)'
c Parameters for the cavity potential
      ACAV(5)=-1.6649500d0
      ACAV(1)= 0.56198800d0
      ACAV(2)=-0.072798148d0
      ACAV(3)= 0.00426122036d0
      ACAV(4)=-0.0000925233817d0
      ACAV(6)= 0.0840D0
      ACAV(7)=15.39333D0
      BCAV(1)= 1.319978287d0
      BCAV(2)=-0.840953501d0
      BCAV(3)=-0.001602388122d0
      BCAV(4)=-8.392886499d0
      BCAV(5)=BCAV(2)+BCAV(4)
      BCAV(6)= 1.6D0
      BCAV(7)=-8.4751210228D0
      endif

      if (qhsr) then
      write(stdo,'(A)')' Hard-Sphere-Restriction contribution'
      WRITE(stdo,'(2(A,F10.5))') ' Pressure        ', PRESI,
     &                           ' Surface tension ', STENS
      endif

      if (qangu) then
      write(stdo,'(2A)') ' Angular potential for H2O (for ',
     &                   ' isotropic orientation near the boundary)'
c Parameters for the angular potential
      DRHA   = -1.D0
      CANG(5)=0.840661d0
      CANG(4)=- 1.20064327666d0
      CANG(3)=- 3.067139576698d0
      CANG(2)=1.766839115555d0
      CANG(1)= 2.4085919002d0
      endif

       RETURN
       END
