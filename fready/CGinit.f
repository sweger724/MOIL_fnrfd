      subroutine CGinit()
      
      implicit none
            
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c      local variables
      integer i,j,l,Nsum, p1,p2,p3,p4,m1,m2,m3,m4 
      integer id1,id2,getCGid, ibond, iang

      double precision x, E1, E2, E12, E3, dE
      double precision r,rx,ry,rz

      character*3 Rname(20)
      integer ir

c
c  Iniailize side chain masses
c
       call initCGmass()

       do i= 1,npt
         if (more(i).eq. 0.0) then
            Structure2nd(i)='C'
         else
            Structure2nd(i)='H'
         endif  
       enddo
        

C       Prepare backbone HB energy parameters 
        App(1) = 3.73d0 * 14.4d0
        App(2) = 0.d0   * 14.4d0
        App(3) = 5.13d0 * 14.4d0

       Bpp(1) = 1.306d0 * 14.4d0
       Bpp(2) = 1.129d0 * 14.4d0
       Bpp(3) = 0.335d0 * 14.4d0

       epp(1) = 0.305d0 * 14.4d0
       epp(2) = 0.365d0 * 14.4d0
       epp(3) = 0.574d0 * 14.4d0

       rpp(1) = 4.75d0
       rpp(2) = 4.54d0
       rpp(3) = 4.48d0


C
C  Prepare LJ datastructure
C
      do i =1, Namino
        do l=1, Namino
C       x^-12  + x^-6 for CMCM interactions
            LJr((i-1)*Namino+l,4)=1.d7 * 0.6
            LJr((i-1)*Namino+l,5)=-5.d3 * 0.6
C       x^-12 for CACM & CACA interactions
            LJr((i-1)*Namino+l+Namino**2,4)=1.d6 * 0.6
            LJr((i-1)*Namino+l+2*Namino**2,4)=1.d6 * 0.6
        end do
      end do


C Adjust LJ terms due that each parameter of LJ interaction is 
C readed only with precision of 3 digits... and thus we must here 
C extra ensure that lim_{x-> join+} LJ(x) = lim_{x-> join-} LJ(x)
C otherwise mini_pwl is complining
C we also ensure that LJ(13.5A) = 0 for each pair
      do i =1, 3*Namino**2
          call eCG_NB(LJa(i,10),i,E1,dE)
          call eCG_NB(LJa(i,10)+1d-6,i,E2,dE)
          LJa(i,0) = LJa(i,0) + E1 - E2
          
          call eCG_NB(CG_cutoff,i,E1,dE)
          LJa(i,0) = LJa(i,0) -E1
          LJr(i,3) = LJr(i,3) -E1
      enddo


      Rname(1) = "HIS"
      Rname(2) = "ARG"
      Rname(3) = "ASP"
      Rname(4) = "ASN"
      Rname(5) = "GLY"
      Rname(6) = "ALA"
      Rname(7) = "PRO"
      Rname(8) = "TYR"
      Rname(9) = "TRP"
      Rname(10) = "CYS"
      Rname(11) = "LEU"
      Rname(12) = "ILE"
      Rname(13) = "MET"
      Rname(14) = "VAL"
      Rname(15) = "LYS"
      Rname(16) = "GLU"
      Rname(17) = "GLN"
      Rname(18) = "SER"
      Rname(19) = "THR"
      Rname(20) = "PHE" 

c
c   Prepare bonds datastructure (for shake)
c
      do i = 1, nb
C       prepare cys - trans backbone types
C       
         if (ptnm(ib1(i)).eq.'CA  '.and.ptnm(ib2(i)).eq.'CA  ') then
            rx = coor(1,ib1(i)) - coor(1,ib2(i))
            ry = coor(2,ib1(i)) - coor(2,ib2(i))
            rz = coor(3,ib1(i)) - coor(3,ib2(i))
            r = dsqrt(rx**2 + ry**2 + rz**2)
            if (r .lt. 3.3d0) then
              bondType(i) = 21
            else
              bondType(i) = 0
            endif
         else
           bondType(i) = CGid(ib1(i))
         end if
C        write(6,*)"bondType(i): ", bondType(i)
      end do

C fill the list of S-S bonds (SSbond)
        Ncyscys = 0
        do i = 1,totmon-1
          do j = i+1,totmon
            if (moname(i).eq."CYSZ" .and. moname(j).eq."CYSZ") then
               rx = coor(1,poipt(i)) - coor(1,poipt(j))
               ry = coor(2,poipt(i)) - coor(2,poipt(j))
               rz = coor(3,poipt(i)) - coor(3,poipt(j))
               r = dsqrt(rx**2 + ry**2 + rz**2)
               if (r .lt. 2.9d0) then
C                 this is an S-S bond            
                  Ncyscys = Ncyscys + 1
C                  write(6,*)"SS:",Ncyscys,i,j,r
                  SSbond(Ncyscys,1) = i
                  SSbond(Ncyscys,2) = j
C                 add a bond to bonds data structure
                  nb = nb + 1
                  ib1(nb) = poipt(i)
                  ib2(nb) = poipt(j)
                  bondtype(nb) = 22
               endif
            end if
          enddo
        enddo


C  set minimum of potential to zero      
      do ibond = 0, 20
        if(CGBond(ibond,1) .eq. 1) then
            CGBond(ibond,3) = 0.d0
        elseif(CGBond(ibond,1) .eq. 2) then
           x = CGBond(ibond,2)
           call eCG_WellBondPot(x,ibond,E1,dE)
           x = CGBond(ibond,5)
           call eCG_WellBondPot(x,ibond,E2,dE)
           if (E1 .lt. E2) then
C             write(6,*)"Correct by:",ibond, E1
              CGBond(ibond,3) = CGBond(ibond,3) - E1
              CGBond(ibond,6) = CGBond(ibond,6) - E1
           else
C             write(6,*)"Correct by:",ibond, E2
              CGBond(ibond,3) = CGBond(ibond,3) - E2
              CGBond(ibond,6) = CGBond(ibond,6) - E2
           endif
        else
            x = CGBond(ibond,2)
            call eCG_WellBondPot(x,ibond,E1,dE)
            x = CGBond(ibond,5)
            call eCG_WellBondPot(x,ibond,E2,dE)
            x = CGBond(ibond,9)
            call eCG_WellBondPot(x,ibond,E3,dE)
            if (E1 .lt. E2) then
              E12 = E1
            else
              E12 = E2
            endif

            if (E12 .lt. E3) then
C             write(6,*)"Correct by:",ibond, E12            
              CGBond(ibond,3) = CGBond(ibond,3) - E12
              CGBond(ibond,6) = CGBond(ibond,6) - E12
              CGBond(ibond,10) = CGBond(ibond,10) - E12
           else
C             write(6,*)"Correct by:",ibond, E3    
              CGBond(ibond,3) = CGBond(ibond,3) - E3
              CGBond(ibond,6) = CGBond(ibond,6) - E3
              CGBond(ibond,10) = CGBond(ibond,10) - E3
           endif 
            
        endif
      end do

      write(6,*)"# of bonds: ",nb
c
c   Prepare angles datastructure
c
      do i = 1, nangl
         p1 = iangl1(i)
         p2 = iangl2(i)
         p3 = iangl3(i)
         m1 = poimon(p1)
         m2 = poimon(p2)
         m3 = poimon(p3)
         if ((ptnm(p1).eq."CM  ".and. m3.eq.m2+1) .or.
     &       (ptnm(p3).eq."CM  ".and. m1.eq.m2+1)) then
             anglType(i) = CGid(p2) + 1*Namino
         elseif((ptnm(p1).eq."CM  ".and. m3.eq.m2-1) .or.
     &          (ptnm(p3).eq."CM  ".and. m1.eq.m2-1)) then
             anglType(i) = CGid(p2) + 2*Namino
         else 
             anglType(i) = CGid(p2)
         endif
C        write(6,*)"Structure2nd:",m2,Structure2nd(p2),
C     &             anglType(i)
C        write(6,*)"anglType(i):",anglType(i),ptnm(p1),m1,m2,m3
      end do

C  set minimum of potential to zero
      do iang = 1, 3*Namino
C        write(6,*)"STored angles:",CGAngle(iang,1),CGAngle(iang,2),
C     & CGAngle(iang,3),CGAngle(iang,4),CGAngle(iang,5)
        if(CGAngle(iang,1) .eq. 1) then
            CGAngle(iang,3) = 0.d0
        elseif(CGAngle(iang,1) .eq. 2) then
           x = CGAngle(iang,2)
           call eCG_WellAnglePot(x,iang,E1,dE)
           x = CGAngle(iang,5)
           call eCG_WellAnglePot(x,iang,E2,dE)
           if (E1 .lt. E2) then
C             write(6,*)"Correct by:",iang, E1
              CGAngle(iang,3) = CGAngle(iang,3) - E1
              CGAngle(iang,6) = CGAngle(iang,6) - E1
           else
C             write(6,*)"Correct by:",iang, E2
              CGAngle(iang,3) = CGAngle(iang,3) - E2
              CGAngle(iang,6) = CGAngle(iang,6) - E2
           endif
        elseif (CGAngle(iang,1) .eq. 3) then 
            x = CGAngle(iang,2)
            call eCG_WellAnglePot(x,iang,E1,dE)
            x = CGAngle(iang,5)
            call eCG_WellAnglePot(x,iang,E2,dE)
            x = CGAngle(iang,9)
            call eCG_WellAnglePot(x,iang,E3,dE)
            if (E1 .lt. E2) then
              E12 = E1
            else
              E12 = E2
            endif

            if (E12 .lt. E3) then
C             write(6,*)"Correct by:",iang, E12
              CGAngle(iang,3) = CGAngle(iang,3) - E12
              CGAngle(iang,6) = CGAngle(iang,6) - E12
              CGAngle(iang,10) = CGAngle(iang,10) - E12
           else
C             write(6,*)"Correct by:",iang, E3
              CGAngle(iang,3) = CGAngle(iang,3) - E3
              CGAngle(iang,6) = CGAngle(iang,6) - E3
              CGAngle(iang,10) = CGAngle(iang,10) - E3
           endif

        endif
      end do

      write(6,*)"# of angles", nangl
     
     
c
c   Init torsions datastructure
c
      do i = 1, ntors
        p1 = itor1(i) 
        p2 = itor2(i)
        p3 = itor3(i)
        p4 = itor4(i)
        
        id1 = CGid(p2)
        id2 = CGid(p3)
        torType(i) = Namino*(id1-1)+id2
        if (ptnm(p1) .eq."CM  ") torType(i) = torType(i)+1*Namino**2
        if (ptnm(p4) .eq."CM  ") torType(i) = torType(i)+2*Namino**2 
C       write(6,*)"Torsion init: ",ptnm(p1),ptnm(p2),ptnm(p3),ptnm(p4)
C       write(6,*)"Torsion init: ",i,torType(i)
      end do
      write(6,*)"# of torsions ",ntors

      end
