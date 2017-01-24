        program torstat
c 
c Torstat: A subroutine to calculate torsion statistics
c
        
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/CONVERT.BLOCK'

        character*4 name
        character*4 ctype
        integer namel
        integer of,geti
        logical find,empty,side,side2

        integer ucon,ucor,uwtors
        integer lp,lpend
        integer level,nptx,i,j
        integer rbin

        integer i1,i2,i3,i4
        integer iphi
        integer ires1,ires2
        double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
        double precision ux,uy,uz,vx,vy,vz,uv,uu,vv

        stdi = 5
        stdo = 6
        stderr = 0
        rbin = 1

        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        lpstr  = 1
        lpend  = 1

        jnkf = 25
        open (unit=jnkf,status='scratch')

        name  = 'Tors'
        namel = 4
        ucon  = -1
        ucor  = -1
        debug = .false.
        norew = .false.

        call init_var()

        ctype = 'unkw'
1       continue

        call rline(name,namel,stdi)

        if (find('debu')) debug = .true.
        if (find('norw')) norew = .true.

        if (find('file')) then
         if (find('rcon').or. find('conn')) then
                 ucon = of()
                 call rconn(ucon)
                 ires1 = 1
                 ires2 = totmon
         end if
         if (find('rcrd')) then
                 ucor = of()
                 rewind ucor
         end if
         if (find('tors')) uwtors = of()
        end if

        ires1 = geti('mon1',ires1)
        ires2 = geti('mon2',ires2)
        if (find('coor')) then
         call get4c(ctype,empty)
         if (ctype.eq.'CHAR') then
          lpstr = 1
          lpend = 1
         else if (ctype.eq.'PATH' .or. ctype.eq.'DYNA') then
          lpstr = geti('lpst',lpstr)
          lpend = geti('lpen',lpend)
         end if
        end if

        if (find('acti')) go to 2
        go to 1
2       continue
        
        iphi = -2
        do 5 i=ires1,ires2
         if (moname(i).eq.'PRO ' .or. moname(i).eq.'NTER'
     1          .or. moname(i).eq.'NTRG' .or. moname(i).eq.'CTER'
     2          .or. moname(i).eq.'CTRG' .or. moname(i).eq.'TIP3'
     3          .or. moname(i).eq.'HEME' .or. moname(i).eq.'HEM1'
     4          ) go to 5
c
c pick phi/psi/chi relevant atoms
c At iphi we put phi. At iphi+1 we put psi. At iphi+2 we put chi
c At iphi+3 chi2. If chi does not exist we put value of -1
c
          iphi = iphi + 4
          do 3 j=poipt(i-1)+1,poipt(i)
           if (ptnm(j).eq.'H   ') then
                itor1(iphi) = j
           else if (ptnm(j).eq.'N   ') then
                itor2(iphi) = j
                itor1(iphi+1) = j
                itor1(iphi+2) = j
           else if (ptnm(j).eq.'CA  ') then
                itor3(iphi) = j
                itor2(iphi+1) = j
                itor2(iphi+2) = j
                itor1(iphi+3) = j
           else if (ptnm(j).eq.'C   ') then
                itor4(iphi) = j
                itor3(iphi+1) = j
           else if (ptnm(j).eq.'O   ') then
                itor4(iphi+1) = j
          end if
3        continue
         side = .not. ((moname(i).eq.'GLY ').or.(moname(i).eq.'ALA '))
         side2 = side .and. (.not. ((moname(i).eq.'VAL ') .or.
     1          (moname(i).eq.'THR ') .or. (moname(i).eq.'SER ')))
         if (.not.side) then
                itor3(iphi+2) = -1
         else
          do 4 j=poipt(i-1)+1,poipt(i)
           if (ptnm(j).eq.'CB  ') then
                itor3(iphi+2) = j
           else if (ptnm(j)(2:2).eq.'G') then
                itor4(iphi+2) = j
                go to 45
           end if
4         continue
45        continue
         end if 
         if (.not.side2) then
                 itor2(iphi+3) = -1
         else
          do 41 j=poipt(i-1)+1,poipt(i)
           if (ptnm(j).eq.'CB  ') then
                itor2(iphi+3) = j
           else if (ptnm(j)(2:2).eq.'G') then
                itor3(iphi+3) = j
           else if (ptnm(j)(2:3).eq.'D1') then
                itor4(iphi+3) = j
                go to 415
           end if
41         continue
415        continue
         end if 

         if (side .and. (itor3(iphi+2).eq.-1)) then
          level = 1
          call alert(name,namel,'Unknown side chain',18,level)
         end if
5       continue

        if (ctype.eq.'unkw') then
         level = 1
         call alert(name,namel,'Missing coor. type',18,level)
        end if

        if (lpstr.eq.0 .or. lpend.eq.0) then
         level = 1
         call alert(name,namel,'Illegal range of structures',27,level)
        end if

        nptx = npt
        do 12 lp=lpstr,lpend
         if (ctype.eq.'CHAR') then
          call getcrd(ucor,'CHARM')
         else if (ctype.eq.'PATH') then
          rewind ucor
          call rpath(ucor,lp)
         else if (ctype.eq.'DYNA') then
          if (.not.norew) rewind ucor
          call rdyncrd(ucor,lp,nptx,0,rbin)
         end if

c make a run on phi psi and chi1 torsions
c This is while working through (first the wcon file). Each residue
c which is not first or last must have phi and psi. Each residue
c which is gly, ala or proline must have a chi1 torsion
c

         do 6 i=1,4*iphi
c PHI, PSI, CHI1 and CHI2
          if (itor3(i).lt.0) then
                phase1(i) = -999.d0
                phase2(i) = -999.d0
                phase3(i) = -999.d0
                go to 6
          else if (itor2(i).lt.0) then
                phase1(i) = -999.d0
                phase2(i) = -999.d0
                phase3(i) = -999.d0
                go to 6
          end if
          i1 = itor1(i)
          i2 = itor2(i)
          i3 = itor3(i)
          i4 = itor4(i)

          dx1 = coor(1,i2) - coor(1,i1)
          dy1 = coor(2,i2) - coor(2,i1)
          dz1 = coor(3,i2) - coor(3,i1)

          dx2 = coor(1,i3) - coor(1,i2)
          dy2 = coor(2,i3) - coor(2,i2)
          dz2 = coor(3,i3) - coor(3,i2)

          dx3 = coor(1,i4) - coor(1,i3)
          dy3 = coor(2,i4) - coor(2,i3)
          dz3 = coor(3,i4) - coor(3,i3)

          ux  = dy1*dz2 - dz1*dy2
          uy  = dz1*dx2 - dx1*dz2
          uz  = dx1*dy2 - dy1*dx2

          vx  = dy2*dz3 - dz2*dy3
          vy  = dz2*dx3 - dx2*dz3
          vz  = dx2*dy3 - dy2*dx3

          uu  = (ux*ux+uy*uy+uz*uz)
          vv  = (vx*vx+vy*vy+vz*vz)
          uv  = (ux*vx+uy*vy+uz*vz)/dsqrt(uu*vv)

          phase1(i) = dacos(uv)*pi180

          dx1 = uy*vz - uz*vy
          dy1 = uz*vx - ux*vz
          dz1 = ux*vy - uy*vx

          if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) phase1(i) = -phase1(i)
          phase1(i) = phase1(i) - 180.d0
          if (phase1(i).lt.-180) phase1(i) = phase1(i) + 360
          phase2(i) = phase1(i)
          phase3(i) = phase1(i)

6        continue


         iphi = -2
         do 7 i=ires1,ires2
          if (moname(i).eq.'PRO ' .or. moname(i).eq.'NTER'
     1          .or. moname(i).eq.'NTRG' .or. moname(i).eq.'CTER'
     2          .or. moname(i).eq.'CTRG' .or. moname(i).eq.'TIP3'
     3          .or. moname(i).eq.'HEME' .or. moname(i).eq.'HEM1'
     4          ) go to 7
          iphi = iphi + 4
           write(uwtors,101) i,moname(i),phase1(iphi),phase1(iphi+1)
     1  ,phase1(iphi+2),phase1(iphi+3),lesid(poipt(i))
101        format(1x,i7,1x,a4,1x,4(f8.3,1x),1x,i5)
7         continue
12       continue
          stop
          end
