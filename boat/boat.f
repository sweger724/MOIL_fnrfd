        program boat
c 
c BOAT: A program to calculate the BOnd, Angles and Torsions of one
c               or a series of structures.
c The calculated values are stored in the CONNECT array and are
c written over the equilibrium values. Therefore do not combine this
c routine with energy calculations
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
        logical find,empty

        integer ucon,ucor,uwboat
        integer lp,lpend
        integer level,nptx,i
        integer rbin

        integer i1,i2,i3,i4
        double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
        double precision ux,uy,uz,vx,vy,vz,uv,uu,vv
        double precision r12,r23,cost

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

        name  = 'boat'
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
         if (find('rcon').or.find('conn')) ucon = of()
         if (find('rcrd')) ucor = of()
         if (find('boat')) uwboat = of()
        end if

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
        
        call rconn(ucon)

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

         do 3 i=1,nb
          i1 = ib1(i)
          i2 = ib2(i)
          dx1 = coor(1,i2) - coor(1,i1)
          dy1 = coor(2,i2) - coor(2,i1)
          dz1 = coor(3,i2) - coor(3,i1)
          kbond(i) = dsqrt(dx1*dx1+dy1*dy1+dz1*dz1)
          ktors1(i) = kbond(i) - req(i)
3        continue

         do 4 i=1,nangl
          i1 = iangl1(i)
          i2 = iangl2(i)
          i3 = iangl3(i)

          dx1 = coor(1,i2) - coor(1,i1)
          dy1 = coor(2,i2) - coor(2,i1)
          dz1 = coor(3,i2) - coor(3,i1)

          dx2 = coor(1,i2) - coor(1,i3)
          dy2 = coor(2,i2) - coor(2,i3)
          dz2 = coor(3,i2) - coor(3,i3)

          r12 = dsqrt(dx1*dx1+dy1*dy1+dz1*dz1)
          r23 = dsqrt(dx2*dx2+dy2*dy2+dz2*dz2)

          cost = (dx1*dx2+dy1*dy2+dz1*dz2)/(r12*r23)

          kangl(i) = dacos(cost)*pi180
          ktors2(i) = kangl(i) - angleq(i)*pi180
4        continue

         do 5 i=1,ntors
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

          phase2(i) = phase1(i)
          phase3(i) = phase1(i)

5        continue


         do 6 i=1,nimp
          i1 = iimp1(i)
          i2 = iimp2(i)
          i3 = iimp3(i)
          i4 = iimp4(i)

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

          kimp(i) = dacos(uv)*pi180

          dx1 = uy*vz - uz*vy
          dy1 = uz*vx - ux*vz
          dz1 = ux*vy - uy*vx

          if (dx1*dx2+dy1*dy2+dz1*dz2 .lt. 0) kimp(i) = - kimp(i)

          ktors3(i) = kimp(i) - impeq(i)*pi180

6        continue
         
         write(uwboat,100)lp
100      format(//,1x,
     1    ' BOnd Angle Torsion (BOAT) list for structure # ',i5,//,
     2    ' For the bond, angle, improper torsion we provide the ',
     2    'atomic indices,',/,' current value and ',
     3    'deviation from equilibrium. Deviation not given for ',
     4    'torsion ')
         write(uwboat,101)
101      format(/,1x,'**** BONDS')
         do 8 i=1,nb
          write(uwboat,102)ib1(i),ib2(i),kbond(i),ktors1(i)
102       format(1x,i7,1x,i7,2(2x,f10.5))
8        continue

         write(uwboat,103)
103      format(/,1x,'**** ANGLES')
         do 9 i=1,nangl
          write(uwboat,104)iangl1(i),iangl2(i),iangl3(i),
     1           kangl(i),ktors2(i)
104       format(3(1x,i7),2(2x,f10.5))
9        continue

         write(uwboat,105)
105      format(/,1x,'**** TORSIONS')
         do 10 i=1,ntors
          write(uwboat,106)itor1(i),itor2(i),itor3(i),itor4(i),
     1                     phase1(i),phase2(i),phase3(i)
106       format(4(1x,i7),2x,3(f10.5))
10       continue

         write(uwboat,107)
107      format(/,1x,'**** IMPROPER TORSIONS')
         do 11 i=1,nimp
          write(uwboat,108)iimp1(i),iimp2(i),iimp3(i),
     1            iimp4(i),kimp(i),ktors3(i)
108       format(4(1x,i7),2(2x,f10.5))
11       continue

12      continue

        write(stdo,109)
109     format(//,1x,' **** BOAT COMPLETED **** ')
        stop
        end
