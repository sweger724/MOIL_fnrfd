        subroutine e_field()

        implicit none
c
c This subroutine computes the contribution to the total
c energy and to the force of a constant (in space and time)
c external electric field
c Input
c 1) DV = potential energy (by default 0.001 V)
c 3) nfield = parallel to the electric field (by default 3, so z)
c 2) DVtype = switchin function between two types of potential profiles (by default 2)
c    DVtype = 1 ==> V(z)=-E*z, where E=-DV/zbox
c    DVtype = 2 ==> V(z)= E(0.5*zbox-z), if z>0
c                   V(z)=-E(0.5*zbox+z), if z<0 
c In a periodic boundary condition sense, these two potentials
c are exactly the same. (1) has a discontinuity at the extremes
c of the main box, (2) in the middle of the box. If the simulation
c involves a membrane in the middle of the box, type (2) is preferable because it puts
c the discontinuity where only uncharged particles are.
c                 ----o= =o----
c                 |   o= =o   |
c                 |   o= =o   |
c                 |   o= =o   |
c                 |   o= =o   |
c                 ----o= =o----
c               -zb/2   0    zb/2
c
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/EFIELD.BLOCK'

        integer i,namel
        double precision q,dx
        logical first
        data first/.TRUE./
        save first
        character*7 name
        namel = 7
        name = 'e_field'

        if (first) then
         if (efield.eq.0.d0 .and. (.not.esymyes)) then
          call alert(name,namel,'must specify efield if no symm',30,1)
         endif
         if (efbox.eq.0.d0 .and. (.not.esymyes)) then
          call alert(name,namel,'must specify efbox if no symm',29,1)
         endif
         if (esymyes) then
          if (nefield.eq.1) efbox = a
          if (nefield.eq.2) efbox = b
          if (nefield.eq.3) efbox = c
         endif
         if (efield.ne.0.d0) then
c convert efield from V/Angstrom to Kcal/mol/uc/Angstrom uc = unit charge
          efield = efield*23.06064815d0
          DV = -efield*efbox
         else if (DV.ne.0.d0) then
c convert delta V from V to Kcal/mol/uc
          DV = DV*23.06064815d0
          efield = -DV/efbox
         endif
         write (stdo,*) 'efield=',efield
         write (stdo,*) 'DV=',DV
         write (stdo,*) 'n=',nefield
         first = .FALSE.
        endif
c
c
        e_efield = 0.d0
        if (DVtype.eq.1) then
        do  i=1,npt
         if (flagchr(i)) then
           q = ptchg(i)
           dx = coor(nefield,i)
c dpot=-force
           dpot(nefield,i) = dpot(nefield,i) - q*efield
           e_efield = e_efield - q*efield*dx 
          endif
        enddo
        else
        do  i=1,npt
         if (flagchr(i)) then
          q = ptchg(i)
          if (coor(nefield,i).ge.0.d0) then
           dx = coor(nefield,i) - efbox*0.5d0
c dpot=-force
           dpot(nefield,i) = dpot(nefield,i) - q*efield
           e_efield = e_efield - q*efield*dx
          else
           dx = coor(nefield,i) + efbox*0.5d0
c dpot=-force 
           dpot(nefield,i) = dpot(nefield,i) - q*efield
           e_efield = e_efield - q*efield*dx
          endif
         endif
        enddo
        endif

        return
        end

