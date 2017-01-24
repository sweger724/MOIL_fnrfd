        subroutine info(unit,jstep,nofreez,inofrz,tpo,tgroup,
     1          ntemp,curtemp)
        implicit none
c
c integer unit
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        character*4 name
        integer unit,namel,ntemp,i,j,l
        integer*4 jstep
        integer inofrz,nofreez(*),tpo(*),tgroup(*)
        double precision curtemp(*),vcm(3),totmass,ekin_cm
        double precision all

        name  = 'info'
        namel = 4

        if (my_pe.eq.0) then
           write(unit,*)
           write(unit,*)
           write(unit,*)
           write(unit,*)'-----------------------------------'
        endif

        if (debug) then
         write(*,*) ' in  INFO velo = '
         do i=1,npt
           write(*,*)'velo',i,'=',(velo(j,i),j=1,3)
         end do       
         write(unit,*)' ptms ',(ptms(i),i=1,npt)
        end if

        call multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     $          curtemp,ntemp)

        if (my_pe.eq.0) then
           write(unit,100)jstep
           write(unit,101)(curtemp(i),i=1,ntemp)
 100       format(1x,' At dynamics step ',i14)
 101       format(1x,' Current temperature(s) are  ',5(1x,f10.2))
        end if


        call wener(unit)

c add the contribution of the center of mass velocity to the energy 

        vcm(1) = 0.d0
        vcm(2) = 0.d0
        vcm(3) = 0.d0
        totmass = 0.d0

        do 15 i=1,inofrz
         l = nofreez(i)
         vcm(1) = vcm(1) + ptms(l)*velo(1,l)
         vcm(2) = vcm(2) + ptms(l)*velo(2,l)
         vcm(3) = vcm(3) + ptms(l)*velo(3,l)
         totmass = totmass + ptms(l)
15      continue

        vcm(1) = vcm(1)/totmass
        vcm(2) = vcm(2)/totmass
        vcm(3) = vcm(3)/totmass

        ekin_cm = 0.5d0*totmass*(vcm(1)*vcm(1)+vcm(2)*vcm(2)+
     1                           vcm(3)*vcm(3))
        ekin_cm = 0.d0

c        write (*,*) 'CM KIN ENERGY',ekin_cm
c        write (*,*) 'rest of kinetic energy',
c     1               0.5d0*kboltzmann*curtemp(1)*tgroup(1)
c        write (*,*) 'POT ENERGY',e_total

        all = 0.d0
        do i=1,ntemp
          all = all + 0.5d0*kboltzmann*curtemp(i)*tgroup(i)
        end do
        all = ekin_cm + all
        all = e_total + all

        if (my_pe.eq.0) then
           write(unit,102)all
 102       format(1x,' current energy (kinetic+potential) is ',f15.3)

           write(unit,*)'-----------------------------------'
           write(unit,*)
           write(unit,*)
           write(unit,*)
        endif

        return 
        end
