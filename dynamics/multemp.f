        subroutine multemp(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
c
c calculate values of multiple temperatures: groups are defined by tpo
c vx vy vz - velocities
c natom    - number of particles (dimension of the vectors)
c inofrz   - number of particles that are not freezed
c nofreez  - vector pointer to the position of the particles
c
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONVERT.BLOCK'

        double precision velo(3,*)
        double precision ptms(*),curtemp(*),vcm(3),totmass
        integer inofrz,ntemp,nofreez(*),tpo(*),tgroup(*)
c
c local
c
        integer i,l,n

c@
c	write(*,*)' Multemp: inofrz = ',inofrz



        do 1 n=1,ntemp
         curtemp(n) = 0.d0
1       continue

c modification to remove the contribution of the center of mass
c to the temperature 

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

        vcm(1) = 0.d0
        vcm(2) = 0.d0
        vcm(3) = 0.d0

        do 2 i=1,inofrz
          l = nofreez(i)
          n = tpo(l)
          curtemp(n) = curtemp(n) + ptms(l) *
     1           ((velo(1,l)-vcm(1))*(velo(1,l)-vcm(1))
     1          + (velo(2,l)-vcm(2))*(velo(2,l)-vcm(2))
     3          + (velo(3,l)-vcm(3))*(velo(3,l)-vcm(3)))
2       continue

        
        do 3 n=1,ntemp
           curtemp(n)=curtemp(n)/(kboltzmann*tgroup(n))
           !write(6,*)"EEE",curtemp(n)
3       continue
        return
        end

        subroutine multemp_prll(velo,ptms,nofreez,inofrz,tpo,tgroup,
     1          curtemp,ntemp)
c
c calculate values of multiple temperatures: groups are defined by tpo
c vx vy vz - velocities
c natom    - number of particles (dimension of the vectors)
c inofrz   - number of particles that are not freezed
c nofreez  - vector pointer to the position of the particles
c
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/CONVERT.BLOCK'

        double precision velo(3,*)
        double precision ptms(*),curtemp(*)
        integer inofrz,ntemp,nofreez(*),tpo(*),tgroup(*)
c
c local
c
        integer l,n

        do 1 n=1,ntemp
         curtemp(n) = 0.d0
1       continue
        do 2 l=pt_start,pt_end
c         l = nofreez(i)
          n = tpo(l)
          curtemp(n) = curtemp(n) + ptms(l) * (velo(1,l)*velo(1,l)
     1          + velo(2,l)*velo(2,l) + velo(3,l)*velo(3,l))
2       continue

        do 3 n=1,ntemp
           call reduce_1(curtemp(n))
           curtemp(n)=curtemp(n)/(kboltzmann*tgroup(n))
3       continue
        return
        end

