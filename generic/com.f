      subroutine com (d, n, coor, ptms, cm)

c     IN
c     d: dimensionality
c     n: number of particles
c     coor: [d n]: coordinates

c     OUT
c     cm: [d 1]: center of mass


      integer d, di, n, ni
      double precision cm(d),ptms(n)
      double precision coor(d,n)
      double precision mass


      mass = 0.d0
        do 300 ni = 1,n
                mass = mass+ ptms(ni)
300     continue
        mass = 1.d0/mass
      do 100 di = 1, d
         cm(di) = 0.d0
         do 200 ni = 1, n
            cm(di) = cm(di) + coor(di, ni)*ptms(ni)
 200     continue
 100  continue

        do 400 di = 1,d
                cm(di) = cm(di)*mass
400     continue
      end
