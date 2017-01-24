        subroutine symcdie_ewald()
        implicit none

c calculating van der Waals and electorstatic direct energies for symmetry 
c related particles (impart contribution to direct sums) within PME algorithm 

      include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/PRESS.BLOCK'

        double precision pick 
        double precision epstmp,spi,ovt,del,term
        double precision xi,yi,zi,rij,xerfc
        double precision rx,ry,rz,r2,s2,aa,bb,b1,e1,e2,qi,df,df1,
     1   df2,ai,bi
        double precision s,s6,derfc,dx,erfcc,aa1,bb1

        double precision rmolx,rmoly,rmolz

        integer istart,i,j,k,jbeg,jend,symidx,ireal,ind

        logical first

        data first/.true./

        save first
        save epstmp,pick,spi,term,ovt,del

        if (first ) then
                if (eelyes) then
                        epstmp = kofdie/eps
                else
                        epstmp =  0.0d0
                end if
        
                if (evdyes) then
                        pick = 1.0d0
                else
                        pick = 0.0d0
                end if
c
                spi=1.0d0/pi
                spi=dsqrt(spi)
                term=2.0d0*spi*ewaldcof
                ovt=1.0d0/3.0d0
                del=1.0d0/erftbdns
                first = .false.
        end if


c
c This first loop is up to the short (vdW) radius

        istart = 1

        if (prll_on_off) then
           jbeg = point1(dpoipt(monp(my_pe+1)))+1
        else   
           jbeg   = point1(npt-1) + 1
        endif   
c 102 is a loop on all symmetry operatio

        do 102 symidx = 1,symnum
c 101 is a loop on all atoms that belong to the current symmetry operation
         do 101 i=istart,iblock1(symidx)

                jend = psym1(i)


                ireal = symreal1(i)
                if (jbeg.le.jend) then
                   if (rotop(symidx)) then
                        xi = arot*coor(1,ireal) + symop(1,symidx)*a
                        yi = brot*coor(2,ireal) + symop(2,symidx)*b
                        zi = crot*coor(3,ireal) + symop(3,symidx)*c
                   else
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        yi = coor(2,ireal) + symop(2,symidx)*b
                        zi = coor(3,ireal) + symop(3,symidx)*c
                   end if

c the symmetry related forces for infinite lattice are only half their
c "real" value
c therefore we multiply by 0.5 below. Note that for finite systems this
c division
c should be omitted




                  
                        ai=epsgm12(ireal)*pick
                        bi=epsgm6(ireal)*pick
                        
                 
                        
                        qi=ptchg(ireal)*epstmp


c 100 is a loop on the neighbors (j) to the atom (i)
                        do 100 k=jbeg,jend
                                j=list1(k)
                                rx = xi - coor(1,j)
                                ry = yi - coor(2,j)
                                rz = zi - coor(3,j)
                                r2=rx*rx+ry*ry+rz*rz
                                s2=1.0d0/r2
                                rij=dsqrt(r2)
                                if (r2.gt.cutvdw2 .and.
     1                                  (.not.nocut)) then
                                        df1 = 0.d0
                                        go to 50
                                end if
                                s6=s2*s2*s2

c       ileana

                                if(.not.arith)then

                                aa=ai*epsgm12(j)*s6*s6
                                bb=bi*epsgm6(j)*s6

                                else
                                aa1 = ai*epsgm12(j)
                                bb1 = 0.5d0*(bi + epsgm6(j))
                                bb1 = bb1*bb1*bb1
                                bb1 = bb1 * bb1
                                bb=aa1*bb1*s6
                                bb1 = bb1*bb1
                                aa = aa1*bb1*s6*s6

                                endif
c       ileana
                                e1 = aa - bb

c                               write(*,*) "1st e1 in symcdie_ewald is ", e1
                                e_vsym = e_vsym + e1 
                                df1 = -6.0d0*s2*(aa+e1)

c pressure
c virial calculation
       if (pressON) then
        virial = virial +
     1  (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*df1*rx +
     2  (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*df1*ry +
     3  (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*df1*rz
             rmolx = xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a
             rmoly = ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b
             rmolz = zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c
             virXX = virXX + df1*rx*rmolx
             virYY = virYY + df1*ry*rmoly
             virZZ = virZZ + df1*rz*rmolz
                  !write (333,*) 'i j dfvdw drmol'
                  !write (333,*) j,ireal,df1*rx,rmolx
                  !write (333,*) j,ireal,df1*ry,rmoly
                  !write (333,*) j,ireal,df1*rz,rmolz
       endif
cccccccccccccccccc

50                              continue
                                s  = 1.0d0/rij
                                xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                                ind = int(erftbdns*xerfc) + 1
                                dx = xerfc - (ind-1)*del
                                derfc = -erf_arr(2,ind) - dx*(  
     $                           erf_arr(3,ind)+0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                                e2 = qi*ptchg(j)*s*erfcc

ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = qi*ptchg(j)*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)

             rmolx = xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a
             rmoly = ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b 
             rmolz = zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c

                  !write (333,*) 'i j dfdir drmol'
                  !write (333,*) j,ireal,V_PIdir*rx,rmolx
                  !write (333,*) j,ireal,V_PIdir*ry,rmoly
                  !write (333,*) j,ireal,V_PIdir*rz,rmolz

c in what follows we need the minus sign as rx is x(i)-x(j), while
c we need x(j)-x(i)

c xx
                  V_PIdirXX = V_PIdirXX + V_PIdir*(-rx)*rmolx
c yy
                  V_PIdirYY = V_PIdirYY + V_PIdir*(-ry)*rmoly
c zz              
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*(-rz)*rmolz
                  endif
cccccccccc

c       ileana
c                               write(*,*) "e2 in symcdie_ewald is  ", e2

                                df2 = -qi*ptchg(j)*s*(s2*erfcc +
     *                            s*ewaldcof*derfc)

                                df = df1 + df2
                                rx = df*rx
                                ry = df*ry
                                rz = df*rz
                                if (.not.rotop(symidx)) then
                                 dpot(1,ireal) = dpot(1,ireal) + rx
                                 dpot(2,ireal) = dpot(2,ireal) + ry
                                 dpot(3,ireal) = dpot(3,ireal) + rz
                                end if
                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_lsym = e_lsym + e2
100                     continue
                end if
        jbeg   = jend + 1

101     continue

        istart = iblock1(symidx) + 1
102     continue


c
c This second loop deals with vdW of uncharged particles


        if (prll_on_off) then
           jbeg = point2(dpoipt(monp(my_pe+1)))+1
        else   
           jbeg   = point2(npt-1) + 1
        endif 


        istart = 1
c 202 is a loop on all symmetry operation
        do 202 symidx = 1,symnum
c 201 is a loop on all atoms that belong to the current symmetry operation
         do 201 i=istart,iblock2(symidx)

                jend = psym2(i)

                ireal = symreal2(i)
                if (jbeg.le.jend) then
                  if (rotop(symidx)) then
                        xi = arot*coor(1,ireal) + symop(1,symidx)*a
                        yi = brot*coor(2,ireal) + symop(2,symidx)*b
                        zi = crot*coor(3,ireal) + symop(3,symidx)*c
                  else
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        yi = coor(2,ireal) + symop(2,symidx)*b
                        zi = coor(3,ireal) + symop(3,symidx)*c
                  end if
c the symmetry related forces for infinite lattice are only half
c their "real" value
c therefore we multiply by 0.5 below. Note that for finite systems
c this division
c should be omitted

c       ileana

                  if(.not.arith) then
                        ai=epsgm12(ireal)*pick
                        bi=epsgm6(ireal)*pick

                        else

                        ai=epsgm12(ireal)*pick
                        bi=epsgm6(ireal)

                        endif


c 200 is a loop on the neighbors (j) to the atom (i)
                        do 200 k=jbeg,jend
                                j=list2(k)
                                rx = xi - coor(1,j)
                                ry = yi - coor(2,j)
                                rz = zi - coor(3,j)
                                r2=rx*rx+ry*ry+rz*rz
                                if (r2.gt.cutvdw2 .and.
     1                                  (.not.nocut)) go to 200
                                s2=1.0d0/r2
                                s6=s2*s2*s2

c       ileana

                                if(.not.arith)then

                                aa=ai*epsgm12(j)*s6*s6
                                bb=bi*epsgm6(j)*s6

                                else
                                aa1 = ai*epsgm12(j)
                                bb1 = 0.5d0*(bi + epsgm6(j))
                                bb1 = bb1*bb1*bb1
                                bb1 = bb1 * bb1
                                bb=aa1*bb1*s6
                                bb1 = bb1*bb1
                                aa = aa1*bb1*s6*s6

                                endif
        
c       ileana
                                e1 = aa - bb

c                               write(*,*)"2nd loop e1 in symcdie_ewald is ",e1

                                df = -6.0d0*s2*(aa+e1)

c pressure
c virial calculation
       if (pressON) then
        virial = virial +
     1  (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*df*rx +
     2  (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*df*ry +
     3  (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*df*rz
             rmolx = xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a
             rmoly = ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b
             rmolz = zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c
             virXX = virXX + df*rx*rmolx
             virYY = virYY + df*ry*rmoly
             virZZ = virZZ + df*rz*rmolz
                  !write (333,*) 'i j dfvdw drmol'
                  !write (333,*) j,ireal,df*rx,rmolx
                  !write (333,*) j,ireal,df*ry,rmoly
                  !write (333,*) j,ireal,df*rz,rmolz
       endif
cccccccccccccccccc

                                rx = df*rx
                                ry = df*ry
                                rz = df*rz
                                if (.not.rotop(symidx)) then
                                 dpot(1,ireal) = dpot(1,ireal) + rx
                                 dpot(2,ireal) = dpot(2,ireal) + ry
                                 dpot(3,ireal) = dpot(3,ireal) + rz
                                end if
                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_vsym = e_vsym + e1 
200                     continue
                end if
        jbeg   = jend + 1
201     continue
        istart = iblock2(symidx) + 1
202     continue


c This third loop is on longer range electrostatic

        istart = 1
        if (prll_on_off) then
           jbeg = point3(dpoipt(monp(my_pe+1)))+1
        else   
           jbeg   = point3(npt-1) + 1
        endif 

c 302 is a loop on all symmetry operation
        do 302 symidx = 1,symnum
         if (symidx.le.4 .or. (.not.metalyes)) then
c 301 is a loop on all atoms that belong to the current symmetry operation
         do 301 i=istart,iblock3(symidx)

                jend = psym3(i)
                ireal = symreal3(i)
                if (jbeg.le.jend) then
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        yi = coor(2,ireal) + symop(2,symidx)*b
                        zi = coor(3,ireal) + symop(3,symidx)*c
c the symmetry related forces for infinite lattice are only half
c their "real" value
c therefore we multiply by 0.5 below. Note that for finite systems
c this division
c should be omitted
                        qi=ptchg(ireal)*epstmp


c 100 is a loop on the neighbors (j) to the atom (i)
                        do 300 k=jbeg,jend
                                j=list3(k)
                                rx = xi - coor(1,j)
                                ry = yi - coor(2,j)
                                rz = zi - coor(3,j)
                                r2=rx*rx+ry*ry+rz*rz
                                if (r2.gt.cutele2 .and.
     1                                  (.not.nocut)) go to 300
                                rij=dsqrt(r2)
                                s2 = 1.0d0/r2
                                s  = 1.0d0/rij
                                xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                                ind = int(erftbdns*xerfc) + 1
                                dx = xerfc - (ind-1)*del
                                derfc = -erf_arr(2,ind) - dx*(  
     $                           erf_arr(3,ind)+0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                                
                                e2 = qi*ptchg(j)*s*erfcc
                                df = -qi*ptchg(j)*s*(s2*erfcc +
     *                            s*ewaldcof*derfc)

ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = qi*ptchg(j)*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)

             rmolx = xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a
             rmoly = ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b  
             rmolz = zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c

                  !write (333,*) 'i j dfdir drmol'
                  !write (333,*) j,ireal,V_PIdir*rx,rmolx
                  !write (333,*) j,ireal,V_PIdir*ry,rmoly
                  !write (333,*) j,ireal,V_PIdir*rz,rmolz

c in what follows we need the minus sign as rx is x(i)-x(j), while
c we need x(j)-x(i)

c xx
                  V_PIdirXX = V_PIdirXX + V_PIdir*(-rx)*rmolx
c yy
                  V_PIdirYY = V_PIdirYY + V_PIdir*(-ry)*rmoly
c zz              
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*(-rz)*rmolz
                  endif

cccccccccccccccccc


c non-interpolated erfc
c                               call erfcfun(xerfc,erfc)
c                               e2 = qi*ptchg(j)*s*erfc
c                               df = - qi*ptchg(j)*s* (s2*erfc + 
c     *                            s*term*exp(-xerfc**2))


                                rx = df*rx
                                ry = df*ry
                                rz = df*rz

                                dpot(1,ireal) = dpot(1,ireal) + rx
                                dpot(2,ireal) = dpot(2,ireal) + ry
                                dpot(3,ireal) = dpot(3,ireal) + rz

                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_lsym = e_lsym + e2
300                     continue
                end if
        jbeg   = jend + 1

301     continue

c Here comes the metal image charges part
c Note cutoff MUST be smaller than half of the box size.
c
        else

         do 401 i=istart,iblock3(symidx)

                jend = psym3(i)
                ireal = symreal3(i)
                if (jbeg.le.jend) then
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        zi = coor(3,ireal) + symop(3,symidx)*c
                        if (coor(2,ireal).lt.0.d0) then
                                yi = -b - coor(2,ireal)
                        else
                                yi = b - coor(2,ireal)
                        end if
c
c an image charge has the reverse sign
c
                        qi=-ptchg(ireal)*epstmp


c 100 is a loop on the neighbors (j) to the atom (i)
                        do 400 k=jbeg,jend
                                j=list3(k)
                                rx = xi - coor(1,j)
                                ry = yi - coor(2,j)
                                rz = zi - coor(3,j)
                                r2=rx*rx+ry*ry+rz*rz
                                if (r2.gt.cutele2 .and.
     1                                  (.not.nocut)) go to 400
                                rij=dsqrt(r2)
                                s2 = 1.0d0/r2
                                s  = 1.0d0/rij
                                xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                                ind = int(erftbdns*xerfc) + 1
                                dx = xerfc - (ind-1)*del
                                derfc = -erf_arr(2,ind) - dx*(  
     $                           erf_arr(3,ind)+0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                                e2 = qi*ptchg(j)*s*erfcc
                                df = -qi*ptchg(j)*s*(s2*erfcc +
     *                            s*ewaldcof*derfc)

c non-interpolated erfc
c                               call erfcfun(xerfc,erfc)
c                               e2 = qi*ptchg(j)*s*erfc
c                               df = - qi*ptchg(j)*s* (s2*erfc + 
c     *                            s*term*exp(-xerfc**2))


                                rx = df*rx
                                ry = df*ry
                                rz = df*rz
c
c no forces on image particles in
c contrast to usual symmetry operation
                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_lmet = e_lmet + e2
400                     continue
                end if
        jbeg   = jend + 1

401     continue

c end of metal branch
        end if
        istart = iblock3(symidx) + 1
302     continue


        return 
        end
