        Subroutine symwat_ewald()
        implicit none

C       This subrouitne calculates the impart water contributions to
C       the direct sum within PME algorithm 

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'

        double precision rx,ry,rz,r2,s2,aoo,boo,e1,e2,df
        double precision s,s6,spi,rij,xerfc,erfcc
        double precision xio,yio,zio,dxio,dyio,dzio,dx,derfc
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
        double precision qo,qh,qoh,term,ovt,del


        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,jh1,jh2
        integer symidx,iwater,istart,ind
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,ohchrg,oochrg,hhchrg
        save spi,qo,qh,qoh,term,ovt,del

        if (first) then
                if (moname(realmono(nwaters)).eq.'TIP3') then
                        oosgm6 = 595.05497170775496395324d0
                        oosgm12 = 582002.66166028445200793331d0
                        ohchrg  = -115.4871969048d0
                        oochrg  = 230.9743938096d0
                        hhchrg  = 57.7435984524d0
                else if (moname(realmono(nwaters)).eq.'SPCE') then
                        oosgm6  = 625.50000300044132791676
                        oosgm12 = 629400.00388394183256163044
                        ohchrg  = -119.2843958022
                        oochrg  = 238.5687916044
                        hhchrg  = 59.642197901104
                else 
                        write(*,*)' Water was not identified '
                        stop
                end if
                spi=1.0d0/pi
                spi=dsqrt(spi)
                qo=oochrg
                qh=hhchrg
                qoh=ohchrg
                term=2.0d0*spi*ewaldcof
                ovt=1.0d0/3.0d0
                del=1.0d0/erftbdns
                first = .false.
        end if
c
c
c To speed up water calculations, we do not support for
c TIP3-TIP3 interactions other dielectric constant than 1
c or applying LES. If you need LES water create a water monomer
c with a name different from TIP3.
c


        istart= 1
        if (prll_on_off) then
                jbeg1 = poinwat1(my_nwat+1)+1
        else
                jbeg1 = poinwat1(nwaters) + 1
        end if
c loop on symmetry operations
c
        do 400 symidx = 1,symnum
         do 300 iwater=istart,iblckwt1(symidx)
                jend1  = psymwt1(iwater)
                i      = symmwat1(iwater)


                io = dpoipt(i)-2
                ih1=io+1
                ih2=io+2

                xio  = coor(1,io) + symop(1,symidx)*a
                yio  = coor(2,io) + symop(2,symidx)*b
                zio  = coor(3,io) + symop(3,symidx)*c

                xih1 = coor(1,ih1) + symop(1,symidx)*a
                yih1 = coor(2,ih1) + symop(2,symidx)*b
                zih1 = coor(3,ih1) + symop(3,symidx)*c

                xih2 = coor(1,ih2) + symop(1,symidx)*a
                yih2 = coor(2,ih2) + symop(2,symidx)*b
                zih2 = coor(3,ih2) + symop(3,symidx)*c

                dxio = 0.d0
                dyio = 0.d0
                dzio = 0.d0

                dxih1 = 0.d0
                dyih1 = 0.d0
                dzih1 = 0.d0

                dxih2 = 0.d0
                dyih2 = 0.d0
                dzih2 = 0.d0

                
                if (jbeg1.le.jend1) then

                do 200 k=jbeg1,jend1
                        j  = listwt1(k)
                        jo = dpoipt(j)-2
                        jh1= jo + 1
                        jh2= jo + 2

c take care of oxygen - oxygen van der Waals
                        rx = xio - coor(1,jo)
                        ry = yio - coor(2,jo)
                        rz = zio - coor(3,jo)

                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2

                        if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
                                df  = 0.d0
                                go to 100
                        end if

                        s6=s2*s2*s2
                        aoo = oosgm12*s6*s6
                        boo = oosgm6*s6
                        e1     = aoo - boo
                        df     = -6.0d0*s2*(aoo+e1)
                        e_vsym = e_vsym + e1 


c the electrostatic part should be done without buffering
100                     continue

c Oxygen - Oxygen part
                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = oochrg*s*erfcc

cc Luca for virial (vdW only: take only df1)
c  only here because all the rest is electro
c
c NB SIGN: dpot(j) = -df*rx => force(j) = df*rx = f_ji
c = -f_ij, => f_ij = -df*rx

                        virvdw2 = virvdw2 - (df*rx*rx + df*ry*ry
     &                             + df*rz*rz)
cc

                        df = df - qo*s*(s2*erfcc + s*ewaldcof*derfc)

c non-interpolated erfc
c                       call erfcfun(xerfc,erfc)
c                       e2 = oochrg*s*erfc
c                       df = df - qo*s*(s2*erfc+s*term*exp(-xerfc**2))

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz
                        e_lsym = e_lsym + e2

c electrostatic interaction of Oi1-Hj1
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        
                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        
                        e_lsym = e_lsym + e2
                        
c electrostatic interaction of Oi1-Hj2
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Oj
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Oj
                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Hj1
                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) -dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Hj2
                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Hj1
                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

200             continue
                end if
                jbeg1 = jend1 + 1
                dpot(1,io)  = dpot(1,io)  + dxio
                dpot(2,io)  = dpot(2,io)  + dyio
                dpot(3,io)  = dpot(3,io)  + dzio
                dpot(1,ih1) = dpot(1,ih1) + dxih1
                dpot(2,ih1) = dpot(2,ih1) + dyih1
                dpot(3,ih1) = dpot(3,ih1) + dzih1
                dpot(1,ih2) = dpot(1,ih2) + dxih2
                dpot(2,ih2) = dpot(2,ih2) + dyih2
                dpot(3,ih2) = dpot(3,ih2) + dzih2
300         continue
         istart = iblckwt1(symidx) + 1
400     continue


c start second loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c

         istart = 1
         if (prll_on_off) then
          jbeg2 = poinwat2(my_nwat+1) + 1
         else
          jbeg2 = poinwat2(nwaters)+1
         end if
c loop on symmetry operations
c
        do 700 symidx = 1,symnum
             if (symidx.le.4 .or. (.not.metalyes)) then
         do 600 iwater=istart,iblckwt2(symidx)
                jend2  = psymwt2(iwater)
                i      = symmwat2(iwater)


                io = dpoipt(i)-2
                ih1=io+1
                ih2=io+2

                xio  = coor(1,io) + symop(1,symidx)*a
                yio  = coor(2,io) + symop(2,symidx)*b
                zio  = coor(3,io) + symop(3,symidx)*c

                xih1 = coor(1,ih1) + symop(1,symidx)*a
                yih1 = coor(2,ih1) + symop(2,symidx)*b
                zih1 = coor(3,ih1) + symop(3,symidx)*c

                xih2 = coor(1,ih2) + symop(1,symidx)*a
                yih2 = coor(2,ih2) + symop(2,symidx)*b
                zih2 = coor(3,ih2) + symop(3,symidx)*c

                dxio = 0.d0
                dyio = 0.d0
                dzio = 0.d0

                dxih1 = 0.d0
                dyih1 = 0.d0
                dzih1 = 0.d0

                dxih2 = 0.d0
                dyih2 = 0.d0
                dzih2 = 0.d0




              if (jbeg2.le.jend2) then

              do 500 k=jbeg2,jend2
                      j  = listwt2(k)
                      jo = dpoipt(j)-2
                      jh1= jo + 1
                      jh2= jo + 2

                      rx = xio - coor(1,jo)
                      ry = yio - coor(2,jo)
                      rz = zio - coor(3,jo)
                      r2=rx*rx+ry*ry+rz*rz

c if distance of O-O larger than buffer cutoff
c exclude the whole water-water interaction
c

                      if (r2.gt.cutele2) go to 500

c Oxygen - Oxygen part
                        s2  = 1.d0/r2
                        rij=dsqrt(r2)
                        s = 1.0d0/rij

                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = oochrg*s*erfcc
                        df = - qo*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz
                        e_lsym = e_lsym + e2

c electrostatic interaction of Oi1-Hj1
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        e_lsym = e_lsym + e2
                        
c electrostatic interaction of Oi1-Hj2
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Oj
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Oj

                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
                        df = -qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Hj1

                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi1-Hj2
                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Hj1

                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lsym = e_lsym + e2

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
                        df = -qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

500           continue
              end if

                jbeg2 = jend2 + 1

                dpot(1,io)  = dpot(1,io)  + dxio
                dpot(2,io)  = dpot(2,io)  + dyio
                dpot(3,io)  = dpot(3,io)  + dzio
                dpot(1,ih1) = dpot(1,ih1) + dxih1
                dpot(2,ih1) = dpot(2,ih1) + dyih1
                dpot(3,ih1) = dpot(3,ih1) + dzih1
                dpot(1,ih2) = dpot(1,ih2) + dxih2
                dpot(2,ih2) = dpot(2,ih2) + dyih2
                dpot(3,ih2) = dpot(3,ih2) + dzih2


600     continue

c
c HERE COMES THE METAL IMAGE PART
c

        else

         do 6000 iwater=istart,iblckwt2(symidx)
                jend2  = psymwt2(iwater)
                i      = symmwat2(iwater)


                io = dpoipt(i)-2
                ih1=io+1
                ih2=io+2

                        if (coor(2,io) .lt.0) then
                                yio = -b -coor(2,io)
                                yih1= -b -coor(2,ih1)
                                yih2= -b -coor(2,ih2)
                        else
                                yio = b -coor(2,io)
                                yih1= b -coor(2,ih1)
                                yih2= b -coor(2,ih2)
                        end if

                xio  = coor(1,io) + symop(1,symidx)*a
                zio  = coor(3,io) + symop(3,symidx)*c

                xih1 = coor(1,ih1) + symop(1,symidx)*a
                zih1 = coor(3,ih1) + symop(3,symidx)*c

                xih2 = coor(1,ih2) + symop(1,symidx)*a
                zih2 = coor(3,ih2) + symop(3,symidx)*c


              if (jbeg2.le.jend2) then

              do 5000 k=jbeg2,jend2
                      j  = listwt2(k)
                      jo = dpoipt(j)-2
                      jh1= jo + 1
                      jh2= jo + 2

                      rx = xio - coor(1,jo)
                      ry = yio - coor(2,jo)
                      rz = zio - coor(3,jo)
                      r2=rx*rx+ry*ry+rz*rz

c if distance of O-O larger than buffer cutoff
c exclude the whole water-water interaction
c

                      if (r2.gt.cutele2) go to 5000

c Oxygen - Oxygen part
                        s2  = 1.d0/r2
                        rij   = dsqrt(r2)

                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -oochrg*s*erfcc
                        df =  qo*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz
                        e_lmet = e_lmet + e2

c electrostatic interaction of Oi1-Hj1
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -ohchrg*s*erfcc
                        df = qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        e_lmet = e_lmet + e2
                        
c electrostatic interaction of Oi1-Hj2
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -ohchrg*s*erfcc
                        df = qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi1-Oj
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -ohchrg*s*erfcc
                        df = qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi2-Oj

                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = -ohchrg*s*erfcc
                        df = qoh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi1-Hj1

                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -hhchrg*s*erfcc
                        df = qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi1-Hj2
                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -hhchrg*s*erfcc
                        df = qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi2-Hj1

                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -hhchrg*s*erfcc
                        df = qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lmet = e_lmet + e2

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = -hhchrg*s*erfcc
                        df = qh*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lmet = e_lmet + e2

5000           continue
              end if

                jbeg2 = jend2 + 1

6000    continue
         end if

        istart = iblckwt2(symidx) + 1

700     continue

        return 
        end
