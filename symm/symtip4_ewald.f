        subroutine symtip4_ewald()
        implicit none

C       cdie 

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/PRESS.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/EWALD.BLOCK'

        double precision rx,ry,rz,r,r2,s2,aoo,boo,e1,e2,df
        double precision s,s6,spi,xerfc,erfcc
        double precision xio,yio,zio,dxio,dyio,dzio,dx,derfc
        double precision xjo,yjo,zjo,xjh1,yjh1,zjh1
        double precision xjh2,yjh2,zjh2
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
        double precision mscale,drmo,drmh,term,ovt,del

        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,jh1,jh2
        integer symidx,iwater,istart,ind
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,ohchrg,oochrg,hhchrg
        save drmo,drmh,term.ovt.del

        if (first) then
                spi=1.0d0/pi
                spi=dsqrt(spi)
                term=2.0d0*spi*ewaldcof
                ovt=1.0d0/3.0d0
                del=1.0d0/erftbdns

                i = realmono(idxtip3(nwaters))

                   io = dpoipt(i)-2
                   ih1 = dpoipt(i)

                   if (.not.arith) then
                    oosgm6 = epsgm6(io)**2
                    oosgm12 = epsgm12(io)**2
                   else
                    oosgm6 = epsgm12(io)**2*epsgm6(io)**6
                    oosgm12 = epsgm12(io)**2*epsgm6(io)**(12)
                   endif

                   oochrg = kofdie/eps * ptchg(io)*ptchg(io)
                   ohchrg = kofdie/eps * ptchg(io)*ptchg(ih1)
                   hhchrg = kofdie/eps * ptchg(ih1)*ptchg(ih1) 
                   mscale = 0.125d0/(dcos(104.52d0/(2.d0*pi180))
     1                    *sqrt(reqoh2))
                   drmo = 1.d0-mscale
                   drmh = 0.5d0*mscale
		   first = .false.
        end if

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

                        xjo = coor(1,jo)
                        yjo = coor(2,jo)
                        zjo = coor(3,jo)
                        xjh1 = coor(1,jh1)
                        yjh1 = coor(2,jh1)
                        zjh1 = coor(3,jh1)
                        xjh2 = coor(1,jh2)
                        yjh2 = coor(2,jh2)
                        zjh2 = coor(3,jh2)

c take care of oxygen - oxygen van der Waals
                        rx = xio - xjo
                        ry = yio - yjo
                        rz = zio - zjo

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

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

c the electrostatic part should be done without buffering
100                     continue
c M - M part
c compute Rm = (1-c)*Ro + [(Rh1+Rh2)*c/2] and store it at Ro.
c Note that M projects forces on the oxygen as well as on the hydrogen which ocmplicates the code
c and make it somehwat more expensive.
c After the Van der Waals interaction was computed no need
c for Ro anymore. This requires however that the force and energies
c for the van der Waals interaction will be accumulated for the vdW now (than later).
c
                        xio = drmo*xio + drmh*(xih1+xih2)
                        yio = drmo*yio + drmh*(yih1+yih2)
                        zio = drmo*zio + drmh*(zih1+zih2)
                        xjo = drmo*xjo + drmh*(xjh1+xjh2)
                        yjo = drmo*yjo + drmh*(yjh1+yjh2)
                        zjo = drmo*zjo + drmh*(zjh1+zjh2)
                        rx = xio - xjo
                        ry = yio - yjo
                        rz = zio - zjo
                        r2 = rx*rx + ry*ry + rz*rz
                        r = dsqrt(r2)
                        s = 1.0d0/r
                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = oochrg*s*erfcc

                        df = df - oochrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz
                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz
                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz

                        dpot(1,jo) = dpot(1,jo) - drmo*rx
                        dpot(2,jo) = dpot(2,jo) - drmo*ry
                        dpot(3,jo) = dpot(3,jo) - drmo*rz
                        dpot(1,jh1) = dpot(1,jh1) - drmh*rx
                        dpot(2,jh1) = dpot(2,jh1) - drmh*ry
                        dpot(3,jh1) = dpot(3,jh1) - drmh*rz
                        dpot(1,jh2) = dpot(1,jh2) - drmh*rx
                        dpot(2,jh2) = dpot(2,jh2) - drmh*ry
                        dpot(3,jh2) = dpot(3,jh2) - drmh*rz
                        e_lsym = e_lsym + e2


c virial calculation (not ready for tip4p)
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(io))-symop(1,symidx)*a)*rx+
     2         (ymol(pmol(jo))-ymol(pmol(io))-symop(2,symidx)*b)*ry+
     3         (zmol(pmol(jo))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif

c electrostatic interaction of Mi1-Hj1
                        rx = xio - xjh1
                        ry = yio - yjh1
                        rz = zio - zjh1
                        
                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = 1/r
                        
                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohcgrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz
                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz
                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(io))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(io))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif

                        
c electrostatic interaction of Mi1-Hj2
                        rx = xio - xjh2
                        ry = yio - yjh2
                        rz = zio - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = 1.d0/r

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc


                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz

                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz

                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz


                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(io))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(io))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif

c electrostatic interaction of Hi1-Mj
                        rx = xih1 - xjo
                        ry = yih1 - yjo
                        rz = zih1 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        r  = dsqrt(r2)

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)


                        e2  = ohchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jo) = dpot(1,jo) - dmro*rx
                        dpot(2,jo) = dpot(2,jo) - dmro*ry
                        dpot(3,jo) = dpot(3,jo) - dmro*rz

                        dpot(1,jh1) = dpot(1,jh1) - dmrh*rx
                        dpot(2,jh1) = dpot(2,jh1) - dmrh*ry
                        dpot(3,jh1) = dpot(3,jh1) - dmrh*rz

                        dpot(1,jh2) = dpot(1,jh2) - dmrh*rx
                        dpot(2,jh2) = dpot(2,jh2) - dmrh*ry
                        dpot(3,jh2) = dpot(3,jh2) - dmrh*rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jo))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jo))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif

c electrostatic interaction of Hi2-Mj
                        rx = xih2 - xjo
                        ry = yih2 - yjo
                        rz = zih2 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = 1.d0/r

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jo) = dpot(1,jo) - drmo*rx
                        dpot(2,jo) = dpot(2,jo) - drmo*ry
                        dpot(3,jo) = dpot(3,jo) - drmo*rz
                        dpot(1,jh1) = dpot(1,jh1) - drmh*rx
                        dpot(2,jh1) = dpot(2,jh1) - drmh*ry
                        dpot(3,jh1) = dpot(3,jh1) - drmh*rz
                        dpot(1,jh2) = dpot(1,jh2) - drmh*rx
                        dpot(2,jh2) = dpot(2,jh2) - drmh*ry
                        dpot(3,jh2) = dpot(3,jh2) - drmh*rz

                        e_lsym = e_lsym + e2

c virial calculation (not ready for tip4)
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jo))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jo))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi1-Hj1
                        rx = xih1 - xjh1
                        ry = yih1 - yjh1
                        rz = zih1 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi1-Hj2
                        rx = xih1 - xjh2
                        ry = yih1 - yjh2
                        rz = zih1 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi2-Hj1
                        rx = xih2 - xjh1
                        ry = yih2 - yjh1
                        rz = zih2 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = 1.d0/r

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih2,jh1,e2

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - xjh2
                        ry = yih2 - yjh2
                        rz = zih2 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s2 = 1.0d0/r2
                        s  = 1.0d0/r

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih2,jh2,e2
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

                        xjo = coor(1,jo)
                        yjo = coor(2,jo)
                        zjo = coor(3,jo)
                        xjh1 = coor(1,jh1)
                        yjh1 = coor(2,jh1)
                        zjh1 = coor(3,jh1)
                        xjh2 = coor(1,jh2)
                        yjh2 = coor(2,jh2)
                        zjh2 = coor(3,jh2)
                        xio = drmo*xio + drmh*(xih1+xih2)
                        yio = drmo*yio + drmh*(yih1+yih2)
                        zio = drmo*zio + drmh*(zih1+zih2)
                        xjo = drmo*xjo + drmh*(xjh1+xjh2)
                        yjo = drmo*yjo + drmh*(yjh1+yjh2)
                        zjo = drmo*zjo + drmh*(zjh1+zjh2)
                        rx = xio - xjo
                        ry = yio - yjo
                        rz = zio - zjo

                      r2=rx*rx+ry*ry+rz*rz

c if distance of M-M larger than buffer cutoff
c exclude the whole water-water interaction
c
                      if (r2.gt.cutele2) go to 500

c Mi - Mj part
                        r   = dsqrt(r2)
                        s   = 1.d0/r
                        s2  = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = oochrg*s*erfcc
                        df = -oochrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz
                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz
                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz

                        dpot(1,jo) = dpot(1,jo) - drmo*rx
                        dpot(2,jo) = dpot(2,jo) - drmo*ry
                        dpot(3,jo) = dpot(3,jo) - drmo*rz
                        dpot(1,jh1) = dpot(1,jh1) - drmh*rx
                        dpot(2,jh1) = dpot(2,jh1) - drmh*ry
                        dpot(3,jh1) = dpot(3,jh1) - drmh*rz
                        dpot(1,jh2) = dpot(1,jh2) - drmh*rx
                        dpot(2,jh2) = dpot(2,jh2) - drmh*ry
                        dpot(3,jh2) = dpot(3,jh2) - drmh*rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(io))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jo))-ymol(pmol(io))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jo))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif

c electrostatic interaction of Mi1-Hj1
                        rx = xio - xjh1
                        ry = yio - yjh1
                        rz = zio - zjh1
                        
                        r2 = rx*rx + ry*ry + rz*rz
                        r = dsqrt(r2)
                        s = 1.d0/r
                        s2 = s*s
                        
                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz
                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz
                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(io))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(io))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif

c electrostatic interaction of Mi1-Hj2
c420                    continue
                        rx = xio - xjh2
                        ry = yio - yjh2
                        rz = zio - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r = dsqrt(r2)
                        s = 1.d0/r
                        s2 = s*s
                        
                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + drmo*rx
                        dyio = dyio + drmo*ry
                        dzio = dzio + drmo*rz

                        dxih1 = dxih1 + drmh*rx
                        dyih1 = dyih1 + drmh*ry
                        dzih1 = dzih1 + drmh*rz

                        dxih2 = dxih2 + drmh*rx
                        dyih2 = dyih2 + drmh*ry
                        dzih2 = dzih2 + drmh*rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(io))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(io))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(io))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(io))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(io))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(io))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi1-Mj
c430                    continue
                        rx = xih1 - xjo
                        ry = yih1 - yjo
                        rz = zih1 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        r = dsqrt(r2)
                        s = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jo) = dpot(1,jo) - dmro*rx
                        dpot(2,jo) = dpot(2,jo) - dmro*ry
                        dpot(3,jo) = dpot(3,jo) - dmro*rz

                        dpot(1,jh1) = dpot(1,jh1) - dmrh*rx
                        dpot(2,jh1) = dpot(2,jh1) - dmrh*ry
                        dpot(3,jh1) = dpot(3,jh1) - dmrh*rz

                        dpot(1,jh2) = dpot(1,jh2) - dmrh*rx
                        dpot(2,jh2) = dpot(2,jh2) - dmrh*ry
                        dpot(3,jh2) = dpot(3,jh2) - dmrh*rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jo))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jo))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih1,jo,e2

c electrostatic interaction of Hi2-Oj

c440                    continue

                        rx = xih2 - xjo
                        ry = yih2 - yjo
                        rz = zih2 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s  = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = ohchrg*s*erfcc
                        df = -ohchrg*s*(s2*erfcc + s*ewaldcof*derfc)

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jo) = dpot(1,jo) - drmo*rx
                        dpot(2,jo) = dpot(2,jo) - drmo*ry
                        dpot(3,jo) = dpot(3,jo) - drmo*rz
                        dpot(1,jh1) = dpot(1,jh1) - drmh*rx
                        dpot(2,jh1) = dpot(2,jh1) - drmh*ry
                        dpot(3,jh1) = dpot(3,jh1) - drmh*rz
                        dpot(1,jh2) = dpot(1,jh2) - drmh*rx
                        dpot(2,jh2) = dpot(2,jh2) - drmh*ry
                        dpot(3,jh2) = dpot(3,jh2) - drmh*rz

                        e_lsym = e_lsym + e2

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jo))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jo))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jo))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jo))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jo))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jo))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi1-Hj1
c450                    continue

                        rx = xih1 - xjh1
                        ry = yih1 - yjh1
                        rz = zih1 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s  = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)


                        e2  = hhchrg*s
                        df  =  -e2*s2

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi1-Hj2
c460                    continue
                        rx = xih1 - xjh2
                        ry = yih1 - yjh2
                        rz = zih1 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s  = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(ih1))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(ih1))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(ih1))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi2-Hj1
c470                    continue

                        rx = xih2 - xhj1
                        ry = yih2 - yjh1
                        rz = zih2 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s  = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)

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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh1))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh1))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh1))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh1))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh1))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh1))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif


c electrostatic interaction of Hi2-Hj2
c480                    continue
                        rx = xih2 - xjh2
                        ry = yih2 - yjh2
                        rz = zih2 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        r  = dsqrt(r2)
                        s  = 1.d0/r
                        s2 = s*s

                        xerfc = ewaldcof*r
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind)
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        e2 = hhchrg*s*erfcc
                        df = -hhchrg*s*(s2*erfcc + s*ewaldcof*derfc)


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

c virial calculation
        if (pressON) then
         virXX = virXX +
     1         (xmol(pmol(jh2))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(jh2))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(jh2))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(jh2))-xmol(pmol(ih2))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(jh2))-ymol(pmol(ih2))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(jh2))-zmol(pmol(ih2))-symop(3,symidx)*c)*rz
        endif


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
c Here come the metal image charges (THIS IS NOT WORKING WITH TIP4. Ron)
c
        else

          do 6000 iwater=istart,iblckwt2(symidx)
                jend2  = psymwt2(iwater)
                i      = symmwat2(iwater)


                io = dpoipt(i)-2
                ih1=io+1
                ih2=io+2

c
c the positions of the image charges of the whole water
c molecules are determined
c by the oxygen position only. It is hoped
c that all water atoms are sufficiently closed.
c
                        if (coor(2,io) .lt. 0.d0 ) then
                         yio  = -b - coor(2,io)
                         yih1 = -b - coor(2,ih1)
                         yih2 = -b - coor(2,ih2)
                        else
                         yio  = b - coor(2,io)
                         yih1 = b - coor(2,ih1)
                         yih2 = b - coor(2,ih2)
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
                        s   = dsqrt(s2)
c
c note the inverse charge on the image particle
c
                        e2  = -oochrg*s
                        df  = -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz
                e_lmet = e_lmet + e2

c electorstatic interaction of Oi1 & Hj1

                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        
c note change in sign for image charge

                        e2  = -ohchrg*s
                        df  =  -e2*s2
                        
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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -ohchrg*s
                        df  =  -e2*s2

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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -ohchrg*s
                        df  =  -e2*s2

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
                        s  = dsqrt(s2)


c note change in sign for image charge

                        e2  = -ohchrg*s
                        df  =  -e2*s2

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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -hhchrg*s
                        df  =  -e2*s2

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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -hhchrg*s
                        df  = -e2*s2

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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -hhchrg*s
                        df  =  -e2*s2

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
                        s  = dsqrt(s2)

c note change in sign for image charge

                        e2  = -hhchrg*s
                        df  =  -e2*s2

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
