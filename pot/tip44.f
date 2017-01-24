        subroutine tip44()
        implicit none

c tip4 - tip4 interaction refined for ewald Horn et al JCP 120,9665(2004)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/PRESS.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'

        double precision rx,ry,rz,r2,s2,a,b,e1,e2,df
        double precision s,s6
        double precision xio,yio,zio,dxio,dyio,dzio
        double precision xjo,yjo,zjo,xjh1,yjh1,zjh1
        double precision xjh2,yjh2,zjh2
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
        double precision mscale,drmo,drmh

        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,h1,jh2
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,ohchrg,oochrg,hhchrg
        save drmo,drmh

        if (first) then
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
c
c water - water interactions for TIP4P model.
c


        do 400 i=1,my_nwat
                jbeg1  = poinwat1(i)+1
                jbeg2  = poinwat2(i)+1
                jend1  = poinwat1(i+1)
                jend2  = poinwat2(i+1)


                io = dpoipt(indxwat(i))-2
                ih1=io+1
                ih2=io+2

                xio  = coor(1,io)
                yio  = coor(2,io)
                zio  = coor(3,io)

                xih1 = coor(1,ih1)
                yih1 = coor(2,ih1)
                zih1 = coor(3,ih1)

                xih2 = coor(1,ih2)
                yih2 = coor(2,ih2)
                zih2 = coor(3,ih2)

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
                        xjo = coor(1,jo)
                        yjo = coor(2,jo)
                        zjo = coor(3,jo)
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
                        a = oosgm12*s6*s6
                        b = oosgm6*s6
                        e1 = a - b
                        df  = -6.0d0*s2*(a+e1)
                        e_vdw = e_vdw + e1 
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
c compute Rm = Ro + [(Rh1+Rh2)/2 -Ro]*c and store it at Ro.
c Note that M projects forces on the oxygen as well as on the hydrogen which ocmplicates the code
c and make it somehwat more expensive.
c After the Van der Waals interaction was computed no need
c for Ro anymore. This requires however that the force and energies
c for the van der Waals interaction will be accumulated for the vdW now (than later).
c
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
                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.d0/r2
                        s  = dsqrt(s2)
                        e2  = oochrg*s
                        df = df - e2*s2

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
                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
c NOT ADAPTED FOR TIP4P!! Ron.
                       virXX = virXX -(xmol(pmol(io))-xmol(pmol(jo)))*rx
                       virYY = virYY -(ymol(pmol(io))-ymol(pmol(jo)))*ry
                       virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jo)))*rz
                       virial= virial -
     1                   (xmol(pmol(io))-xmol(pmol(jo)))*rx -
     2                   (ymol(pmol(io))-ymol(pmol(jo)))*ry -
     3                   (zmol(pmol(io))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Mi-Hj1
                        rx = xio - xjh1
                        ry = yio - yjh1
                        rz = zio - zjh1
                        
                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        
                        e2  = ohchrg*s
                        df  =  -e2*s2
                        

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
                        
                        e_el = e_el + e2
                        
c virial calculation
c Not verified for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Mi-Hj2
                        rx = drmo*xio+drmh*(xih1+xih2) - xjh2
                        ry = drmo*yio+drmh*(yih1+yih2) - yjh2
                        rz = drmo*zio+drmh*(zih1+zih2) - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

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

                        e_el = e_el + e2

c virial calculation
c Not adjusted for TIP4P
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi1-Mj
                        rx = xih1 - xjo
                        ry = yih1 - yjo
                        rz = zih1 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jo) = dpot(1,jo) - drmo*rx
                        dpot(2,jo) = dpot(2,jo) - drmo*ry
                        dpot(3,jo) = dpot(3,jo) - drmo*rz
                        dpot(1,jh1) = dpot(1,jh1) - drmh*rx
                        dpot(2,jh1) = dpot(2,jh1) - drmh*ry
                        dpot(3,jh1) = dpot(3,jh1) - drmh*rz
                        dpot(1,jh2) = dpot(1,jh2) - drmh*rx
                        dpot(2,jh2) = dpot(2,jh2) - drmh*ry
                        dpot(3,jh2) = dpot(3,jh2) - drmh*rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Hi2-Mj
                        rx = xih2 - xjo
                        ry = yih2 - yjo
                        rz = zih2 - zjo

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

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

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Hi1-Hj1
                        rx = xih1 - xjh1
                        ry = yih1 - yjh1
                        rz = zih1 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

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

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Hi1-Hj2
                        rx = xih1 - xjh2
                        ry = yih1 - yjh2
                        rz = zih1 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi2-Hj1
                        rx = xih2 - xjh1
                        ry = yih2 - yjh1
                        rz = zih2 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - xjh2
                        ry = yih2 - yjh2
                        rz = zih2 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jh2)))*rz
                    endif

200             continue
                end if
  

c start second loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c



              if (jbeg2.le.jend2) then

              do 210 k=jbeg2,jend2
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
                      xjo = drmo*xjo+drmh*(xjh1+xjh2)
                      yjo = drmo*yjo+drmh*(yjh1+yjh2)
                      zjo = drmo*zjo+drmh*(zjh1+zjh2)

                      rx = xio - xjo
                      ry = yio - yjo
                      rz = zio - zjo
                      r2=rx*rx+ry*ry+rz*rz

c if distance of Mi-Mj larger than buffer cutoff
c exclude the whole water-water interaction
c
                      if (r2.gt.cutele2) go to 210

c Mi - Mj part
                        s2  = 1.d0/r2
                        s   = dsqrt(s2)
                        e2  = oochrg*s
                        df  = -e2*s2

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
                        e_el = e_el + e2


c virial calculation
c Not adjusted for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Mi1-Hj1
10                      continue
                        rx = xio - xjh1
                        ry = yio - yjh1
                        rz = zio - zjh1
                        
                        r2 = rx*rx + ry*ry + rz*rz
c                       if (r2.gt.cutele2) go to 20
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        
                        e2  = ohchrg*s
                        df  =  -e2*s2
                        
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
                        
                        e_el = e_el + e2
                        
c virial calculation
c Not adjusted for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Mi1-Hj2
20                      continue
                        rx = xio - xjh2
                        ry = yio - yjh2
                        rz = zio - zjh2

                        r2 = rx*rx + ry*ry + rz*rz


                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

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

                        e_el = e_el + e2

c virial calculation
c Not validated for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi1-Mj
30                      continue
                        rx = xih1 - xjo
                        ry = yih1 - yjo
                        rz = zih1 - zjo

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 40

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

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

                        e_el = e_el + e2

c virial calculation
c Sorry, not ready for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Hi2-Mj
40                      continue
                        rx = xih2 - xjo
                        ry = yih2 - yjo
                        rz = zih2 - zjo

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 50

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jo) = dpot(1,jo) - dmro*rx
                        dpot(2,jo) = dpot(2,jo) - dmro*ry
                        dpot(3,jo) = dpot(3,jo) - dmro*rz

                        dpot(1,jh1) = dpot(1,jh1) - dmrh*rx
                        dpot(2,jh1) = dpot(2,jh1) - dmrh*ry
                        dpot(3,jh1) = dpot(3,jh1) - dmrh*rz

                        dpot(1,jh2) = dpot(1,jh2) - dmrh*rx
                        dpot(2,jh2) = dpot(2,jh2) - dmrh*ry
                        dpot(3,jh2) = dpot(3,jh2) - dmrh*rz

                        e_el = e_el + e2

c virial calculation
c Not ready for TIP4P
c
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Hi1-Hj1
50                      continue

                        rx = xih1 - xjh1
                        ry = yih1 - yjh1
                        rz = zih1 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) got o 60

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

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

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Hi1-Hj2
60                      continue

                        rx = xih1 - xjh2
                        ry = yih1 - yjh2
                        rz = zih1 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 70

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih1 = dxih1 + rx
                        dyih1 = dyih1 + ry
                        dzih1 = dzih1 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih1))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(ih1))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(ih1))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih1))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(ih1))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(ih1))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi2-Hj1
70                      continue

                        rx = xih2 - xjh1
                        ry = yih2 - yjh1
                        rz = zih2 - zjh1

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 80

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Hi2-Hj2
80                      continue

                        rx = xih2 - xjh2
                        ry = yih2 - yjh2
                        rz = zih2 - zjh2

                        r2 = rx*rx + ry*ry + rz*rz


c                       if (r2.gt.cutele2) go to 210

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = hhchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxih2 = dxih2 + rx
                        dyih2 = dyih2 + ry
                        dzih2 = dzih2 + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(ih2))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(ih2))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(ih2))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(ih2))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(ih2))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(ih2))-zmol(pmol(jh2)))*rz
                    endif

210           continue
              end if


                dpot(1,io)  = dpot(1,io)  + dxio
                dpot(2,io)  = dpot(2,io)  + dyio
                dpot(3,io)  = dpot(3,io)  + dzio
                dpot(1,ih1) = dpot(1,ih1) + dxih1
                dpot(2,ih1) = dpot(2,ih1) + dyih1
                dpot(3,ih1) = dpot(3,ih1) + dzih1
                dpot(1,ih2) = dpot(1,ih2) + dxih2
                dpot(2,ih2) = dpot(2,ih2) + dyih2
                dpot(3,ih2) = dpot(3,ih2) + dzih2


400     continue
Cyef
C       write(*,*)'watwat contrib',e_el - el_yef

        return 
        end
