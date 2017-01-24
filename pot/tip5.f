        Subroutine watwat()
        implicit none

C       cdie 

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

        double precision rx,ry,rz,r2,s2,a,b,e1,e2,df
        double precision s,s6
        double precision xio,yio,zio,dxio,dyio,dzio
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,llchrg,lhchrg,hhchrg

        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,jh1,jh2
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,lhchrg,llchrg,hhchrg

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

                   hhchrg = kofdie/eps * ptchg(ih1)*ptchg(ih1) 
                   llchrg = hhchrg
                   lhchrg = -lhchrg
                   choh = cos(104.52/2)
                   shoh = sin(104.52/2)
                   clol = cos(109.47/2)
                   slol = sin(109.47/2)
                   rol  = 0.7
                   roh  = 0.9572
                   inv_axis_length = 1.d0/2.d0*roh*roh*(1+choh)
                   inv_vprod_length  = 1.d0/roh*roh*shoh
                first = .false.
        end if


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

c
c construct lone pair coordinates for one water molecule
c
                dxioh1 = xio - xih1
                dyioh1 = yio - yih1
                dzioh1 = zio - zih1

                dxioh2 = xio - xih2
                dyioh2 = yio - yih2
                dzioh2 = zio - zih2

                qxhat = (dxioh1 + dxioh2))*inv_axis_length
                qyhat = (dyioh1 + dyioh2))*inv_axis_length
                qzhat = (dzioh1 + dzioh2))*inv_axis_length

		txhat = (dyioh1*dzioh2-dzioh1*dyioh2)*inv_vprod_length
                tyhat = (dxioh1*dzioh2-dzioh1*dxioh2)*inv_vprod_length
                tzhat = (dxioh1*dyioh2-dyioh1*dxioh2)*inv_vprod_length

                lpx1 = qxhat*clol + txhat*slol
                lpy1 = qyhat*clol + tyhat*slol
                lpz1 = qzhat*clol + tzhat*slol

                lpx2 = qxhat*clol - txhat*slol
                lpy2 = qyhat*clol - tyhat*slol
                lpz2 = qzhat*clol - tzhat*slol

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

                        s = dsqrt(s2)

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


c the electrostatic part should be done without buffering
100                     continue
c Oxygen - Oxygen part
                        e2  = oochrg*s
                        df = df - e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz
                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                       virXX = virXX -(xmol(pmol(io))-xmol(pmol(jo)))*rx
                       virYY = virYY -(ymol(pmol(io))-ymol(pmol(jo)))*ry
                       virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jo)))*rz
                       virial= virial -
     1                   (xmol(pmol(io))-xmol(pmol(jo)))*rx -
     2                   (ymol(pmol(io))-ymol(pmol(jo)))*ry -
     3                   (zmol(pmol(io))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Oi1-Hj1
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        
                        e2  = ohchrg*s
                        df  =  -e2*s2
                        
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        e_el = e_el + e2
                        
c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Oi1-Hj2
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi1-Oj
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

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

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

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

c electrostatic interaction of Hi2-Oj
                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

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

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

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
                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

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
                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

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
                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

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
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

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

                      rx = xio - coor(1,jo)
                      ry = yio - coor(2,jo)
                      rz = zio - coor(3,jo)
                      r2=rx*rx+ry*ry+rz*rz

c if distance of O-O larger than buffer cutoff
c exclude the whole water-water interaction
c
                      if (r2.gt.cutele2) go to 210

c Oxygen - Oxygen part
                        s2  = 1.d0/r2
                        s   = dsqrt(s2)
                        e2  = oochrg*s
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
                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jo)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jo)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jo)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jo)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jo)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jo)))*rz
                    endif

c electrostatic interaction of Oi1-Hj1
10                      continue
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz
c                       if (r2.gt.cutele2) go to 20
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)
                        
                        e2  = ohchrg*s
                        df  =  -e2*s2
                        
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh1) = dpot(1,jh1) - rx
                        dpot(2,jh1) = dpot(2,jh1) - ry
                        dpot(3,jh1) = dpot(3,jh1) - rz
                        
                        e_el = e_el + e2
                        
c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh1)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh1)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh1)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh1)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh1)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh1)))*rz
                    endif

c electrostatic interaction of Oi1-Hj2
20                      continue
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 30

                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        e2  = ohchrg*s
                        df  =  -e2*s2

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz

                        dxio = dxio + rx
                        dyio = dyio + ry
                        dzio = dzio + rz

                        dpot(1,jh2) = dpot(1,jh2) - rx
                        dpot(2,jh2) = dpot(2,jh2) - ry
                        dpot(3,jh2) = dpot(3,jh2) - rz

                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                     virXX = virXX -(xmol(pmol(io))-xmol(pmol(jh2)))*rx
                     virYY = virYY -(ymol(pmol(io))-ymol(pmol(jh2)))*ry
                     virZZ = virZZ -(zmol(pmol(io))-zmol(pmol(jh2)))*rz
                     virial= virial -
     1                 (xmol(pmol(io))-xmol(pmol(jh2)))*rx -
     2                 (ymol(pmol(io))-ymol(pmol(jh2)))*ry -
     3                 (zmol(pmol(io))-zmol(pmol(jh2)))*rz
                    endif

c electrostatic interaction of Hi1-Oj
30                      continue
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

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

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

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

c electrostatic interaction of Hi2-Oj
40                      continue
                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

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

                        dpot(1,jo) = dpot(1,jo) - rx
                        dpot(2,jo) = dpot(2,jo) - ry
                        dpot(3,jo) = dpot(3,jo) - rz

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
50                      continue

                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

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

                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

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

                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

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

                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

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
