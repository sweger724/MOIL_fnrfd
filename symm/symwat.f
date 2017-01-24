        Subroutine symwat()
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

        double precision rx,ry,rz,r2,s2,aoo,boo,e1,e2,df
        double precision s,s6
        double precision xio,yio,zio,dxio,dyio,dzio
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
        double precision epstmp,pick

        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,jh1,jh2
        integer symidx,iwater,istart
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,ohchrg,oochrg,hhchrg
        save epstmp,pick

        if (first) then
                i = realmono(idxtip3(nwaters))
                   io = dpoipt(i)-2
                   ih1 = dpoipt(i)

        if (eelyes) then
                epstmp =  1.0d0
        else
                epstmp =  0.0d0
        end if

        if (evdyes) then
                pick = 1.0d0
        else
                pick = 0.0d0
        end if

                   if (.not.arith) then
                    oosgm6 = pick*epsgm6(io)**2
                    oosgm12 = pick*epsgm12(io)**2
                   else
                    oosgm6 = pick*epsgm12(io)**2*epsgm6(io)**6
                    oosgm12 = pick*epsgm12(io)**2*epsgm6(io)**(12)
                   endif

                   oochrg = epstmp*kofdie/eps * ptchg(io)*ptchg(io)
                   ohchrg = epstmp*kofdie/eps * ptchg(io)*ptchg(ih1)
                   hhchrg = epstmp*kofdie/eps * ptchg(ih1)*ptchg(ih1) 
                first = .false.
        end if
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

                        s = dsqrt(s2)

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
                        e2  = oochrg*s
                        df = (df - e2*s2)

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

c virial calculation
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

c@@                     watene = watene + e2

c@@             write(*,*)' symidx #1 i j e2 ',symidx,io,jo,e2
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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,io,jh1,e2
                        
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

c@@                     watene = watene + e2

c@@             write(*,*)' symidx #1 i j e2 ',symidx,io,jh2,e2
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

c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih1,jo,e2
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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih2,jo,e2

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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih1,jh1,e2

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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #1 i j e2 ',symidx,ih1,jh2,e2

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
c@@                     write(*,*)' #2 i j symidx ',i,j,symidx
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
c@@                     write(*,*)' io jo r2  cutele2 ',io,jo,r2,cutele2
                      if (r2.gt.cutele2) go to 500

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

c@@                     watene = watene + e2

c@@             write(*,*)' symidx #3 i j e2 ',symidx ,io,jo,e2
c electrostatic interaction of Oi1-Hj1
c410                    continue
                        rx = xio - coor(1,jh1)
                        ry = yio - coor(2,jh1)
                        rz = zio - coor(3,jh1)
                        
                        r2 = rx*rx + ry*ry + rz*rz
c                       if (r2.gt.cutele2) go to 420

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

c@@                     watene = watene + e2
                        
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,io,jh1,e2
c electrostatic interaction of Oi1-Hj2
c420                    continue
                        rx = xio - coor(1,jh2)
                        ry = yio - coor(2,jh2)
                        rz = zio - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        
c                       if (r2.gt.cutele2) go to 430

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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,io,jh2,e2

c electrostatic interaction of Hi1-Oj
c430                    continue
                        rx = xih1 - coor(1,jo)
                        ry = yih1 - coor(2,jo)
                        rz = zih1 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 440
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

                        rx = xih2 - coor(1,jo)
                        ry = yih2 - coor(2,jo)
                        rz = zih2 - coor(3,jo)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 450

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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih2,jo,e2

c electrostatic interaction of Hi1-Hj1
c450                    continue

                        rx = xih1 - coor(1,jh1)
                        ry = yih1 - coor(2,jh1)
                        rz = zih1 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 460

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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih1,jh1,e2

c electrostatic interaction of Hi1-Hj2
c460                    continue
                        rx = xih1 - coor(1,jh2)
                        ry = yih1 - coor(2,jh2)
                        rz = zih1 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 470
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

c@@                     watene = watene + e2
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih1,jh2,e2

c electrostatic interaction of Hi2-Hj1
c470                    continue

                        rx = xih2 - coor(1,jh1)
                        ry = yih2 - coor(2,jh1)
                        rz = zih2 - coor(3,jh1)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 480
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
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih2,jh1,e2

c electrostatic interaction of Hi2-Hj2
c480                    continue
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz

c                       if (r2.gt.cutele2) go to 500

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
c@@             write(*,*)' symidx #3 i j e2 ',symidx ,ih2,jh2,e2

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
c Here come the metal image charges
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