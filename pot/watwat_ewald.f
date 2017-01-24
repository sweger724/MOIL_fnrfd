        Subroutine watwat_ewald()
        implicit none

C       This subrouitne calculates the rpart water contributions to
C       the direct sum within PME algorithm 

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
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/PRESS.BLOCK'

        double precision rx,ry,rz,r2,s2,a,b,e1,e2,df
        double precision s,s6,spi,rij,xerfc,erfcc
        double precision xio,yio,zio,dxio,dyio,dzio,dx,derfc
        double precision xih1,yih1,zih1,dxih1,dyih1,dzih1
        double precision xih2,yih2,zih2,dxih2,dyih2,dzih2
        double precision oosgm6,oosgm12,oochrg,ohchrg,hhchrg
        double precision qo,qh,qoh,term,ovt,del

        double precision rmolx,rmoly,rmolz
        double precision epstmp,pick

        integer i,j,k,jbeg1,jend1,jbeg2,jend2
        integer io,ih1,ih2,jo,jh1,jh2,ind
        logical first
        data first/.true./
        save first
        save oosgm6,oosgm12,ohchrg,oochrg,hhchrg
        save spi,qo,qh,qoh,term,ovt,del
        save epstmp,pick

        if (first) then
                i = realmono(idxtip3(nwaters))
                spi=1.0d0/pi
                spi=dsqrt(spi)
                term=2.0d0*spi*ewaldcof
                ovt=1.0d0/3.0d0
                del=1.0d0/erftbdns
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
                   
                   qo=oochrg
                   qh=hhchrg
                   qoh=ohchrg
                   
                first = .false.
        end if

c
c To speed up water calculations, we do not support for
c TIP3-TIP3 interactions other dielectric constant than 1
c or applying LES. If you need LES water create a water monomer
c with a name different from TIP3.
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



c
c now electrostatic and vdw contributions to the direct sum
c note - current e_el does not include receiprocal cont., whenever
c they are splited
c               
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

c                       s = dsqrt(s2)

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

c virial calculation
        if (pressON) then
c -df*rx is the force applied on particle i by particle j in the x direction
         virXX = virXX -
     1      (xmol(pmol(io))-xmol(pmol(jo)))*df*rx
         virYY = virYY -
     1      (ymol(pmol(io))-ymol(pmol(jo)))*df*ry
         virZZ = virZZ -
     1      (zmol(pmol(io))-zmol(pmol(jo)))*df*rz
         virial = virial -
     1      (xmol(pmol(io))-xmol(pmol(jo)))*df*rx -
     2      (ymol(pmol(io))-ymol(pmol(jo)))*df*ry -
     3      (zmol(pmol(io))-zmol(pmol(jo)))*df*rz
        endif

c the electrostatic part should be done without buffering
c
100                     continue
c
c Oxygen - Oxygen part
c
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
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = oochrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jo))
                  rmoly = ymol(pmol(io))-ymol(pmol(jo))
                  rmolz = zmol(pmol(io))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
                        df = df - qo*s*(s2*erfcc + s*ewaldcof*derfc)

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(io))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(io))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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
                        
                        e_el = e_el + e2
                        
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
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(io))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(io))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = - erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jo))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jo))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                       
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jo))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jo))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

c electrostatic interaction of Hi2-Hj2
                        rx = xih2 - coor(1,jh2)
                        ry = yih2 - coor(2,jh2)
                        rz = zih2 - coor(3,jh2)

                        r2 = rx*rx + ry*ry + rz*rz
                        s2 = 1.0d0/r2
                        s  = dsqrt(s2)

                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2
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
                        rij=dsqrt(r2)
                        s = 1.0d0/rij
                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = oochrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = oochrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jo))
                  rmoly = ymol(pmol(io))-ymol(pmol(jo))
                  rmolz = zmol(pmol(io))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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
                        e_el = e_el + e2

c electrostatic interaction of Oi1-Hj1
10                      continue
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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(io))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(io))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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
                        
                        e_el = e_el + e2
                        
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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(io))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(io))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(io))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jo))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jo))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                       
                        e2 = ohchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = ohchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jo))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jo))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jo))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                        
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih1))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(ih1))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(ih1))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jh1))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jh1))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jh1))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
                        if (ind.lt.0) write(6,*) "ERROR:",rij
                        dx = xerfc - (ind-1)*del
                        derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
     
                        e2 = hhchrg*s*erfcc
ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = hhchrg*(erfcc*s*s*s +
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(ih2))-xmol(pmol(jh2))
                  rmoly = ymol(pmol(ih2))-ymol(pmol(jh2))
                  rmolz = zmol(pmol(ih2))-zmol(pmol(jh2))
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc
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

                        e_el = e_el + e2

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
c
        return 
        end
