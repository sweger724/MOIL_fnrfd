        subroutine symcdie()

c calculating van der Waals and electorstatic energies for symmetry related
c particles. 
C       cdie 
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
        include 'COMMON/PRESS.BLOCK'

        double precision pick 
        double precision epstmp
        double precision xi,yi,zi
        double precision rx,ry,rz,r2,s2,aa1,aa,bb1,bb,e1,e2,qi,df
     1  ,df1,df2,ai,bi
        double precision s,s6


        integer istart,i,j,k,jbeg,jend,symidx,ireal


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
                  xi = coor(1,ireal) + symop(1,symidx)*a
                  yi = coor(2,ireal) + symop(2,symidx)*b
                  zi = coor(3,ireal) + symop(3,symidx)*c

c the symmetry related forces for infinite lattice are only half their
c "real" value
c therefore we multiply by 0.5 below. Note that for finite systems this
c division
c should be omitted

c       ileana
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
c                               write(*,*)" 1st loop:aa is  ",aa
c                               write(*,*)" 1st loop:bb is  ",bb

                                e1 = aa - bb
                                e_vsym = e_vsym + e1 
                                df1 = -6.0d0*s2*(aa+e1)

c       ileana

c                               write(*,*)"1st loop:df1 is ", df1

50                              continue
                                s = dsqrt(s2)
                                e2 = qi*ptchg(j)*s
                                df2 = -e2*s2

c       ileana

c                               write(*,*)"1st loop: df2 is  ",df2
                                df = df1 + df2
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

c virial calculation
        if (pressON) then
         virXX = virXX + 
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
c       write (*,*) '1 pmol(',i,')=',pmol(i),
c     1 'xmol(',pmol(i),')=',xmol(pmol(i)),'rx=',rx
c       write (*,*) '1 pmol(',i,')=',pmol(i),
c     1 'ymol(',pmol(i),')=',ymol(pmol(i)),'ry=',ry
c       write (*,*) '1 pmol(',i,')=',pmol(i),
c     1 'zmol(',pmol(i),')=',zmol(pmol(i)),'rz=',rz
c       write (*,*) '1 pmol(',j,')=',pmol(j),
c     1 'xmol(',pmol(j),')=',xmol(pmol(j)),'rx=',rx
c       write (*,*) '1 pmol(',j,')=',pmol(j),
c     1 'ymol(',pmol(j),')=',ymol(pmol(j)),'ry=',ry
c       write (*,*) '1 pmol(',j,')=',pmol(j),
c     1 'zmol(',pmol(j),')=',zmol(pmol(j)),'rz=',rz
        endif

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
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        yi = coor(2,ireal) + symop(2,symidx)*b
                        zi = coor(3,ireal) + symop(3,symidx)*c

c       ileana

                        
                        ai=epsgm12(ireal)*pick
                        bi=epsgm6(ireal)*pick


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

c                               write(*,*)"2nd loop: e1 is ",e1

                                df = -6.0d0*s2*(aa+e1)
                                rx = df*rx
                                ry = df*ry
                                rz = df*rz
                                dpot(1,ireal) = dpot(1,ireal) + rx
                                dpot(2,ireal) = dpot(2,ireal) + ry
                                dpot(3,ireal) = dpot(3,ireal) + rz
                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_vsym = e_vsym + e1 

c virial calculation
        if (pressON) then
         virXX = virXX + 
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
c       write (*,*) '2 pmol(',i,')=',pmol(i),
c     1 'xmol(',pmol(i),')=',xmol(pmol(i)),'rx=',rx
c       write (*,*) '2 pmol(',i,')=',pmol(i),
c     1 'ymol(',pmol(i),')=',ymol(pmol(i)),'ry=',ry
c       write (*,*) '2 pmol(',i,')=',pmol(i),
c     1 'zmol(',pmol(i),')=',zmol(pmol(i)),'rz=',rz
c       write (*,*) '2 pmol(',j,')=',pmol(j),
c     1 'xmol(',pmol(j),')=',xmol(pmol(j)),'rx=',rx
c       write (*,*) '2 pmol(',j,')=',pmol(j),
c     1 'ymol(',pmol(j),')=',ymol(pmol(j)),'ry=',ry
c       write (*,*) '2 pmol(',j,')=',pmol(j),
c     1 'zmol(',pmol(j),')=',zmol(pmol(j)),'rz=',rz
        endif

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
                                s2 = 1.0d0/r2
                                s  = dsqrt(s2)
                                e2 = qi*ptchg(j)*s
                                df = -e2*s2

c       ileana

c                               write(*,*)"3rd loop: df is ", df

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

c virial calculation
        if (pressON) then
         virXX = virXX + 
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx
         virYY = virYY +
     1         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry
         virZZ = virZZ +
     1         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
         virial = virial +
     1         (xmol(pmol(j))-xmol(pmol(ireal))-symop(1,symidx)*a)*rx +
     2         (ymol(pmol(j))-ymol(pmol(ireal))-symop(2,symidx)*b)*ry +
     3         (zmol(pmol(j))-zmol(pmol(ireal))-symop(3,symidx)*c)*rz
        endif

300                     continue
                end if
        jbeg   = jend + 1

301     continue

c Here comes the metal image charges part
c Note that cutoff MUST be smaller than half of the box size.
c
        else
         do 401 i=istart,iblock3(symidx)
                jend = psym3(i)
                ireal = symreal3(i)
                if (jbeg.le.jend) then
                        xi = coor(1,ireal) + symop(1,symidx)*a
                        if (coor(2,ireal).lt.0.d0) then
                                yi = -b - coor(2,ireal)
                        else
                                yi =  b - coor(2,ireal)
                        end if
                        zi = coor(3,ireal) + symop(3,symidx)*c
c
c An image charge has the reverse sign
c
                        qi=-ptchg(ireal)*epstmp

                        do 400 k=jbeg,jend
                                j=list3(k)
                                rx = xi - coor(1,j)
                                ry = yi - coor(2,j)
                                rz = zi - coor(3,j)
                                r2=rx*rx+ry*ry+rz*rz
                                if (r2.gt.cutele2 .and.
     1                                  (.not.nocut)) go to 400
                                s2 = 1.0d0/r2
                                s  = dsqrt(s2)
                                e2 = qi*ptchg(j)*s
                                df = -e2*s2

                                rx = df*rx
                                ry = df*ry
                                rz = df*rz
c
c No forces on the image particle
c in contrast to usual symmetry operation
c
                                dpot(1,j) = dpot(1,j) - rx
                                dpot(2,j) = dpot(2,j) - ry
                                dpot(3,j) = dpot(3,j) - rz
                                e_lmet = e_lmet + e2
400                     continue
                end if
        jbeg   = jend + 1

401     continue

c End of Metal(y/n) if
c
        end if

        istart = iblock3(symidx) + 1
302     continue


        return 
        end


















