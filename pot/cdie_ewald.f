        Subroutine cdie_ewald()
        implicit none

C       this is version of cdie subroutine for PME
c       vdw interactions are not modified 

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

        double precision pick 
        double precision epstmp
        double precision rx,ry,rz,r2,s2,a,b,b1,e1,e2,q,df,df1,df2,ai,bi
        double precision s,s6,tmp,e_les_corr,df_les,cutoff,fuse
        double precision qi,spi,rij,xerfc,erfc,erf,del,term,ovt
        double precision xi,yi,zi,dxi,dyi,dzi,derfc,dx,erfcc

        double precision rmolx,rmoly,rmolz

        integer i,j,k,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
        integer ptbeg,ptend,l,iexc,ind,ipair,maxpair,kk



        if (eelyes) then
                epstmp = kofdie/eps
        else
                epstmp =  0.0d0
        end if
        if (nocut) then
                cutvdw2 = 10000.d0
                cutele2 = 10000.d0
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
        e_les_corr=0.d0
        maxpair=maxlespt
        fuse=50.d0
        cutoff=dsqrt(cutele2)

c       
c -------------------------------------------------------------
c       yael

        if (prll_on_off) then   
           ptbeg = dpoipt(monp(my_pe))+1
           if (my_pe.eq.(num_pes-1)) then
              ptend = npt-1
           else
              ptend = dpoipt(monp(my_pe+1))
           endif   
           if (my_pe.eq.0) then
              jbeg1 = 1
              jbeg2 = 1
              jbeg3 = 1
           else
              jbeg1=point1(ptbeg-1)+1
              jbeg2=point2(ptbeg-1)+1
              jbeg3=point3(ptbeg-1)+1
           endif
        else
           ptbeg = 1
           ptend = npt-1
           jbeg1=1
           jbeg2=1
           jbeg3=1
        end if  
c
        if (ptbeg.eq.1) then
          l=0
        else
          l=exc1(ptbeg-1)
        end if
c
        do 400 i=ptbeg,ptend

                jend1 = point1(i)
                tmp = 1.d0/ptwei(i)
                ai  = epsgm12(i)*pick
                if (arith) then
                bi  = epsgm6(i)
                else
                bi  = epsgm6(i)*pick
                endif
                qi  = ptchg(i)*epstmp
                xi  = coor(1,i)
                yi  = coor(2,i)
                zi  = coor(3,i)
                dxi = 0.d0
                dyi = 0.d0
                dzi = 0.d0
                ipair=0
c
c now for each particle collect correction terms to rs
c over the exclusion list 1-2, 1-3, 1-4
c
                iexc=exc1(i)-l
                if (iexc.gt.0) then

c       write(*,*)'EXC1(I) = ',exc1(i)

                do 180 k=l+1,exc1(i)
                        j=exc2(k)
                    rx = xi - coor(1,j)
                    ry = yi - coor(2,j)
                    rz = zi - coor(3,j)
                    r2=rx*rx+ry*ry+rz*rz
                    s2=1.0d0/r2
                    rij=dsqrt(r2)
                    q = qi*ptchg(j)
                    s  = 1.0d0/rij

                    if ((lesid(i).ne.0) .and.
     *                      (lesid(j).ne.0)) then
                       if (lesid(i).eq.lesid(j)) then
                        if (rij.gt.cutoff) then
                          ipair=ipair+1
                          if (ipair.gt.maxpair) then
                             write(stdo,*)' MAXPAIR too small'
                             stop
                          end if
                          pair(ipair)=j
                        end if  
                       end if     
                       go to 180
c now all the LES correction terms are in separate loop
                    end if   

                    xerfc = ewaldcof*rij
c cubic spline on erfc,derf
c                       write(*,*)' ERFTBDNS XERFC ',erftbdns,xerfc
                    ind = int(erftbdns*xerfc) + 1
                    dx = xerfc - (ind-1)*del
c                       if(ind.gt.100000) write(*,*)' IND = ',ind

                    derfc = - erf_arr(2,ind) - dx * ( erf_arr(3,ind) 
     $                              + 0.5d0 * dx * erf_arr(4,ind) )
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
                    erf = 1.d0 - erfcc
                    e2 = - q*s *erf
c pressure
c there is no correction term if the relative coordinate is not scaled in the
c constant pressure algorithm
c        if (pressON) then
c         V_PIcor = -q*erf*s*s*s+2.d0*ewaldcof/spi*dexp(-xerfc)*s*s
c         V_PIcorXX = V_PIcorXX + V_PIcor*rx*rx
c         V_PIcorYY = V_PIcorYY + V_PIcor*ry*ry
c         V_PIcorZZ = V_PIcorZZ + V_PIcor*rz*rz
c        endif            

c       ileana

c                   write(*,*) " the first e2 is  ",e2
                    df2 = q*s*(s2*erf - s*ewaldcof*derfc)

c                   call erfcfun(xerfc,erfc)
c                   erf = 1.d0 - erfc
c                   e2 = - q*s*erf
c                   df2 = q*s*(s2*erf - s*term*exp(-xerfc**2))

c  #### note - dispersion not supported as yet !!!
c  ##               e_corr_vdw = e_corr_vdw + e1 
                    df1=0.d0

                    df = df1 + df2

                    rx = df*rx
                    ry = df*ry
                    rz = df*rz
                    dxi = dxi + rx
                    dyi = dyi + ry
                    dzi = dzi + rz
                    dpot(1,j) = dpot(1,j) - rx
                    dpot(2,j) = dpot(2,j) - ry
                    dpot(3,j) = dpot(3,j) - rz
                    e_corr = e_corr + e2

 180            continue
                end if

                l=exc1(i)


c
c now for each LES particle collect correction terms to rs
c over the rest of LES particles
c

c               go to 195
                if (lesflag) then

                if (lesid(i).gt.0) then

                do 190 k=1,numles
                    j=lespt(k)
                    if (.not.(i.lt.j)) go to 190
                    rx = xi - coor(1,j)
                    ry = yi - coor(2,j)
                    rz = zi - coor(3,j)
                    r2=rx*rx+ry*ry+rz*rz
c                   s2=1.0d0/r2
                    rij=dsqrt(r2)

c if two LES copies are very close one to each other, they are replaced
c by one particle located in the center of mass in the receip part - hence
c there is no need for correction term for this pair
                    if (rij.lt.leshole) go to 190
c in ewald_init.f it is assumed that for direct and correction terms
c the max distance is not larger than 1.5*cutoff;
c hence if for LES prts rij>1.5*cutoff (as it may happen) one has to use
c non-interpolated erfc - it is also assumed that if rij>fuse
c they are at infinite distance
                    if (rij.gt.fuse) go to 190
                    s2=1.0d0/r2
                    q = qi*ptchg(j)
                    s  = 1.0d0/rij
                    xerfc = ewaldcof*rij

c cubic spline on erfc,derf
c                   ind = int(erftbdns*xerfc) + 1
c                   dx = xerfc - (ind-1)*del
c                   derfc = - erf_arr(2,ind) - dx*( erf_arr(3,ind) 
c     $                              + 0.5d0*dx*erf_arr(4,ind)) 
c              erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
c     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
c                   erf = 1.d0 - erfcc

                    call erfcfun(xerfc,erfc)
                    erf = 1.d0 - erfc

                    df_les=0.d0
                    if (lesid(i).eq.lesid(j))  then
                      if ((rij.gt.cutoff).and.
     &                   (cplbl(i).eq.cplbl(j)))  then
c and if non-bonded just add tmp*receip = -tmp* e2, df
c this is done for other prts in loop over non-bonded pairs but 
c some LES prts may be far enough one from another 
c

                         do 185 kk=1,ipair
                            if (pair(kk).eq.j) go to 187
 185                     continue   
c
c above loop over exclusion list, containing for a given LES prt i all LES prts j
c having distance larger than cutoff and being covalently related
c checks whether they are bonded or not

                        e2 = tmp*q*s*erf
                        df_les = -tmp*q*s*(s2*erf - 
     &                           s*term*exp(-xerfc**2))

                        e_les_corr = e_les_corr + e2
                      end if      
                    end if
 187                continue

                    e2 = - q*s*erf
                    df2 = q*s*(s2*erf - s*term*exp(-xerfc**2))

                    df = df2 + df_les
                    rx = df*rx
                    ry = df*ry
                    rz = df*rz
                    dxi = dxi + rx
                    dyi = dyi + ry
                    dzi = dzi + rz
                    dpot(1,j) = dpot(1,j) - rx
                    dpot(2,j) = dpot(2,j) - ry
                    dpot(3,j) = dpot(3,j) - rz
                    e_les_corr = e_les_corr + e2

 190            continue

                end if
c               with respect to lesid
                end if
c               with respect to lesflag

 195            continue
c
c        end of correction terms 
c
c------------------------------------------------------
c
c compute the direct contributions for non-bonded pairs
c
       !return
c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
c
                if (jbeg1.le.jend1) then

                do 200 k=jbeg1,jend1
                        j=list1(k)
                        rx = xi - coor(1,j)
                        ry = yi - coor(2,j)
                        rz = zi - coor(3,j)
                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2
                        rij=dsqrt(r2)
                        q = qi*ptchg(j)
                        s  = 1.0d0/rij
                        if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
                          df1 = 0.d0
                          xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                          ind = int(erftbdns*xerfc) + 1
                          dx = xerfc - (ind-1)*del
                         derfc = - erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind))
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                          df_les = 0.d0
                          if ((lesid(i).ne.0) .and.
     *                     (lesid(i) .eq.lesid(j))) then
                            q = q*tmp
c                           add les non-bonded recip contribution
                            erf = 1.d0 - erfcc
                            e2 = q*s*erf
                            df_les = -q*s*(s2*erf - s*ewaldcof*derfc)
                            e_les_corr = e_les_corr + e2
                          end if

                          e2 = q*s*erfcc
                          df2 = -q*s*(s2*erfcc + s*ewaldcof*derfc)

ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = q*(erfcc*s*s*s + 
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(i))-xmol(pmol(j))
                  rmoly = ymol(pmol(i))-ymol(pmol(j)) 
                  rmolz = zmol(pmol(i))-zmol(pmol(j))
                  !write (333,*) 'i j df drmol'
                  !write (333,*) i,j,V_PIdir*rx,rmolx
                  !write (333,*) i,j,V_PIdir*ry,rmoly
                  !write (333,*) i,j,V_PIdir*rz,rmolz
c xx
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
c yy
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
c zz              
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc

                          go to 100
                        end if



c       ileana
                        if (.not.arith) then
                        a= ai*epsgm12(j)
                        b= bi*epsgm6(j)
                        else
                        a= ai*epsgm12(j)
                        b= 0.5d0*(bi+epsgm6(j))
                     endif

                        xerfc = ewaldcof*rij
c cubic spline on erfc,derf
                        ind = int(erftbdns*xerfc) + 1
                        dx = xerfc - (ind-1)*del
                        derfc = - erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                        df_les = 0.d0
                        if ((lesid(i).ne.0) .and.
     *                     (lesid(i) .eq.lesid(j))) then
                           a = a*tmp
                        if(.not.arith) b = b*tmp
                           
                           
                           q = q*tmp
c                          add les non-bonded recip contribution
                           erf = 1.d0 - erfcc
                           e2 = q*s*erf
                           df_les = -q*s*(s2*erf - s*ewaldcof*derfc)
                           e_les_corr = e_les_corr + e2
                        end if

                        e2 = q*s*erfcc
                        df2 = -q*s*(s2*erfcc + s*ewaldcof*derfc)

ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = q*(erfcc*s*s*s +  
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(i))-xmol(pmol(j))
                  rmoly = ymol(pmol(i))-ymol(pmol(j)) 
                  rmolz = zmol(pmol(i))-zmol(pmol(j))
                  !write (333,*) 'i j df drmol'
                  !write (333,*) i,j,V_PIdir*rx,rmolx
                  !write (333,*) i,j,V_PIdir*ry,rmoly
                  !write (333,*) i,j,V_PIdir*rz,rmolz
c xx
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
c yy
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
c zz              
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccci

                        if (.not.arith)then
                        s6=s2*s2*s2
                        a = a*s6*s6
                        b = b*s6
                        e1 = a - b
                        df1 = -6.0d0*s2*(a+e1)
                        e_vdw = e_vdw + e1 

c pressure
c virial calculation
                        if (pressON) then
                         rmolx = xmol(pmol(i))-xmol(pmol(j))
                         rmoly = ymol(pmol(i))-ymol(pmol(j))
                         rmolz = zmol(pmol(i))-zmol(pmol(j))
                         virial = virial -
     1                   (xmol(pmol(i))-xmol(pmol(j)))*df1*rx -
     2                   (ymol(pmol(i))-ymol(pmol(j)))*df1*ry -
     3                   (zmol(pmol(i))-zmol(pmol(j)))*df1*rz
                         virXX = virXX - df1*rx*rmolx
                         virYY = virYY - df1*ry*rmoly
                         virZZ = virZZ - df1*rz*rmolz
                  !write (333,*) 'i j dfvdw drmol'
                  !write (333,*) i,j,-df1*rx,rmolx
                  !write (333,*) i,j,-df1*ry,rmoly
                  !write (333,*) i,j,-df1*rz,rmolz
                        endif
cccccccccccccccccc
                           
                        else

                        s6=s2*s2*s2
                        b1 = b*b*b
                        b1 = b1*b1
                        b = a*b1*s6
                        b1 = b1*b1
                        a = a*b1*s6*s6
                
                        e1 = a - b
                        df1 = -6.0d0*s2*(a+e1)
                        e_vdw = e_vdw + e1

                        endif

c the electrostatic part should be done without buffering
100                     continue

                        df = df1 + df2 + df_les

                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        dxi = dxi + rx
                        dyi = dyi + ry
                        dzi = dzi + rz
                        dpot(1,j) = dpot(1,j) - rx
                        dpot(2,j) = dpot(2,j) - ry
                        dpot(3,j) = dpot(3,j) - rz
                        e_el = e_el + e2

200             continue
                end if
                jbeg1 = jend1 + 1
  
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c


              jend2 = point2(i)

c
c  #### no support for dispersion as yet 
c

              if (jbeg2.le.jend2) then
              do 205 k=jbeg2,jend2
                      j=list2(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutvdw2) go to 205
                        a = ai*epsgm12(j)

c       ileana

                        if(.not.arith) then
                        b = bi*epsgm6(j)
                        else
                           b= 0.5d0*(bi+epsgm6(j))
                           endif

                        if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  then
                         a = a*tmp
                      if(.not.arith) b = b*tmp
                         

                      end if
                      s2=1.0d0/r2
                      s6=s2*s2*s2

                      if(.not.arith)then
                
                        a = a*s6*s6
                        b = b*s6
                        
                        else
                         
                      b1 = b*b*b
                      b1 = b1*b1
                       b = a*b1*s6
                      b1 = b1*b1
                       a = a*b1*s6*s6
                       endif
                     
                      e1 = a - b
                      df = -6.0d0*s2*(a+e1)

c pressure
c virial calculation
                        if (pressON) then
                         rmolx = xmol(pmol(i))-xmol(pmol(j))
                         rmoly = ymol(pmol(i))-ymol(pmol(j))
                         rmolz = zmol(pmol(i))-zmol(pmol(j))
                         virial = virial -
     1                   (xmol(pmol(i))-xmol(pmol(j)))*df*rx -
     2                   (ymol(pmol(i))-ymol(pmol(j)))*df*ry -
     3                   (zmol(pmol(i))-zmol(pmol(j)))*df*rz
                         virXX = virXX - df*rx*rmolx
                         virYY = virYY - df*ry*rmoly
                         virZZ = virZZ - df*rz*rmolz
                  !write (333,*) 'i j dfvdw drmol'
                  !write (333,*) i,j,-df*rx,rmolx
                  !write (333,*) i,j,-df*ry,rmoly
                  !write (333,*) i,j,-df*rz,rmolz
                        endif
cccccccccccccccccc

                      rx = df*rx
                      ry = df*ry
                      rz = df*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz

                      e_vdw = e_vdw + e1
205           continue
              end if
              jbeg2 = jend2 + 1

c start THIRD loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c

              jend3 = point3(i)

              if (jbeg3.le.jend3) then

              do 210 k=jbeg3,jend3
                      j=list3(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      rij=dsqrt(r2)
                      if (r2.gt.cutele2) go to 210

                      q=qi*ptchg(j)
                      s2=1.0d0/r2
                      s  = 1.0d0/rij
                      xerfc = ewaldcof*rij

c cubic spline on erfc,derf
                      ind = int(erftbdns*xerfc) + 1
c@
                        if (ind.lt.0) then
                                write(*,*)' i j rij ',i,j,rij
                        end if
                      dx = xerfc - (ind-1)*del
                      derfc = - erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                              + 0.5d0*dx*erf_arr(4,ind)) 
                erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                   + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)

                      df_les = 0.d0
                      if ((lesid(i).ne.0) .and.
     *                     (lesid(i) .eq.lesid(j))) then
                          q = q*tmp
c                         add les non-bonded recip contribution
                          erf = 1.d0 - erfcc
                          e2 = q*s*erf
                          df_les = -q*s*(s2*erf - s*ewaldcof*derfc)
                          e_les_corr = e_les_corr + e2
                      end if
                      
                      e2 = q*s*erfcc
                      df2 = -q*s*(s2*erfcc + s*ewaldcof*derfc)

ccccccccc pressure in Ewald
                  if (pressON) then
                  V_PIdir = q*(erfcc*s*s*s +  
     1            2.d0*ewaldcof/dsqrt(pi)*dexp(-xerfc**(2.d0))*s*s)
                  rmolx = xmol(pmol(i))-xmol(pmol(j))
                  rmoly = ymol(pmol(i))-ymol(pmol(j)) 
                  rmolz = zmol(pmol(i))-zmol(pmol(j))
                  !write (333,*) 'i j dfd drmol'
                  !write (333,*) i,j,V_PIdir*rx,rmolx
                  !write (333,*) i,j,V_PIdir*ry,rmoly
                  !write (333,*) i,j,V_PIdir*rz,rmolz
c xx
                  V_PIdirXX = V_PIdirXX + V_PIdir*rx*rmolx
c yy
                  V_PIdirYY = V_PIdirYY + V_PIdir*ry*rmoly
c zz              
                  V_PIdirZZ = V_PIdirZZ + V_PIdir*rz*rmolz
                  endif
cccccccccc

                      df = df2 + df_les
                      rx = df*rx
                      ry = df*ry
                      rz = df*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz
                      e_el = e_el + e2

210           continue
              end if
              jbeg3 = jend3 + 1



                dpot(1,i) = dpot(1,i) + dxi
                dpot(2,i) = dpot(2,i) + dyi
                dpot(3,i) = dpot(3,i) + dzi
400     continue
c
        e_corr = e_corr + e_les_corr

c remove electrostatic contributions of the type: VP - its related real prts
        if (vp_flag) call vp_cdie_ewald()


        return 
        end
