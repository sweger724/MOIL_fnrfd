        Subroutine cdie()

C       cdie 
        implicit none

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
c mauro - block for mute
        include 'COMMON/MUTA.BLOCK'
        include 'COMMON/PRESS.BLOCK'

        double precision pick 
        double precision epstmp
        double precision rx,ry,rz,r2,s2,a,b,e1,e2,q,df,df1,df2,ai,bi
        double precision s,s6,tmp
        double precision qi
        double precision xi,yi,zi,dxi,dyi,dzi

        double precision rmolx,rmoly,rmolz

        integer i,j,k,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
        integer ptbeg, ptend


        if (eelyes) then
                epstmp = kofdie/eps
        else
                epstmp =  0.0d0
        end if
        if (nocut) then
                cutvdw2 = 100000.d0
                cutele2 = 100000.d0
        end if

        if (evdyes) then
                pick = 1.0d0
        else
                pick = 0.0d0
        end if

        e_vdw = 0.0d0
        e_el=0.d0

        if (prll_on_off) then
               ptbeg = dpoipt(monp(my_pe))+1
               if (my_pe.eq.(num_pes-1)) then
                       ptend = npt-1
               else
                       ptend = dpoipt(monp(my_pe+1))
               end if
               if (my_pe.eq.0) then
                       jbeg1 = 1
                       jbeg2 = 1
                       jbeg3 = 1
               else
                       jbeg1 = point1(ptbeg-1) + 1
                       jbeg2 = point2(ptbeg-1) + 1
                       jbeg3 = point3(ptbeg-1) + 1
               end if
       else
               ptbeg = 1
               ptend = npt -1
               jbeg1=1
               jbeg2=1
               jbeg3=1
       end if


c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
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

                if (jbeg1.le.jend1) then

                do 200 k=jbeg1,jend1
                        j=list1(k)
                        rx = xi - coor(1,j)
                        ry = yi - coor(2,j)
                        rz = zi - coor(3,j)
                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2

                        q = qi*ptchg(j)
                        s = dsqrt(s2)

                        if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
                                df1 = 0.d0
                                if ((lesid(i).ne.0) .and.
     *                         (lesid(i) .eq.lesid(j))) then
                                 q = q*tmp
                                end if
c mauro - mute - count once lambda in interactions within mutants
                                if (mutaid(i).eq.1 .and. 
     &                              mutaid(j).eq.1)then
                                  q = q/lambda
                                else if (mutaid(i).eq.2 .and. 
     &                                   mutaid(j).eq.2) then
                                  q = q/(1.0d0-lambda)
                                endif

                                e2  = q*s
                                df2 = -e2*s2
                                go to 100
                        end if
c a is 4epsilon(ij)=2sqrt(epsi)*2sqrt(epsj)
c b is 0.5d0*(sigma_i+sigma_j) i.e. sigma(ij)
c
                        a= ai*epsgm12(j)
                          if(arith) then
                        b= 0.5d0*(bi+epsgm6(j))
                          else
                        b= bi*epsgm6(j)
                          endif
                        if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  then
                           a = a*tmp
                        if(.not.arith) b = b*tmp  
                           q = q*tmp
                        end if
c mauro - mute
                        if (mutaid(i).eq.1 .and. mutaid(j).eq.1)then
                            a = a/lambda
                            q = q/lambda
                             if(.not.arith) then
                              b = b/lambda
                             endif
                        else if (mutaid(i).eq.2 .and. 
     &                           mutaid(j).eq.2) then
                             a = a/(1.0d0-lambda)
                             q = q/(1.0d0-lambda)
                             if(.not.arith) then
                              b = b/(1.0d0-lambda)
                             endif
                         endif



                        e2  = q*s
                        df2 = -e2*s2
                        if (arith) then
                        b   = s2*b*b
                        s6  = b*b*b
                        e1  = a*(s6*s6 - s6)
                        df1 = -6.0d0*s2*a*(2.d0*s6*s6-s6)
                        else
                        s6=s2*s2*s2
                        a = a*s6*s6
                        b = b*s6
                        e1 = a - b
                        df1 = -6.0d0*s2*(a+e1)
                        endif                        
               !write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~'
               !write(stdo,*)'i and j are ',i,j,dsqrt(r2)
               !write(stdo,'(a,f18.10)')'e1--vdW is ',e1
               !write(stdo,'(a,f8.4)')'e2--elec is ',e2
               !write(stdo,*)'epsgm12 of i and j are',epsgm12(i),epsgm12(j)
               !write(stdo,*)'epsgm6 of i and j are',epsgm6(i),epsgm6(j)
               !write(stdo,*)'a and b are ',a,b
               !write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~~'

c               write(*,*)"df1 is ",df1
                        e_vdw = e_vdw + e1
c mauro - muta <U2-U1>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1) then
                         tmp_e_lambda=tmp_e_lambda+e1/lambda
                   else if (mutaid(i).eq.2 .or. mutaid(j).eq.2) then
                         tmp_e_lambda=tmp_e_lambda-e1/(1.d0-lambda)
                        endif                       

c the electrostatic part should be done without buffering
100                     continue
                        df = df1 + df2
c                       write(*,*) "df is ", df
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        dxi = dxi + rx
                        dyi = dyi + ry
                        dzi = dzi + rz
                        dpot(1,j) = dpot(1,j) - rx
                        dpot(2,j) = dpot(2,j) - ry
                        dpot(3,j) = dpot(3,j) - rz
c                       write (*,*) "index j is ",j
c                     write (*,*) "dpot1 is",dpot(1,j)
c                     write (*,*) "dpot2 is",dpot(2,j)
c                     write (*,*) "dpot3 is",dpot(3,j)
                        e_el = e_el + e2

c virial calculation
                    if (pressON) then
                       virXX = virXX -(xmol(pmol(i))-xmol(pmol(j)))*rx
                       virYY = virYY -(ymol(pmol(i))-ymol(pmol(j)))*ry
                       virZZ = virZZ -(zmol(pmol(i))-zmol(pmol(j)))*rz
                       virial= virial -
     1                   (xmol(pmol(i))-xmol(pmol(j)))*rx -
     2                   (ymol(pmol(i))-ymol(pmol(j)))*ry -
     3                   (zmol(pmol(i))-zmol(pmol(j)))*rz
                       rmolx = xmol(pmol(i))-xmol(pmol(j))
                       rmoly = ymol(pmol(i))-ymol(pmol(j))
                       rmolz = zmol(pmol(i))-zmol(pmol(j))
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
                  write (333,*) 'i j df drmol'
                  write (333,*) i,j,-rx,rmolx
                  write (333,*) i,j,-ry,rmoly
                  write (333,*) i,j,-rz,rmolz
                        endif

c mauro - muta <U2-U1>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1) then
                         tmp_e_lambda=tmp_e_lambda+e2/lambda
                   else if (mutaid(i).eq.2 .or. mutaid(j).eq.2) then
                         tmp_e_lambda=tmp_e_lambda-e2/(1.d0-lambda)
                        endif


200             continue
                end if
                jbeg1 = jend1 + 1
  
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c


              jend2 = point2(i)


              if (jbeg2.le.jend2) then
              do 205 k=jbeg2,jend2
                      j=list2(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutvdw2) go to 205
                        a = ai*epsgm12(j)
                        if(arith) then
                        b = 0.5d0*(bi + epsgm6(j))
                        else
                        b = bi*epsgm6(j)
                        endif
                        if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  then
                         a = a*tmp
                         if(.not.arith) b = b*tmp
                      end if
c mauro - mute
                        if (mutaid(i).eq.1 .and. mutaid(j).eq.1)then
                            a = a/lambda
                             if(.not.arith) then
                              b = b/lambda
                             endif 
                        else if (mutaid(i).eq.2 .and. 
     &                           mutaid(j).eq.2) then
                             a = a/(1.0d0-lambda)
                             if(.not.arith) then
                              b = b/(1.0d0-lambda)
                             endif
                         endif


                      s2=1.0d0/r2
c a is 4epsilon(ij)
c b is 0.5(sigmai+sigmaj) i.e. sigma(ij)
c
                      if (arith) then
                      b   = s2*b*b
                      s6  = b*b*b
                      e1  = a*(s6*s6 - s6)
                      df = -6.0d0*s2*a*(2.d0*s6*s6-s6)
                      else
                      s6=s2*s2*s2
                      a = a*s6*s6
                      b = b*s6
                      e1 = a - b
                      df = -6.0d0*s2*(a+e1)
                      endif

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
c mauro - muta <U2-U1>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1) then
                         tmp_e_lambda=tmp_e_lambda+e1/lambda
                   else if (mutaid(i).eq.2 .or. mutaid(j).eq.2) then
                         tmp_e_lambda=tmp_e_lambda-e1/(1.d0-lambda)
                        endif

c virial calculation
                        if (pressON) then
                       virXX = virXX -(xmol(pmol(i))-xmol(pmol(j)))*rx
                       virYY = virYY -(ymol(pmol(i))-ymol(pmol(j)))*ry
                       virZZ = virZZ -(zmol(pmol(i))-zmol(pmol(j)))*rz
                         virial = virial -
     1                   (xmol(pmol(i))-xmol(pmol(j)))*rx -
     2                   (ymol(pmol(i))-ymol(pmol(j)))*ry -
     3                   (zmol(pmol(i))-zmol(pmol(j)))*rz
                       rmolx = xmol(pmol(i))-xmol(pmol(j))
                       rmoly = ymol(pmol(i))-ymol(pmol(j))
                       rmolz = zmol(pmol(i))-zmol(pmol(j))
                  write (333,*) 'i j df drmol'
                  write (333,*) i,j,-rx,rmolx
                  write (333,*) i,j,-ry,rmoly
                  write (333,*) i,j,-rz,rmolz
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
                      if (r2.gt.cutele2) go to 210

                      q=qi*ptchg(j)


c chen
                        if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  q = q * tmp
c mauro - mute 
                        if (mutaid(i).eq.1 .and.
     &                      mutaid(j).eq.1)then
                            q = q/lambda
                        else if (mutaid(i).eq.2 .and.
     &                           mutaid(j).eq.2) then
                            q = q/(1.0d0-lambda)
                        endif

                      s2=1.0d0/r2
                      s = dsqrt(s2)
                      e2 = q*s
                      df2 = -e2*s2

                      rx = df2*rx
                      ry = df2*ry
                      rz = df2*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz
                      e_el = e_el + e2
c mauro - muta <U2-U1>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1) then
                         tmp_e_lambda=tmp_e_lambda+e2/lambda
                   else if (mutaid(i).eq.2 .or. mutaid(j).eq.2) then
                         tmp_e_lambda=tmp_e_lambda-e2/(1.d0-lambda)
                        endif

c virial calculation
                        if (pressON) then
                       virXX = virXX -(xmol(pmol(i))-xmol(pmol(j)))*rx
                       virYY = virYY -(ymol(pmol(i))-ymol(pmol(j)))*ry
                       virZZ = virZZ -(zmol(pmol(i))-zmol(pmol(j)))*rz
                         virial = virial -
     1                   (xmol(pmol(i))-xmol(pmol(j)))*rx -
     2                   (ymol(pmol(i))-ymol(pmol(j)))*ry -
     3                   (zmol(pmol(i))-zmol(pmol(j)))*rz
                       rmolx = xmol(pmol(i))-xmol(pmol(j))
                       rmoly = ymol(pmol(i))-ymol(pmol(j))
                       rmolz = zmol(pmol(i))-zmol(pmol(j))
                  write (333,*) 'i j df drmol'
                  write (333,*) i,j,-rx,rmolx
                  write (333,*) i,j,-ry,rmoly
                  write (333,*) i,j,-rz,rmolz
                        endif

210           continue
              end if
              jbeg3 = jend3 + 1

                dpot(1,i) = dpot(1,i) + dxi
                dpot(2,i) = dpot(2,i) + dyi
                dpot(3,i) = dpot(3,i) + dzi
400     continue
c        if (vp_flag) call vp_cdie()

        return 
        end






