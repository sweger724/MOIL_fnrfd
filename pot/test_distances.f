        subroutine test_dis()

C       test_distances: If a distance in the neighbour
C       list is shorter than 1.5A print a warning message
C

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

        double precision rx,ry,rz,r2
        double precision xi,yi,zi

        integer i,j,k,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
        integer ptbeg,ptend


        if (nocut) then
                cutvdw2 = 10000.d0
                cutele2 = 10000.d0
        end if


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
c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
c


        do 400 i=ptbeg,ptend


                jend1 = point1(i)
                xi  = coor(1,i)
                yi  = coor(2,i)
                zi  = coor(3,i)

                if (jbeg1.le.jend1) then

                do 200 k=jbeg1,jend1
                        j=list1(k)
                        rx = xi - coor(1,j)
                        ry = yi - coor(2,j)
                        rz = zi - coor(3,j)
                        r2=rx*rx+ry*ry+rz*rz
                        if (r2.lt.2.25d0) then
                 write(6,*)'************'
                 write(6,1000)i,ptnm(i),j,ptnm(j),dsqrt(r2)
1000     format(1x,' ptid ',i6,' ptnm ',a4,' ptid ',i6,' ptnm ',a4,/
     1          ,1x,' Distance : ',f8.5,' too short! ')
                        end if
200     continue

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
                        if (r2.lt.2.25d0) then
                         write(6,*)'************'
                         write(6,1000)i,ptnm(i),j,ptnm(j),dsqrt(r2)
                        end if
205     continue

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
                        if (r2.lt.2.25d0) then
                         write(6,*)'************'
                         write(6,1000)i,ptnm(i),j,ptnm(j),dsqrt(r2)
                        end if
210           continue

              end if
              jbeg3 = jend3 + 1

400     continue

c
c test also 1-4 distances. 
c
        do 600 i=1,totspe
                j = spec1(i)
                k = spec2(i)
                rx = coor(1,j)-coor(1,k)
                ry = coor(2,j)-coor(2,k)
                rz = coor(3,j)-coor(3,k)
                r2 = rx*rx + ry*ry + rz*rz
                if (r2.lt.2.25d0) then
                        write(6,*)' 1-4 1-4 1-4 1-4 '
                        write(6,1000)i,ptnm(i),j,ptnm(j),dsqrt(r2)
                end if
600     continue

        return 
        end
