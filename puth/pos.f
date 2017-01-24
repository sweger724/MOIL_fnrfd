       subroutine pos(success,ih)
       implicit none
c
c Placing hydrogens
c We made use of known coordinates and the data stored in the connectivity block
c Only bond lengths and bond angles are used (i.e. NO torsions are employed)
c In case of missing data, the direction are picked at random
c (e.g. water molecule)
c No information on possible hydrogen bonds candidate is taken into
c account in this admittedly oversimplifed placement
c
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        integer ih 
        logical success

        character*3 name
        integer namel

        integer level
        integer i,j,nbondj
        integer ibj(maxptbn)
        double precision reqhj
        double precision cos1,cos2,cos3
        double precision sin1,sin2,sin3
        double precision gamma,delta
        double precision e2(3), e4(3)
        real e1(3), e3(3)
        double precision pi

        pi     = 4.d0*datan(1.d0)
        debug  = .false.
        name   = 'pos'
        namel  = 3
        success = .false.
c
c Find out a bonded particle to the neig particle (a neighbor to the pt)
c It is assumed that hydrogen has only one bond (i.e. no Boron)
c
        do 1 i=1,nb_all
         if (ib1(i).eq.ih) then
          if (coor(1,ib2(i)).gt.9998.d0) then
           level = 3
           write(*,*)' ib1 ib2 coor ',ib1(i),ib2(i),coor(1,ib2(i))
           call alert(name,namel,'H bonded to undefined crd',25,level)
           return
          end if
          j     = ib2(i)
          reqhj = req(i)
          go to 2
         else if (ib2(i).eq.ih) then
          if (coor(1,ib1(i)).gt.9998.d0) then
           level = 3
           write(*,*)' ib1 ib2 coor ',ib1(i),ib2(i),coor(1,ib1(i))
           call alert(name,namel,'H bonded to undefined crd',25,level)
           return
          end if
          j     = ib1(i)
          reqhj = req(i)
          go to 2
         end if
1       continue
        level = 1
        write(6,*) " H bonded to nothing", ih, nb_all
        call alert(name,namel,'H bonded to nothing !?',22,level)
        return
2       continue
c
c Find out all particles bonded to j (excluding the hydrogen)
c count only bonds to a particle with defined coordinates
c
        nbondj = 0
        do 3 i=1,nb_all
         if ((ib1(i).eq.j .and. ib2(i).ne.ih) .or.
     1            (ib2(i).eq.j .and. ib1(i).ne.ih)) then
          if (ib1(i).eq.j) then
           if (coor(1,ib2(i)).gt.9998.d0) go to 3
           nbondj = nbondj + 1
           ibj(nbondj)   = ib2(i)
          else
           if (coor(1,ib1(i)).gt.9998.d0) go to 3
           nbondj = nbondj + 1
           ibj(nbondj)   = ib1(i)
          end if
         end if
3       continue
c
c Now check how many particles are bonded to the j-th particle
c
        if (nbondj.eq.0) then
c generating position for the hydrogen connected to an isolated center
c use correct distance but random direction
c
         call RANLUX(e1,3)
         cos1  = e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3)
         cos1  = 1.d0/dsqrt(cos1)
         e1(1) = e1(1)*cos1
         e1(2) = e1(2)*cos1
         e1(3) = e1(3)*cos1
         coor(1,ih) = coor(1,j) + reqhj*e1(1)
         coor(2,ih) = coor(2,j) + reqhj*e1(2)
         coor(3,ih) = coor(3,j) + reqhj*e1(3)
        else if (nbondj.eq.1) then
c search for the angle correspoding to ih-neigh-ibj(1)
         do 4 i=1,nangl_all
          if (iangl2(i).eq.j) then
            if ((ih.eq.iangl1(i) .and. ibj(1).eq.iangl3(i)) .or.
     1        (ih.eq.iangl3(i) .and. ibj(1).eq.iangl1(i))) then
c get the cosine of the angle between the OTHER bond to j and the bond between
c j and the hydrogen
             if (dabs(angleq(i)-pi).lt.1.d-6
     1        .or. dabs(angleq(i)).lt.1.d-6) then
              level = 1
              write(*,*)' i j k ',iangl1(i),iangl2(i),iangl3(i)
              call alert(name,namel,'LINEAR ANGLE',12,level)
             end if
             cos1    = dabs(cos(angleq(i)))
            else
             goto 4
            end if
            sin1    = dsqrt(1.d0-cos1*cos1)
c generate a unit vector along the OTHER bond
            e2(1)   = coor(1,j) - coor(1,ibj(1))
            e2(2)   = coor(2,j) - coor(2,ibj(1))
            e2(3)   = coor(3,j) - coor(3,ibj(1))
            cos2    = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
            cos2    = 1.d0/dsqrt(cos2)
            e2(1)   = e2(1)*cos2
            e2(2)   = e2(2)*cos2
            e2(3)   = e2(3)*cos2
c generate a vector in a random direction, to initiate sampling in the
c allowed circle
            call RANLUX(e3,3)
            cos3    = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
            cos3    = 1.d0/dsqrt(cos3)
            e3(1)   = e3(1)*cos3
            e3(2)   = e3(2)*cos3
            e3(3)   = e3(3)*cos3
            cos2    = e2(1)*e3(1) +  e2(2)*e3(2) + e2(3)*e3(3)
            if ((dabs(cos2-1.d0).lt.1.d-6) .or.
     1       (dabs(cos2+1.d0).lt.1.d-6)) then
             level = 1
        call alert(name,namel,'Vectors linearly dependent',27,level)
            end if
c orthonormalize e3 with respect to the existing bond direction
            e3(1)   = e3(1) - cos2*e2(1)
            e3(2)   = e3(2) - cos2*e2(2)
            e3(3)   = e3(3) - cos2*e2(3)
            cos3    = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
            cos3    = 1.d0/dsqrt(cos3)
            e3(1)   = e3(1)*cos3
            e3(2)   = e3(2)*cos3
            e3(3)   = e3(3)*cos3
c generate the hydrogen position. the displacement vector
c as compared to the j atom coordinateis a linear combination
c of the two orthonormalized vectors e2 and e3
            coor(1,ih)   = coor(1,j) + reqhj*(sin1*e3(1)+cos1*e2(1))
            coor(2,ih)   = coor(2,j) + reqhj*(sin1*e3(2)+cos1*e2(2))
            coor(3,ih)   = coor(3,j) + reqhj*(sin1*e3(3)+cos1*e2(3))
            go to 5
          end if
4        continue
         level = 1
         call alert(name,namel,'ANGLE NOT FOUND',15,level)
5        continue        
        else if (nbondj.eq.2) then
c only two bonds are used in the present verison. A danger that does not exist
c in the present force field (extended atom model) but does exist
c for all atom case is
c the possibility of inverted improper torsion.
         cos1 = 999.d0
         cos2 = 999.d0
         cos3 = 999.d0
         do 6 i=1,nangl_all
          if (iangl2(i).eq.j) then
           if (dabs(angleq(i)-pi).lt.1.d-6 .or.
     1      dabs(angleq(i)).lt.1.d-6) then
            level = 1
            call alert(name,namel,'LINEAR ANGLE',12,level)
           end if
           if ((iangl1(i).eq.ih .and. iangl3(i).eq.ibj(1)) .or.
     1          (iangl1(i).eq.ibj(1) .and. iangl3(i).eq.ih)) then
                cos1 = cos(angleq(i))
           else if ((iangl1(i).eq.ih .and. iangl3(i).eq.ibj(2)) .or.
     1          (iangl1(i).eq.ibj(2) .and. iangl3(i).eq.ih)) then
                cos2 = cos(angleq(i))
           else if ((iangl1(i).eq.ibj(1) .and. iangl3(i).eq.ibj(2)) .or.
     1          (iangl1(i).eq.ibj(2) .and. iangl3(i).eq.ibj(1))) then
                cos3 = cos(angleq(i))
           end if
          end if
6        continue
         if (cos1.gt.1.d0 .or. cos2.gt.1.d0 .or. cos3.gt.1d0) then
          level = 1
          call alert(name,namel,'Undefined cosine values',23,level)
         end if
         e1(1) = coor(1,ibj(1)) - coor(1,j)
         e1(2) = coor(2,ibj(1)) - coor(2,j)
         e1(3) = coor(3,ibj(1)) - coor(3,j)
c sin1 is used here as a buffer for normalization
         sin1  = e1(1)*e1(1) + e1(2)*e1(2) + e1(3)*e1(3)
         sin1  = 1.d0/dsqrt(sin1)
         e1(1) = e1(1)*sin1
         e1(2) = e1(2)*sin1
         e1(3) = e1(3)*sin1

         e2(1) = coor(1,ibj(2)) - coor(1,j)
         e2(2) = coor(2,ibj(2)) - coor(2,j)
         e2(3) = coor(3,ibj(2)) - coor(3,j)
c sin2 is used here as a buffer for normalization
         sin2  = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
         sin2  = 1.d0/dsqrt(sin2)
         e2(1) = e2(1)*sin2
         e2(2) = e2(2)*sin2
         e2(3) = e2(3)*sin2
c orthonormalize e2 with respect to e1
         sin1  = e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)
         e2(1) = e2(1) - sin1*e1(1)
         e2(2) = e2(2) - sin1*e1(2)
         e2(3) = e2(3) - sin1*e1(3)
c normalize the new e2
         sin2  = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
         sin2  = 1.d0/dsqrt(sin2)
         e2(1) = e2(1)*sin2
         e2(2) = e2(2)*sin2
         e2(3) = e2(3)*sin2
c generate the third unit vector by a vector products of e1 and e2
c
         e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
         e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
         e3(3) = e1(1)*e2(2) - e1(2)*e2(1)
c check that e3 is normalized
         sin3 = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
         if (dabs(sin3).lt.1.d-6) then
          level = 1
          call alert(name,namel,'Fishy vector product',20,level)
         end if
c A unit vector along the bond between j and ih is expanded in terms
c of e1,e2 and e3 as follows - R(ih)-R(j) = cos1 e1 + gamma e2 + delta e3
c where gamma=(cos2 -cos1*cos3)/sqrt(1-cos3^2) and delta=sqrt(1-gamma^2-cos1^2)
         gamma = (cos2-cos3*cos1)/dsqrt(1.d0 - cos3*cos3)
         delta = gamma*gamma+cos1*cos1
         if (1.d0-delta.lt.1.d-10) then
          delta = 0.d0
         else
          delta = dsqrt(1.d0-delta)
         end if
         coor(1,ih)=coor(1,j)+reqhj*(cos1*e1(1)+gamma*e2(1)+delta*e3(1))
         coor(2,ih)=coor(2,j)+reqhj*(cos1*e1(2)+gamma*e2(2)+delta*e3(2))
         coor(3,ih)=coor(3,j)+reqhj*(cos1*e1(3)+gamma*e2(3)+delta*e3(3))

        else if (nbondj.eq.3) then

c
c calculate eigenvectors along the exisiting bonds pointing to the central
c bond.
c
         e1(1) = coor(1,j) - coor(1,ibj(1))
         e1(2) = coor(2,j) - coor(2,ibj(1))
         e1(3) = coor(3,j) - coor(3,ibj(1))

         e2(1) = coor(1,j) - coor(1,ibj(2))
         e2(2) = coor(2,j) - coor(2,ibj(2))
         e2(3) = coor(3,j) - coor(3,ibj(2))

         e3(1) = coor(1,j) - coor(1,ibj(3))
         e3(2) = coor(2,j) - coor(2,ibj(3))
         e3(3) = coor(3,j) - coor(3,ibj(3))

c calculate a unit vector in the direction of the missing bond. It is
c assumed that it is in the direction of the average of the e1 to e3
c eigenvectors which (hopefully) is not too bad
c
         e4(1) = e1(1) + e2(1) + e3(1)
         e4(2) = e1(2) + e2(2) + e3(2)
         e4(3) = e1(3) + e2(3) + e3(3)
         delta = 1.d0/dsqrt(e4(1)*e4(1)+e4(2)*e4(2)+e4(3)*e4(3))
         e4(1) = e4(1)*delta
         e4(2) = e4(2)*delta
         e4(3) = e4(3)*delta
         coor(1,ih) = coor(1,j) + reqhj*e4(1)
         coor(2,ih) = coor(2,j) + reqhj*e4(2)
         coor(3,ih) = coor(3,j) + reqhj*e4(3)

        else
         write(stdo,*)' Nbondj = ',nbondj
         call alert('pos',3,' Impossible # of bonds to H',26,1)
        end if
        success = .true.
        return
        end
