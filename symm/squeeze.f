        subroutine squeeze()
c
c a subroutine to squeeze lattice molecules that escaped from
c the primitive unit cell (taken here cubic) back to their place.
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/SYMM.BLOCK'
c
c local
        logical first
        double precision totmas,xcm,ycm,zcm,xnorm
        double precision xlow,xhigh,ylow,yhigh,zlow,zhigh
        integer i,k,istart
        logical cm
        data first/.true./

        save first, xlow, xhigh, ylow, yhigh, zlow, zhigh

        if (a.lt.1.d0) then
                write(*,*) " In Squeeze, a=",a
                stop
        end if
c       totmas = 0.d0
c       xcm    = 0.d0
c       ycm    = 0.d0
c       zcm    = 0.d0
        if (first) then
         xlow   = -a*0.5d0
         xhigh  = a*0.5d0
         ylow   = -b*0.5d0
         yhigh  = b*0.5d0
         zlow   = -c*0.5d0
         zhigh  = c*0.5d0
         first  = .false.
        end if

c       do 1 i=1,npt
c        totmas = totmas + ptms(i)
c        xcm    = xcm + ptms(i)*coor(1,i)
c        ycm    = ycm + ptms(i)*coor(2,i)
c        zcm    = zcm + ptms(i)*coor(3,i)
c1      continue

c       totmas = 1.d0/totmas
c       xcm = xcm*totmas
c       ycm = ycm*totmas
c       zcm = zcm*totmas

c       do 2 i=1,npt
c        coor(1,i) = coor(1,i) - xcm
c        coor(2,i) = coor(2,i) - ycm
c        coor(3,i) = coor(3,i) - zcm
c2      continue

c
c currently water molecules and selected ions are transformed
c back to the box. Decisions are made according to the
c position of the oxygen of the water or the ion.
c
        do 10 i=1,totmon
         if (moname(i).eq.'DMPC'.or.moname(i).eq.'DOPC'
     &       .or.moname(i).eq.'NATA') then
                 cm = .true.
         else
                 cm = .false.
         end if
         if (moname(i).eq.'TIP3' .or. moname(i).eq.'NAAT'
     1  .or. moname(i).eq.'NAO' .or. moname(i).eq.'CL'
     2  .or. moname(i).eq.'MG' .or. moname(i).eq.'CMO'
     3  .or. moname(i).eq.'CO' .or. moname(i).eq.'SUL'
     3  .or.moname(i)(1:3).eq.'SPC' .or. moname(i).eq.'K'
     4  .or.moname(i).eq.'AR'.or. cm) then
          if (i.eq.1) then
           istart = 1
          else
           istart = poipt(i-1)+1
          end if
           if (cm) then
                xcm = 0.d0
                ycm = 0.d0
                zcm = 0.d0
                do 45 k=istart,poipt(i)
                        xcm = xcm + coor(1,k)
                        ycm = ycm + coor(2,k)
                        zcm = zcm + coor(3,k)
45              continue
                xnorm = 1.d0/dfloat(poipt(i)-istart+1)
                xcm = xcm*xnorm
                ycm = ycm*xnorm
                zcm = zcm*xnorm
           else
                xcm = coor(1,istart)
                ycm = coor(2,istart)
                zcm = coor(3,istart)
           end if

           if (.not.metalyes) then
            if (ycm.lt.ylow) then
             do 5 k=istart,poipt(i)
              coor(2,k) = coor(2,k) + b
5            continue
c            write(stdo,*)' Monomer ',i,'  has been Y- transformed '
            end if

            if (ycm.gt.yhigh) then
             do 8 k=istart,poipt(i)
              coor(2,k) = coor(2,k) - b
8            continue
c            write(stdo,*)' Monomer ',i,'  has been Y+ transformed '
            end if
           end if

           if (xcm.lt.xlow) then
            do 4 k=istart,poipt(i)
             coor(1,k) = coor(1,k) + a
4           continue
c           write(stdo,*)' Monomer ',i,'  has been X- transformed '
           end if

           if (zcm.lt.zlow) then
            do 6 k=istart,poipt(i)
             coor(3,k) = coor(3,k) + c
6           continue
c           write(stdo,*)' Monomer ',i,'  has been Z- transformed '
           end if

           if (xcm.gt.xhigh) then
            do 7 k=istart,poipt(i)
             coor(1,k) = coor(1,k) - a
7           continue
c           write(stdo,*)' Monomer ',i,'  has been X+ transformed '
           end if

           if (zcm.gt.zhigh) then
            do 9 k=istart,poipt(i)
             coor(3,k) = coor(3,k) - c
9           continue
c           write(stdo,*)' Monomer ',i,'  has been Z+ transformed '
           end if
         end if
10      continue
        return
        end

