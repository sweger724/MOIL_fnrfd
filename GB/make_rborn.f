
      subroutine make_rborn
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/CONNECT.BLOCK'

      integer i,j,bondtonum
      character atom,bondto
c      write (6,*) 'Setting up GBSA Parameter set from moil types '

c YS
      Fsmax = 0.0
      firstbool = 1
      stepno = 0
      do i = 1,npt

         atom = ptnm(i)(1:1)

         if (atom.eq.'H'.or.atom.eq.'1'.or.atom.eq.'2') then
c         if (atom.eq."H") then
            fs(i) = 0.85d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)

C     Now lets find what this H is bonded to (this is a bit lazy but we only
C     call this once so..

            bondtonum = -1
            do j = 1,nb
               if (ib1(j).eq.i) then
                  bondtonum = ib2(j)
               end if
               if (ib2(j).eq.i) then
                  bondtonum = ib1(j)
               end if
            end do

            if (bondtonum.eq.-1) then
               bondto = 'O'
            else
               bondto = ptnm(bondtonum)(1:1)
            end if
            if (bondto.eq.'O') then
               rborn(i) = 0.8d0
            else if (bondto.eq.'N') then
               rborn(i) = 1.2d0
            else if (bondto.eq.'C') then
               rborn(i) = 1.3d0
            else 
               rborn(i) = 1.2d0
            end if

         else if (atom.eq.'C') then
            fs(i) = 0.72d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)

            rborn(i) = 1.70d0
C            rborn(i) = 0.3

         else if (atom.eq.'N') then
            fs(i) = 0.79d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.55d0
         else if (atom.eq.'O') then
            fs(i) = 0.85d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.50d0
         else if (atom.eq.'F') then
            fs(i) = 0.88d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.47d0
         else if (atom.eq.'P') then
            fs(i) = 0.86d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.85d0
         else if (atom.eq.'S') then 
            fs(i) = 0.96d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.80d0
C United methyl is like carbon?
         else if (atom.eq.'M') then
            fs(i) = 0.72d0
            if ( fs(i) .gt. Fsmax ) Fsmax = fs(i)
            rborn(i) = 1.70d0
          
         else 
         write (6,*) 'BAD ATOM SYMBOL ',i,atom
         end if
      end do
c      do i = 1,npt
c         write (6,*) 'rborn param ',fs(i),rborn(i)
c      end do


      return
      end
