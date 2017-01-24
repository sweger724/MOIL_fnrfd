      subroutine cdie2_wat1()

C second derivatives for the non-bonded pairs of the first non-bonded 
c list (both electrostatic and Van der Walls)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SYMM.BLOCK'
c
      double precision rx,ry,rz,r2,s2,ai,bi,qi
      double precision s,s6
      double precision s8,s3,abq1,abq2,add
      double precision xi,yi,zi
      integer symidx
      integer index,jndex,kndex
c
      integer ii,i,j
c
      if (usesym) then

         do 100 ii=ifirst,ilast
c
           i=pairs1(ii)
           j=pairs2(ii)
           symidx=pairs1s(ii)
c
           if (rotop(symidx)) then
              xi = arot*coor(1,i) + symop(1,symidx)*a
              yi = brot*coor(2,i) + symop(2,symidx)*b
              zi = crot*coor(3,i) + symop(3,symidx)*c
           else
              xi = coor(1,i) + symop(1,symidx)*a
              yi = coor(2,i) + symop(2,symidx)*b
              zi = coor(3,i) + symop(3,symidx)*c
           end if

           rx = xi - coor(1,j)
           ry = yi - coor(2,j)
           rz = zi - coor(3,j)
           r2=rx*rx+ry*ry+rz*rz
c           write(6,*)' dist. ',sqrt(r2)
           qi= q_pair(ii)
           s2=1.0d0/r2
           s = dsqrt(s2)
           s3 = s2*s

           if (tmp_pair(ii).eq.1) then

              if ((r2.gt.cutvdw2).and.(.not.nocut)) then
                 ai = 0.d0
                 bi = 0.d0
              else 
                 ai= a_pair(ii)
                 bi= b_pair(ii)
              end if
              s6=s2*s2*s2
              s8=s2*s6
              abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1               +3.d0*qi*s3*s2
              abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3

           else

              abq1=3.d0*qi*s3*s2
              abq2=-qi*s3

           end if

c           write(*,*)'#1 abq1 abq2 ',abq1,abq2
c

c NOTE - rotations not supported

           index = 6*(i-1)+1
           jndex = 6*(j-1)+1
           kndex = 6*(ii-ifirst)+1
           add         = abq1*rx*rx + abq2
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add         = abq1*ry*rx
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c     
           add         = abq1*rz*rx 
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c     
           add        = abq1*ry*ry + abq2
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add        = abq1*rz*ry
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add        = abq1*rz*rz + abq2
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
           
c
 100     continue

      else

         do 200 ii=ifirst,ilast
c
           i=pairs1(ii)
           j=pairs2(ii)
c
           rx = coor(1,i) - coor(1,j)
           ry = coor(2,i) - coor(2,j)
           rz = coor(3,i) - coor(3,j)
           r2=rx*rx+ry*ry+rz*rz

           qi= q_pair(ii)
           s2=1.0d0/r2
           s = dsqrt(s2)
           s3 = s2*s
c           write(6,*)'wat1  i j r e ',i,j,sqrt(r2),s*qi
           if (tmp_pair(ii).eq.1) then

              if ((r2.gt.cutvdw2).and.(.not.nocut)) then
                 ai = 0.d0
                 bi = 0.d0
              else 
                 ai= a_pair(ii)
                 bi= b_pair(ii)
              end if
              s6=s2*s2*s2
              s8=s2*s6
              abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1               +3.d0*qi*s3*s2
              abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3

           else

              abq1=3.d0*qi*s3*s2
              abq2=-qi*s3

           end if

c           write(*,*)'#1 abq1 abq2 ',abq1,abq2
c
           index = 6*(i-1)+1
           jndex = 6*(j-1)+1
           kndex = 6*(ii-ifirst)+1
           add         = abq1*rx*rx + abq2
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add         = abq1*ry*rx
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add         = abq1*rz*rx 
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add        = abq1*ry*ry + abq2
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add        = abq1*rz*ry
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
           add        = abq1*rz*rz + abq2
           index = index + 1
           jndex = jndex + 1
           kndex = kndex + 1
           diag(index) = diag(index) + add
           diag(jndex) = diag(jndex) + add
           offdiag(kndex) = -add
c
 200     continue

      end if   
         
      return 
      end
c-------------------------------------------------------------------
      subroutine cdie2_wat2()
c
C second derivatives for the non-bonded pairs of the third non-bonded 
c list which includes particles with upper cutoff cutele2 - cutoff 
c appropriate for electrostics - 
c (and lower cutoff - cutvdw2 - the lower cutoff electrostatic was 
c already calculated using the van der Waals (first) loop for charged
c particles). 
c Includes ONLY electrostic forces
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SYMM.BLOCK'

      double precision rx,ry,rz,r2,s2,qi
      double precision s
      double precision s3,abq1,abq2,add
      double precision xi,yi,zi
      integer symidx
      integer index,jndex,kndex
c
      integer i,j,ii,iik
c
c
      if (usesym) then
         
         ii=ifirst-1
         do 100 iik=ifirst,ilast
c
           ii=ii+1 
           i=pairs1(ii)
           j=pairs2(ii)
           symidx=pairs1s(ii)
c
           if (rotop(symidx)) then
              xi = arot*coor(1,i) + symop(1,symidx)*a
              yi = brot*coor(2,i) + symop(2,symidx)*b
              zi = crot*coor(3,i) + symop(3,symidx)*c
           else
              xi = coor(1,i) + symop(1,symidx)*a
              yi = coor(2,i) + symop(2,symidx)*b
              zi = coor(3,i) + symop(3,symidx)*c
           end if

           rx = xi - coor(1,j)
           ry = yi - coor(2,j)
           rz = zi - coor(3,j)
           r2=rx*rx+ry*ry+rz*rz
c           skip the pair of waters if oxygens too far
            if ((tmp_pair(ii).eq.1) .and.
     &           (r2.gt.cutele2) .and. (.not.nocut)) then
                  ii=ii+8
                  go to 100
            end if
            s2=1.0d0/r2
            s = dsqrt(s2)
            s3 = s2*s
            abq1=3.d0*qi*s3*s2
            abq2=-qi*s3
c            write(*,*)'#3 abq1 abq2 ',abq1,abq2
c
            index = 6*(i-1)+1
            jndex = 6*(j-1)+1
            kndex = 6*(ii-ifirst)+1
            add         = abq1*rx*rx + abq2
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*ry*rx
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*rx 
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*ry*ry + abq2
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*ry
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*rz + abq2
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
 100     continue

      else

         ii=ifirst-1
         do 110 iik=ifirst,ilast
c
            ii=ii+1
            i=pairs1(ii)
            j=pairs2(ii)
c
            rx = coor(1,i) - coor(1,j)
            ry = coor(2,i) - coor(2,j)
            rz = coor(3,i) - coor(3,j)
            r2=rx*rx+ry*ry+rz*rz
c           skip the pair of waters if oxygens too far
            if ((tmp_pair(ii).eq.1) .and.
     &           (r2.gt.cutele2) .and. (.not.nocut)) then
                  ii=ii+8
                  go to 110
            end if
            qi=q_pair(ii)   
            s2=1.0d0/r2
            s = dsqrt(s2)
            s3 = s2*s
            abq1=3.d0*qi*s3*s2
            abq2=-qi*s3
c            write(*,*)'#3 abq1 abq2 ',abq1,abq2
c
            index = 6*(i-1)+1
            jndex = 6*(j-1)+1
            kndex = 6*(ii-ifirst)+1
            add         = abq1*rx*rx + abq2
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*ry*rx
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*rx 
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*ry*ry + abq2
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*ry
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
            add         = abq1*rz*rz + abq2
            index = index + 1
            jndex = jndex + 1
            kndex = kndex + 1
            diag(index) = diag(index) + add
            diag(jndex) = diag(jndex) + add
            offdiag(kndex) = -add
c
 110     continue

      end if

      return 
      end










