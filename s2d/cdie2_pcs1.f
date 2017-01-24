      subroutine cdie2_pcs1(sepfast)

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
      double precision s,s6,b6,b12
      double precision s8,s3,abq1,abq2,add
      double precision xi,yi,zi
      integer symidx
      integer index,jndex,kndex
      logical sepfast
c
      integer ii,i,j
c
      if (usesym) then

         do 100 ii=ifirst,ilast
c
           i=pairs1(ii)
           j=pairs2(ii)
           if (sepfast) then
              if ((i.gt.npts).and.(j.gt.npts)) go to 100 
           end if
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
           if ((r2.gt.cutvdw2).and.(.not.nocut)) then
               ai = 0.d0
               bi = 0.d0
            else 
               ai= a_pair(ii)
               bi= b_pair(ii)
           end if
           qi= q_pair(ii)
           if (lesid(i).ne.0) then
              if (lesid(i).eq.lesid(j)) then
                 ai = ai*tmp_pair(ii)
                 bi = bi*tmp_pair(ii)
                 qi = qi*tmp_pair(ii)
              end if
           end if
           s2=1.0d0/r2
           s6=s2*s2*s2
           s8=s2*s6
           s = dsqrt(s2)
           s3 = s2*s
           if (arith) then
           b6 = bi*bi*bi*bi*bi*bi
           b12 = b6*b6
           abq1=168.d0*ai*b12*s8*s8-48.d0*ai*b6*s8*s2
     1               +3.d0*qi*s3*s2
           abq2=-12.d0*ai*b12*s6*s8+6.d0*b6*ai*s8-qi*s3
           else
           abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1               +3.d0*qi*s3*s2
           abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3
           endif
c           write(*,*)'#1 abq1 abq2 ',abq1,abq2
c

c NOTE - rotations not supported

           if (sepfast) then
c
c eliminate slow-fast contributions - assumption: both i and j
c cannot be larger than npts
c
             if (i.gt.npts) then
                index = 6*(j-1)+1
                add         = abq1*rx*rx + abq2
                diag(index) = diag(index) + add
                add         = abq1*ry*rx
                diag(index+1) = diag(index+1) + add
                add         = abq1*rz*rx
                diag(index+2) = diag(index+2) + add
                add         = abq1*ry*ry + abq2
                diag(index+3) = diag(index+3) + add
                add         = abq1*rz*ry
                diag(index+4) = diag(index+4) + add
                add         = abq1*rz*rz + abq2
                diag(index+5) = diag(index+5) + add
                go to 100
             end if
             if (j.gt.npts) then
                index = 6*(i-1)+1
                add         = abq1*rx*rx + abq2
                diag(index) = diag(index) + add
                add         = abq1*ry*rx
                diag(index+1) = diag(index+1) + add
                add         = abq1*rz*rx
                diag(index+2) = diag(index+2) + add
                add         = abq1*ry*ry + abq2
                diag(index+3) = diag(index+3) + add
                add         = abq1*rz*ry
                diag(index+4) = diag(index+4) + add
                add         = abq1*rz*rz + abq2
                diag(index+5) = diag(index+5) + add
                go to 100
             end if

           end if

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
           if ((r2.gt.cutvdw2).and.(.not.nocut)) then
               ai = 0.d0
               bi = 0.d0
            else 
               ai= a_pair(ii)
               bi= b_pair(ii)
           end if
           qi= q_pair(ii)
           if (lesid(i).ne.0) then
              if (lesid(i).eq.lesid(j)) then
                 ai = ai*tmp_pair(ii)
                 bi = bi*tmp_pair(ii)
                 qi = qi*tmp_pair(ii)
              end if
           end if
           s2=1.0d0/r2
           s6=s2*s2*s2
           s8=s2*s6
           s = dsqrt(s2)
           s3 = s2*s
           if (arith) then
           b6 = bi*bi*bi*bi*bi*bi
           b12 = b6*b6
           abq1=168.d0*ai*b12*s8*s8-48.d0*ai*b6*s8*s2
     1               +3.d0*qi*s3*s2
           abq2=-12.d0*ai*b12*s6*s8+6.d0*b6*ai*s8-qi*s3
           else
           abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1               +3.d0*qi*s3*s2
           abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3
           endif
c           write(*,*)'#1 abq1 abq2 ',abq1,abq2
c

           if (sepfast) then
c
c eliminate slow-fast contributions - assumption: both i and j
c cannot be larger than npts
c
             if (i.gt.npts) then
                index = 6*(j-1)+1
                add         = abq1*rx*rx + abq2
                diag(index) = diag(index) + add
                add         = abq1*ry*rx
                diag(index+1) = diag(index+1) + add
                add         = abq1*rz*rx
                diag(index+2) = diag(index+2) + add
                add         = abq1*ry*ry + abq2
                diag(index+3) = diag(index+3) + add
                add         = abq1*rz*ry
                diag(index+4) = diag(index+4) + add
                add         = abq1*rz*rz + abq2
                diag(index+5) = diag(index+5) + add
                go to 200
             end if
             if (j.gt.npts) then
                index = 6*(i-1)+1
                add         = abq1*rx*rx + abq2
                diag(index) = diag(index) + add
                add         = abq1*ry*rx
                diag(index+1) = diag(index+1) + add
                add         = abq1*rz*rx
                diag(index+2) = diag(index+2) + add
                add         = abq1*ry*ry + abq2
                diag(index+3) = diag(index+3) + add
                add         = abq1*rz*ry
                diag(index+4) = diag(index+4) + add
                add         = abq1*rz*rz + abq2
                diag(index+5) = diag(index+5) + add
                go to 200
             end if

           end if

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
