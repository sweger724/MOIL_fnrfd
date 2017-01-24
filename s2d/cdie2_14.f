      subroutine cdie2_14()
c
C second derivatives for the non-bonded pairs of the special 1-4 
c bonded  list 
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/GETVEC.BLOCK'

      double precision rx,ry,rz,r2,s2
      double precision s,s6
      double precision s8,s3,abq1,abq2,add
      integer index,jndex,kndex
c
      integer i,j,k

      do 400 k=ifirst,ilast
            i = spec1(k)
            j = spec2(k)
            rx = coor(1,i) - coor(1,j)
            ry = coor(2,i) - coor(2,j)
            rz = coor(3,i) - coor(3,j)
            r2=rx*rx+ry*ry+rz*rz
            s2=1.0d0/r2
            s6=s2*s2*s2
            s8=s2*s6
            s = dsqrt(s2)
            s3 = s2*s
            abq1=168.d0*p14(1,k)*s8*s8-48.d0*p14(2,k)*s8*s2
     1           +3.d0*p14(3,k)*s3*s2
            abq2=-12.d0*p14(1,k)*s6*s8+6.d0*p14(2,k)*s8-p14(3,k)*s3
c
            index = 6*(i-1)+1
            jndex = 6*(j-1)+1
            kndex = 6*(k-ifirst)+1
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
400   continue

      return 
      end




