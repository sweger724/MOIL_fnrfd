      subroutine bond2_pcs()
c
c calculate second derivatives of the bond potential energy
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c
c
      double precision rx,ry,rz,r2,s,r,add,a1,a2
      integer mm,i,j,index,jndex,kndex

c rx,ry,rz = distance between particles in that axis
c r2 = distance squared between particles
c s = distance between particles
c add, a1,a2 = temporary storage use din the calculations

c      write (6,*) 'BONDPCS: i have ',ifirst,ilast

      do 100 mm=ifirst,ilast
                i=ib1(mm)
                j=ib2(mm)
c                write (6,*) 'BONDPCS LOOP :',mm,i,j
c                write (6,*) 'BONDPCS: COORD1 '
c                write (6,*) coor(1,i),coor(2,i),coor(3,i)
c                write (6,*) 'BONDPCS: COORD2 '
c                write (6,*) coor(1,j),coor(2,j),coor(3,j)
                rx=coor(1,i)-coor(1,j)
                ry=coor(2,i)-coor(2,j)
                rz=coor(3,i)-coor(3,j)
                r2=rx*rx + ry*ry + rz*rz
                r=dsqrt(r2)
                s=1.d0/r
                if (r.lt.0.1) then
                   write (6,*) 'My bonding distance < 0.1.'
                   write (6,*) 'This is problematic for my 2nd deriv.'
                   write (6,*) 'calculations , Please try lowering'
                   write (6,*) 'the TOTAL energy so you get more'
                   write (6,*) 'reasonable structures :) '
                   stop

                end if
                a1=2.d0*kbond(mm)*req(mm)*s*s*s
                a2=2.d0*kbond(mm)*s*(r-req(mm))
c
                index = 6*(i-1)+1
                jndex = 6*(j-1)+1
                kndex = 6*(mm-ifirst)+1
c
                add   = a1*rx*rx + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
c
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*rx
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
c
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*rz*rx
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
c
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*ry + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
c
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*rz
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
c
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*rz*rz + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                offdiag(kndex) = -add
100     continue
      return
      end
