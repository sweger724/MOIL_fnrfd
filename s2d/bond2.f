      subroutine bond2()
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SCNDRV.BLOCK'

c calculate second derivatives of the bond potential energy


      double precision rx,ry,rz,r2,s,r,add,a1,a2
      integer mm,i,j,index,jndex,kndex

c e_bond=total bond energy 
c rx,ry,rz = distance between particles in that axis
c r2 = distance squared between particles
c s = distance between particles
c add, a1,a2 = temporary storage use din the calculations



      do 100 mm=1,nb
                i=ib1(mm)
                j=ib2(mm)
                rx=coor(1,i)-coor(1,j)
                ry=coor(2,i)-coor(2,j)
                rz=coor(3,i)-coor(3,j)
                r2=rx*rx + ry*ry + rz*rz
                r=dsqrt(r2)
                s=1.d0/r
                a1=2.d0*kbond(mm)*req(mm)*s*s*s
                a2=2.d0*kbond(mm)*s*(r-req(mm))
                index = 6*(i-1)+1
                jndex = 6*(j-1)+1
                kndex = 6*(mm-1)+1
                add   = a1*rx*rx + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*rx
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*rz*rx
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*ry + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*ry*rz
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
                index = index + 1
                jndex = jndex + 1
                kndex = kndex + 1
                add   = a1*rz*rz + a2
                diag(index)   = diag(index) + add
                diag(jndex)   = diag(jndex) + add
                d2bond(kndex) = -add
100     continue
      return
      end
