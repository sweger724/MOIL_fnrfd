      subroutine eCG_bond()
      
      implicit none
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c calculate bond energies and forces.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision rx,ry,rz,r2,s,r,db,dE,E
      integer mm,i,j, T

c e_bond=total bond energy 
c rx,ry,rz = distance squared between particles in that axis
c r2 = distance squared between particles
c s = distance between particles
c db,df temporary variables
c e = bond energy (non-acumulated)
c also used : ichunk for vectorizacion length
c xtmp,ytmp,ztemp for vectorization purposes



c initialize eb
      e_bond=0.d0
c initialize loop over bonds in ichunk chunks
      do 100 mm=1,nb

                i = ib1(mm)
                j = ib2(mm)
                rx=coor(1,i)-coor(1,j)
                ry=coor(2,i)-coor(2,j)
                rz=coor(3,i)-coor(3,j)
                r2=rx*rx + ry*ry + rz*rz
                r=dsqrt(r2)
                T = bondType(mm)
C               write(6,*)"bond style: ",T      
                  call eCG_WellBondPot(r,T,E,dE)
                  if(E .gt. E_CG_max) then
                    write(6,*)"Residue: ",poimon(i),poimon(j)
                  endif
c@@@
C               write(6,*)' bond between i j r ebij ',i,j,2.d0/r,e
c@@@
                
                e_bond = e_bond + E
                dE = dE/r
                dpot(1,i) = dpot(1,i) + rx*dE
                dpot(2,i) = dpot(2,i) + ry*dE
                dpot(3,i) = dpot(3,i) + rz*dE
                dpot(1,j) = dpot(1,j) - rx*dE
                dpot(2,j) = dpot(2,j) - ry*dE
                dpot(3,j) = dpot(3,j) - rz*dE

c               by assumption no bonds between slow and fast

100     continue

      return
      end
