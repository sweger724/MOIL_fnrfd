      subroutine ebond()
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
c mauro add for muta
      include 'COMMON/MUTA.BLOCK'

c calculate bond energies and forces.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision rx,ry,rz,r2,s,r,db,df,e
      integer mm,i,j

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
      do 100 mm=bonds_start,bonds_end

                i = ib1(mm)
                j = ib2(mm)
                rx=coor(1,i)-coor(1,j)
                ry=coor(2,i)-coor(2,j)
                rz=coor(3,i)-coor(3,j)
                r2=rx*rx + ry*ry + rz*rz
                s=dsqrt(r2)
                r=2.d0/s
                db=s-req(mm)
                df=kbond(mm)*db
                e=df*db
c@@@
		if (e.gt.1.e2) then
              write(6,*)'bond',i,j,s,req(mm),e
                end if
c@@@

c mauro - <U1-U2>l
                        if (mutaid(i).eq.1 .or. mutaid(j).eq.1) then
                           tmp_e_lambda=tmp_e_lambda+e/lambda
                        endif
                        if (mutaid(i).eq.2 .or. mutaid(j).eq.2) then
                           tmp_e_lambda=tmp_e_lambda-e/(1.0d0-lambda)
                        endif
c end mauro - <U1-U2>l

                e_bond=e_bond + e
                df=df*r
                dpot(1,i) = dpot(1,i) + rx*df
                dpot(2,i) = dpot(2,i) + ry*df
                dpot(3,i) = dpot(3,i) + rz*df
                dpot(1,j) = dpot(1,j) - rx*df
                dpot(2,j) = dpot(2,j) - ry*df
                dpot(3,j) = dpot(3,j) - rz*df

c               by assumption no bonds between slow and fast

100     continue

      return
      end
