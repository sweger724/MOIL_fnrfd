      subroutine etether()
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/UNITS.BLOCK'

c calculate tether energies and forces.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      double precision rx,ry,rz,tf2,r2
      integer mm,i

c e_tether=total tether energy 
c rx,ry,rz = distance squared between particles in that axis
c r2 = distance squared between particles
c s = distance between particles
c db,df temporary variables
c e = tether energy (non-acumulated)
c also used : ichunk for vectorizacion length
c xtmp,ytmp,ztemp for vectorization purposes


c initialize e_tether
      e_tether=0.d0
c initialize loop over tethers in ichunk chunks
      do 100 mm=1,n_tether

      		i=tether(mm)
      		rx=coor(1,i)-coor2(1,i)
      		ry=coor(2,i)-coor2(2,i)
      		rz=coor(3,i)-coor2(3,i)
      		r2=rx*rx + ry*ry + rz*rz

		e_tether = e_tether + tf(mm)*r2

		tf2 = 2.d0*tf(mm)

		dpot(1,i) = dpot(1,i) + rx*tf2
		dpot(2,i) = dpot(2,i) + ry*tf2
		dpot(3,i) = dpot(3,i) + rz*tf2

100	continue


      return
      end
