	subroutine ener14()


      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'

	double precision rx,ry,rz,r2,s2,a,b,e1,e2,df,df1,df2
	double precision s,s6

	integer i,j,k


	e_vdw14 = 0.0d0
	e_el14  = 0.0d0


		 do 100 i=1,totspe
			j = spec1(i)
			k = spec2(i)
			rx = coor(1,j) - coor(1,k)
			ry = coor(2,j) - coor(2,k)
			rz = coor(3,j) - coor(3,k)
			r2=rx*rx+ry*ry+rz*rz
			s2=1.0d0/r2

			if (p14(1,i).gt.1.d-12) then
			 s6=s2*s2*s2
                         a = p14(1,i)*s6*s6
                         b = p14(2,i)*s6
                         e1 = (a - b)
                         df1 = -6.0d0*s2*(a+e1)
			else
			 e1  = 0.d0
			 df1 = 0.d0
			end if

			if (dabs(p14(3,i)).gt.1.d-12) then
			 s = dsqrt(s2)
                         e2 = p14(3,i)*s
                         df2 = -e2*s2
			else
			 e2  = 0.d0
			 df2 = 0.d0
			end if
			 
cc Luca for virial (vdW only: take only df1)
c here j,k pair is as i,j in cdie_ew
c NB SIGN: dpot(j) = -df*rx => force(j) = df*rx = f_ji
c    here:      k                    k              kj
c = -f_ij, => f_ij = -df*rx
c      jk       jk

                        virvdw5 = virvdw5 - (df1*rx*rx + df1*ry*ry
     &                             + df1*rz*rz)
cc

                        df = df1 + df2
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        dpot(1,j) = dpot(1,j) + rx
                        dpot(2,j) = dpot(2,j) + ry
                        dpot(3,j) = dpot(3,j) + rz
                        dpot(1,k) = dpot(1,k) - rx
                        dpot(2,k) = dpot(2,k) - ry
                        dpot(3,k) = dpot(3,k) - rz
			e_vdw14 = e_vdw14 + e1 
			e_el14 = e_el14 + e2

100		 continue

	return 
	end
