        subroutine metal_wall()
c
c calculating metal "wall" as flat exponential repulsion from +/- B/
c along the y axis (in accord with symmetry, that must apply)
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
	  include 'COMMON/METAL.BLOCK'

        integer i
        double precision y,y6,grad_v,q_e,y_start,a0_6
	  logical first
	  
	  data first/.true./

	  save grad_v, y_start, a0_6, first

	  if (first) then
	  	grad_v  = v_elec/b_wall
	  	y_start = 0.5d0*b_wall
		a0_6   = 6.d0*a0_metal
		do 10 i=1,npt
		 if (coor(2,i).lt.-y_start) then
			call alert('wall',4,'Wall too high',13,1)
		 else if (coor(2,i).gt.y_start) then
			write(*,*)' i y y_start ',i,coor(2,i),y_start
			call alert('wall',4,'Wall too low',13,1)
		 end if
10		continue
		first = .false.
	  end if

	  if (dabs(grad_v).gt.1.d-10) then
         do 1 i=1,npt
		if (flagchr(i)) then
c
c Electrode plate part
c
		    q_e       = ptchg(i)*grad_v
		    e_lctd    = e_lctd + q_e*(coor(2,i) + y_start)
		    dpot(2,i) = dpot(2,i) + q_e
		end if
1	    continue
	   end if
c
c Disperssion replusion part
c

	  do 2 i=1,npt
		    if (coor(2,i).gt.0.d0) then
                 y   = coor(2,i)-y_start
		    else
		     y   = coor(2,i)+y_start
		    end if
		    y6 = y*y
		    y6 = y6*y6*y6
                e_wall    = e_wall + a0_metal/y6
                dpot(2,i) = dpot(2,i) - a0_6/(y6*y)
2       continue

        return
        end

