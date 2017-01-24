c -------------------------------------------------------------------------
c Gather the full list of shaked bonds in each processor

	subroutine gather_shake(nshak,ishak1,ishak2,dissq,shak_diff)

	integer ishak1(*),ishak2(*),nshak
	double precision dissq(*),shak_diff(*)

           return
           end

c ************************************************************
	subroutine gather_shake_balance (shake_bal)

	integer shake_bal(*)

	return
	end

c *******************************************************************	
c get the shared particles from each processor
	subroutine get_shared_pt(shak_l,count,pts_num,
     *                    RBpts,p_dissq,p_shak_diff)

	integer shak_l(2,*)
	double precision RBpts(3,*)
	double precision p_dissq(*),p_shak_diff(*)
	integer pts_num,count

	return
	end


	subroutine gatherRB(pts,velocs,mass,my_pts,my_vel,
     *                                   pts_num,nptsh,displace)

	integer pts_num,nptsh
	integer displace
	double precision pts(3,*),velocs(3,*),mass(*)
	double precision my_pts(3,*),my_vel(3,*)

	return
	end


	subroutine gather_int (my_val,list)

	integer my_val,list(*)

	return 
	end
	
	subroutine gather_int_pairs(my_array,array,displace,nptsh)

	integer array(2,*),my_array(2,*)
	integer nptsh,displace
	
	return
	end

	subroutine gather_double(array,displace,nptsh)

	double precision array(2,*)
	integer nptsh,displace

	return
	end

	subroutine gather_double1(array,displace,nptsh)

	double precision array(*)
	integer nptsh,displace

	return
	end

	subroutine gather_velo(array,displace,nptsh)

	double precision array(3,*)
	integer nptsh,displace

	return
	end
  
	subroutine gather_int1(array,displace,nptsh)

	integer array(*)
	integer nptsh,displace

	return
	end


c Get integer array from the right neighbor, and send to the left.
	
	subroutine send_right_int (send,recv,size)

	integer send(*),recv(*),size

	return 
	end

	subroutine send_right_double (send,recv,size)

	double precision send(*),recv(*)
	integer size

	return 
	end

	subroutine send_left_int (send,recv,size)

	integer send(*),recv(*),size

	return 
	end

	subroutine send_left_double (send,recv,size)

	double precision send(*),recv(*)
	integer size

	return 
	end


	subroutine gather_coord (crd,npt,pt_start,pt_end,disp,nptsh)
	
	double precision crd(3,*)
	integer npt,pt_start,pt_end,disp,nptsh

	return
	end
