	subroutine add_v()
c	The updates velocities vectors/vx/y/z




        integer i,error,length
        character*80 message

	return
	end	




c-----------------------------------------------------------------------------
	subroutine update_ex_shared()
c	The updates the ex_shared data structure




        integer i,error,length
        character*80 message

	return
	end	
c---------------------------------------------------------------------------
	subroutine shak_unhappy()
c	Tell everybody that I failed to converge




        integer i,error,length,my_pe_p1
        character*80 message

	return
	end	
c---------------------------------------------------------------------------
	subroutine shak_happy()
c	Tell everybody that I failed to converge



        integer i,error,length,my_pe_p1
        character*80 message

	return
	end	

c---------------------------------------------------------------------------
	subroutine get_nshakes(nshakes)
	
	integer nshakes(*)




        integer i,error,length,my_pe_p1
        character*80 message

	return
	end	
c---------------------------------------------------------------------------
	subroutine get_ishakes(nshakes, g_ishak1, g_ishak2)
	
	integer nshakes(*), g_ishak1(*), g_ishak2(*)


        integer i,error,length,my_pe_p1, dis,temp
        character*80 message

	return
	end	

		
c---------------------------------------------------------------------------
	subroutine get_dissqs(nshakes, g_dissq, g_shak_diff)
	
	integer nshakes(*)
	double precision g_dissq(*), g_shak_diff(*)


        integer i,error,length, dis
        character*80 message

	end	

		
	
