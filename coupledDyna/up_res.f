c A subroutine for updating restart files.
	subroutine 	up_res(step)
	include 'COMMON/RESTART.BLOCK'	
	integer step
	
 	if (urst_crd .eq. -1) return	

c  	UPDATE THE RESTART FILES
 	  call putcrd(urst_crd,'CHARM')
	  call close_open(urst_crd)
c	A parameter  file for restart currently only step #
 	   write(urst,*) step
 	   call close_open(urst)

	   return
	   end

