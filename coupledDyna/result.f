	program result
c  This is just a unnecessarily too complicated  program which  
c  computes the rmsd between two given coordinates of a same
c  protein...
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/SHAKE.BLOCK'
	include 'COMMON/MSHAKE.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/DYNA.BLOCK'
	include 'COMMON/SWITCH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/SYMM.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONSTRAN.BLOCK'
	include 'COMMON/TETHER.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/RESTART.BLOCK'
	include 'COMMON/SSBP.BLOCK'
	include 'COMMON/EWALD.BLOCK'
	include 'COMMON/METAL.BLOCK'
	include 'COMMON/SGB.BLOCK'
	include 'COMMON/EBALL.BLOCK'
	include 'COMMON/COUPLED_DYNA.BLOCK'
	include 'mpif.h'

	character*8 name
	integer namel 

	logical failure

	integer j,k,l,iat1,iat2,istep
	integer kk,mat_idx(10*maxshak)

	double precision tmp,gradf
	double precision mat_val(10*maxshak)

	integer iii
	double precision rr_wfly,rx_wfly,ry_wfly,rz_wfly
        logical w_stop(maxmono)

	name = 'result'
        namel= 6

c
c     Here a parallel run is initialized
c
c
	call MPI_INIT(MPIerr)
C	write (*,*) 'Error', MPIerr
	if (MPIerr.ne.0) then
            write(*,*)' *** Error while starting MPI'
            call alert(name,namel,' Error while starting MPI',25,1)
        end if
	call MPI_COMM_RANK(MPI_COMM_WORLD,procID,MPIerr)
	write (*,*) 'Error', MPIerr
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,MPIerr)
	write (*,*) 'Error', MPIerr
	write (*,*) 'numprocs', numprocs
	write (*,*) 'procID', procID
	if (numprocs.ne.2) then
	    write(*,*)' *** Started ',numprocs,'processes'
            write(*,*)' *** but expecting... ', 2
            call alert(name,namel,
     &		' Not correct # of procceses started',35,1)
	end if
	 
c initialize variables
c     proportion of homology energy
	E_homology_const = 10
	call init_dyna()

        call compare_structures()

C  Shut down MPI
	call MPI_FINALIZE(MPIerr)	   
	
	end

