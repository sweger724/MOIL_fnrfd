	subroutine init_mpi_var()
c
c   Initialize all MPI parameters after the data from 
c   the input files are read
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COUPLED_DYNA.BLOCK'
	integer j,k,l

c
c       Initialize the indexes of CA atoms in the array CAindex
c
	do 100 j=1,alignN
	  l=align(j)+1
	  write(*,*) 'J= ',j,' align(j)= ',align(j)
C    up there is +1, because the 1st particle is a N-terminus
C	  

c      Assign indeces of CA atoms to array CAindex()
c
	  if (l.eq.1) then
	    do 105 k=1,poipt(l)
              if (ptnm(k).eq.'CA  ')  then
                CAindex(j) = k
		write(*,*) 'CAindex(',j,')= ',CAindex(j)
              end if
105         continue

	  else
	    do 110 k=poipt(l-1)+1,poipt(l)
	      if (ptnm(k).eq.'CA  ')  then
		CAindex(j) = k
		write(*,*) 'CAindex(',j,')= ',CAindex(j)
	      end if
110	    continue
	  end if
100	continue

	end

