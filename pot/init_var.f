        subroutine init_var()
c
c  subroutine to initialize variables
c
        include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/LINE.BLOCK'

        pi180 = 180.d0/(4.d0*datan(1.d0))
	silent = .false.

        return
        end

