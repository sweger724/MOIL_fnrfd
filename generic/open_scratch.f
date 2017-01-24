        subroutine open_scratch()
        implicit none
        include 'COMMON/LINE.BLOCK'

        jnkf = 25
        open (unit=jnkf,status='scratch')
        return 
        end

       
