        subroutine comm_init()
        include "mpif.h"
        integer ierror
        call MPI_INIT(ierror)
        if (ierror .ne. MPI_SUCCESS) then
                call alert(comm_init,9,"Fail to init MPI",16,1)
        end if
        return
