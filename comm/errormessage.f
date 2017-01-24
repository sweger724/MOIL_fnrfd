        subroutine errormessage(error,message)
        integer error, nlen, ierr
        character*(*) message
czva        include '../include/Tcommlib_f.h'
czva        include '../include/Trtlib_f.h'
c
        include 'mpif.h'
c
        nlen=80
czva    call tcomm lib error message(error,message)

c        call MPI_Error_string(error,message,nlen,ierr)
c
        return
        end


