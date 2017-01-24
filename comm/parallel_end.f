        subroutine parallel_end()

        integer ierr

        call MPI_FINALIZE(ierr)

        return
        end
