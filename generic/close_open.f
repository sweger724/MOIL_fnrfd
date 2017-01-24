        subroutine close_open(unit)
c
c a subroutine that close, rewind and reopen a file with unit
c number - unit. The file is opened with status unknown.
c
        integer unit
c local
        character *80 file_name
        logical exists
        inquire(unit=unit,name=file_name,opened=exists)

        if (.not.exists) then
         call alert('close_open',9,'File does not exist',20,1)
         return
        end if

        close(unit)
        open(unit=unit,file=file_name,status='unknown')
        rewind unit

        return
        end
