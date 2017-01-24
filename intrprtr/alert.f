        subroutine alert(func,length1,message,length2,level)
c
c       SEND AN ALERT MESSAGE WHEN ERRORS ARE DETECTED
c
c       FUNC is the name of the subroutine with the error
c       MESSAGE is explanation of the error
c
        include 'COMMON/UNITS.BLOCK'
        integer length1,length2,level
        character*(*) func,message

        if (level.gt.0) then
         write(stdo,1)
1        format(//,1x,'******** RED ALERT ********')
        else
         write(stdo,15)
15       format(//,1x,'******** YELLOW ALERT ********')
        end if

        write(stdo,2)func,message
2       format(1x,'******** SURPRISE IN ROUTINE ',a,/,1x,
     1'*** DESCRIPTION: ',a,//)
        if (level.gt.0) stop
        return
        end
