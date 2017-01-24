        function number()
c
c number return the first floating point number in a line (separated
c by spaces from previous expressions) and erase the number from the
c line
c
        character*12 legal
        double precision number 

        integer i,j,k,l,level
        logical success

        include 'COMMON/LINE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        data legal/'1234567890-.'/

        success = .false.
        do 2 i=0,nexp-1
                if (i.eq.0) then
                        k=1
                        l=point(1)
                else
                        k=point(i)+1
                        l=point(i+1)
                end if
                if (debug) then
                  write(stdo,*)' k l line (k:l) ',k,l,line(k:l)
                end if
                do 1 j=1,12
                 if (line(k:k).eq.legal(j:j)) then
                        rewind (jnkf)
                        write(jnkf,*)line(k:l)
                        rewind (jnkf)
                        read(jnkf,*,err=999)number
                        success = .true.
                        line(k:l) = ' '
                        return
                 end if
1               continue
2       continue
        if (.not. success) then
         number = 0.d0
         level = 1
         call alert('number',6,' No number found',16,level)
         return
        end if
999     continue
        number = 0.d0
        level = 1
        call alert('number',6,'Format error in number number',29,level)
        return
        end
