        subroutine load_balance(n,my_pe,n_pes,istart,iend)
        implicit none
c
c Given a number n of compute tasks done by independent iterations of some kind.
c give a reasonable partition of n between n_pes processors
c on output n is modified to be the number of compute tasks on "my_pe"
c

        integer n,my_pe,n_pes,istart,iend
c local
        integer j

C if the number of elements to be distributed is zero
c complain and return
c
        if (n.eq.0) then
          write(6,*)' No elements to distribute in load_balance '
          istart = 1
          iend   = 0
          return
        end if
C@
C       write(6,*)' n my_pe n_pes ',n,my_pe,n_pes

c if the number of eleements is smaller than the current my_pe+1 (the first
c pe is zero). The current pe has nothing to do. Set variables appropriately
c so that any con reading will be excuted correctly.
c
        if (n.lt.my_pe+1) then
          istart = n+1
          iend   = n
          n      = 0
          return
        end if

c j is the reminder when n is divided by the number of processing elements
c  (num_pes)
c
        j = n - (n/n_pes)*n_pes

c n is now an estimate for the new number of elements to be computed at
c the present pe.
c
        n = n/n_pes
C@
c       write(*,*)' after some processing j n ',j,n

        if (n.eq.0) then
                write(*,*)' N is very small '
                write(*,*)' N = ... ; n_pes = ... ',n,n_pes
        end if

C and now a correction for the present n due to the reminder
c

        if (j.ne.0 .and. my_pe.lt.j) then
                n = n + 1
        end if
C now find from where to start the reading from con file or another
c full vector representation
c
C@
C       write(*,*)' j my_pe n ',j,my_pe,n
        
        if (j.eq.0) then
                istart = my_pe*n + 1
        else if (j.ne.0 .and. my_pe.lt.j) then
                istart = my_pe*n + 1
        else if (my_pe.ge.j .and. n.ge.1) then
                istart = my_pe*n + j + 1
        else
                write(*,*)' WARNING! istart was not assigned '
        end if
C@
C       write(*,*)' istart = ',istart

c if n eq 0 (i.e. total number of elements smaller than number of pe-s)
c only some of the processors will have a single elements to
c comute. For these pe-s add 1 to istart to iend
c
        if (n.eq.0 .and. my_pe.lt.j) then
         iend = istart + 1
        else
c otherwise iend will be n more elements starting from istart
c
         iend = istart + n - 1
        end if

        return
        end
        

