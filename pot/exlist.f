        subroutine exlist(l,m,flag)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        integer j,k,l,m,jbeg,jend
        logical flag
        flag = .false.
        jbeg=exc1(l-1)+1
        jend=exc1(l)
        if (jbeg.le.jend) then
                do 100 j=jbeg,jend
                        k=exc2(j)
                        if(k.eq.m) then
                                flag=.true.
                                goto 101
                        end if
100             continue
101             continue
        end if
        return
        end

        subroutine incl14(l,m,flag14)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        integer j,k,l,m,jbeg,jend
        logical flag14
        flag14 = .false.
        jbeg=spec1(l-1)+1
        jend=spec1(l)
        if (jbeg.le.jend) then
                do 100 j=jbeg,jend
                        k=spec2(j)
                        if(k.eq.m) then
                                flag14=.true.
                                goto 101
                        end if
100             continue
101             continue
        end if
        return
        end
