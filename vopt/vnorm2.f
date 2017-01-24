       double precision function vnorm2(v,length)
       double precision v(*)
       integer length
c
c local
       integer i
       vnorm2 = 0.d0
       do 1 i=1,length
        vnorm2 = vnorm2 + v(i) * v(i)
1       continue
        vnorm2 = dsqrt(vnorm2)
       return
       end
