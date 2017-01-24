       double precision function vnorm(v,length)
       double precision v(*)
       integer length
c
c local
       integer i
       vnorm = 0.d0
       do 1 i=1,length
        vnorm = vnorm + v(i) * v(i)
1       continue
       vnorm = dsqrt (vnorm/length)
       return
       end
