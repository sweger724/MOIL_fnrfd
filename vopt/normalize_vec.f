
       subroutine normalize_vec(v,length)
       double precision v(*)
       integer length
c
c local
       double precision vnorm
       integer i
       vnorm = 0.d0
       do 1 i=1,length
        vnorm = vnorm + v(i) * v(i)
1       continue
       vnorm = dsqrt (1.d0/vnorm)
       do 2 i=1,length
        v(i) = v(i)*vnorm
2      continue
       return
       end
