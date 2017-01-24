       subroutine deselect(r,pointr,nselec)

       implicit none
       double precision r(3,*)
       integer nselec,pointr(*),i,j,k

       do 1 i = 1,nselec
          j = pointr(i)
          do 2 k = 1,3
             r(k,i) = r(k,j)
 2          continue
 1       continue
       
       return
       end
