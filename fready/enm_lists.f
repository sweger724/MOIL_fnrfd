      subroutine enm_lists(r_initial,r_final)
      
      implicit none
            
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ELASTIC.BLOCK'

c      local variables

      integer i,j,counter
      double precision rx,ry,rz,r2,cutoff2
      double precision r_initial(3,maxpt),r_final(3,maxpt)

      cutoff2= enm_cutoff*enm_cutoff
      counter =0
      g1list1(0)=0 
      g2list1(0)=0

c
c.....generate lists for the first network
c       
        write(6,*) "I'm here!",enm_cutoff
      do 100 i = 1, npt-1
         do 200 j =i+1,npt
            
            rx = r_initial(1,i) - r_initial(1,j)
            ry = r_initial(2,i) - r_initial(2,j)
            rz = r_initial(3,i) - r_initial(3,j)
            
            r2 = rx*rx + ry*ry + rz*rz
            
            if (r2.lt.cutoff2) then
                counter = counter + 1
                g1list2(counter) = j
                req1(counter) = dsqrt(r2)
C                write(6,*)'Pair i,j',i,j,req1(counter)
            endif
            
200      continue
         
         g1list1(i)=counter
c         write(6,*)'g1list1(i) ',i,g1list1(i) 
100   continue
      write(6,*)'Contacts in structure1: ',counter,enm_cutoff
      counter=0
c
c.....generate lists for the second network
c
      do 300 i = 1, npt-1
         do 400 j =i+1,npt

            rx = r_final(1,i) - r_final(1,j)
            ry = r_final(2,i) - r_final(2,j)
            rz = r_final(3,i) - r_final(3,j)

            r2 = rx*rx + ry*ry + rz*rz

            if (r2.lt.cutoff2) then
                counter = counter + 1
                g2list2(counter) = j
                req2(counter) = dsqrt(r2)
            endif

400      continue

         g2list1(i)=counter

300   continue
      write(6,*)'Contacts in structure2: ',counter,enm_cutoff
      end
