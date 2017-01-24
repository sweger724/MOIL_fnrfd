       subroutine SortCenters(Csize,closed)

       implicit none

       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/SEARCH.BLOCK'

        integer min,Csize,closed,i,j,k1,k2,k, l
        double precision tmp

        ! sort centers by f values
        do i = 1, Csize
          if (ord(i).ge. closed+2) ord(i) = Csize +1
        end do
        
        do i = closed+2, Csize
          k=1
          
          do while (ord(k).lt.i)
            k = k + 1
          end do

          min = k
          
          do j = 1 , Csize
            if (ord(j).ge. i) then
               if (f(j).lt.f(min)) min = j
             end if
          end do
          
          ! exchange i and min
C            tmp = ord(min)
            ord(min) = i
C            ord(i) = tmp
          
        end do

       return
       end

       subroutine WriteState(Csize,closed)
        implicit none

       include 'COMMON/LENGTH.BLOCK'
       include 'COMMON/PATH.BLOCK'
       include 'COMMON/CONNECT.BLOCK'
       include 'COMMON/SEARCH.BLOCK'

        integer min,Csize,closed,i,j,k1,k2,k, l
        double precision rms
        
          write(6,*)"========="
          write(6,*)"              f          g        h      rmsd
     &     conf.  ord T"
          do j = 1 , Csize

            do i = 1, Csize
              if (ord(i) .eq. j) then
                call Distance(centers(1,(i-1)*npt+1),r_final(1,1),rms)
       write(6,'(a,i4,1x,f9.3,1x,f9.3,1x,f9.3,1x,f6.3,
     &           1x,f9.3,i5,i2)')
     &            "Node ",i,f(i),g(i),h(i), rms, conform(i)
     &            ,ord(i),nodeType(i)
                goto 13
              end if
            end do
13          continue
          end do 
          
C          do i = 1 , Csize-1
C            do j = i+1 , Csize
C              call Distance(centers(1,(i-1)*npt+1)
C     &                     ,centers(1,(j-1)*npt+1),rms)
C              if(rms .lt. cutoff) write(6,*) "DDD:",i,j,rms
C            end do
C          end do

C          do i=1, closed
C            write(6,*)"========="
C            write(6,*)"Transitions from node", i
C            do j=1, neighbors(i)
C              write(6,*)" ",node(i,j), P(i,j), dlog(P(i,j))
C            end do
C          end do
C          
C          write(6,*)"========="
C         return
       end
