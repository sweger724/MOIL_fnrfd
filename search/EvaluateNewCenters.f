        subroutine  EvaluateNewCenters(Csize,oldsize,closed,expandNode)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/PATH2.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/DYNA.BLOCK'

        integer Csize, oldsize, closed, expandNode
        integer FindNeighbour, i, k,ki,kj,j

        double precision tmp, hhMin, hh

        do i = oldsize+1,Csize
          k = FindNeighbour(expandNode,i)
          if (k .eq. -1) then 
            g(i) = 100000000
            write(6,*)"warning: neighbor NOT found",expandNode,i
           else 
            g(i) = g(expandNode) - dlog(P(expandNode,k))
          end if
        end do

        do i = closed+2, oldsize
          k = FindNeighbour(expandNode,ord(i))
          if (k .ne. -1) then
            tmp = g(expandNode) - dlog(P(expandNode,k));
            if (tmp < g(ord(i))) g(ord(i)) = tmp
          end if
        end do

         do i = closed+2,Csize
          
          hhMin = 100000
          do j =1, Csize
            if (nodeType(j) .ne. nodeType(ord(i))) then
              ki = (ord(i)-1)*npt+1
              kj = (j-1)*npt+1
              call Distance(centers(1,ki),centers(1,kj),rms)
              if (rms .lt. cutoff) rms = cutoff
                hh = - 50.d0 * (rms-cutoff)/cutoff * dlog(0.15d0) +g(j);
              if (hh .lt. hhMin) hhMin = hh
            end if
          end do
          
          h(ord(i)) = hhMin
          
          if (crdstyl .eq. "PATH") h(ord(i)) = 0.d0

        end do
!        h(2) = 100000

        do i = closed+2, Csize
          f(ord(i)) = g(ord(i)) + h(ord(i))
        end do

        return
        end

        function FindNeighbour(from,to)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/SEARCH.BLOCK'

        integer FindNeighbour,a,b,c, from,to

        a = 1

        do while (neighbors(from).ge.a .and. node(from,a).lt.to) 
          a = a + 1
        end do
        
        if (node(from,a) .eq. to ) then 
           FindNeighbour = a
          else 
           FindNeighbour = -1
        end if

        return

C        a = 1
C        b = maxNeighbours/2
C        c = maxNeighbours

C        while (node(from,b) .ne. to ) do
C          if (node(from,b) .gt. to) then 
C            b = a + (b-a)/2
C            c = b
C          else if (node(from,b) .lt. to) then
C            b = b + (c-b)/2
C            a = b
C          else
C            ! found
C          end if
C        end do

C        return

        end
