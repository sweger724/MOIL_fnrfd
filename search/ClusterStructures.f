      subroutine ClusterStructures(Csize,oldsize,count,expandNode)
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SEARCH.BLOCK'

      integer Csize, oldsize, count, expandNode

      integer relevant, oldRelevant, jump, itt
      double precision rms,sum,diff(3,maxpt),minRms
      double precision Pt(maxNeighbors,count)

      integer i,j,k,kk,l,loop,picked(maxNeighbors),minI

      relevant = 0

      if (count .eq. 0) then
        neighbors(expandnode) = 0
      end if

      do i =1,Csize
       call Distance(centers(1,npt*(i-1)+1),
     &               centers(1,npt*(expandNode-1)+1),rms)
        if (rms .lt. 2.d0*cutoff .and. relevant .lt. maxNeighbors) then
          relevant = relevant + 1
          picked(relevant) = i
        end if
        if (i .eq. oldsize) oldRelevant = relevant
      end do

      do loop=1,10

      do i=1,count
        sum = 0.d0
        do j = 1,relevant
          call Distance(structures(1,(i-1)*npt+1)
     &                 ,centers(1,(picked(j)-1)*npt+1),rms)
          Pt(j,i) = exp(-(rms/(0.25*cutoff))**2)
          sum = sum + Pt(j,i)
        end do

        do j=1,relevant
          Pt(j,i) = Pt(j,i)/sum
        end do

      end do

C     Find new centers
      do j = oldRelevant+1,relevant
        sum = 0
        kk = (picked(j)-1)*npt
        do k = 1, npt
          coor(1,k) = 0.0
          coor(2,k) = 0.0
          coor(3,k) = 0.0
        end do

        do i = 1, count
         
!    align the center with the structure
          call Distance(structures(1,(i-1)*npt+1)
     &                 ,centers(1,(picked(j)-1)*npt+1),rms)
          do k = 1, npt
            coor(1,k) = coor(1,k) + Pt(j,i)*structures(1,npt*(i-1)+k)
            coor(2,k) = coor(2,k) + Pt(j,i)*structures(2,npt*(i-1)+k)
            coor(3,k) = coor(3,k) + Pt(j,i)*structures(3,npt*(i-1)+k)
          end do

          sum = sum + Pt(j,i)
        end do

        do k = 1, npt
          coor(1,k) = coor(1,k) / sum 
          coor(2,k) = coor(2,k) / sum
          coor(3,k) = coor(3,k) / sum
        end do 

        ! find the closest structure to the median
        minRms = 10000
        do i = 1, count
           call Distance(structures(1,1+(i-1)*npt),coor(1,1),rms)
           if (rms .lt. minRms) then
             minRms = rms
             minI = i
           end if
        end do

        do k = 1, npt
          do l = 1,3
            coor(l,k) = structures(l,k+(minI-1)*npt)
            diff(l,k) = (centers(l,kk+k) - coor(l,k))/100
          end do
        end do

        do l = 1,relevant
          if (l .ne. j) then
            call Distance(coor(1,1),centers(1,(picked(l)-1)*npt+1),rms)
            itt = 0
            do while ( rms .lt. cutoff)
              ! move the center to avoid clash with other center
              itt = itt + 1
              do k = 1, npt
                coor(1,k) = coor(1,k) + diff(1,k)
                coor(2,k) = coor(2,k) + diff(2,k)
                coor(3,k) = coor(3,k) + diff(3,k)
              end do
             call Distance(coor(1,1),centers(1,(picked(l)-1)*npt+1),rms)
              if (itt .eq. 101) then
C                  if (loop .eq. 10)write(6,*)"ERROR!!!", rms
                  rms = cutoff
C                  stop
              end if
            end do

          end if
        end do

        ! find the closest structure to the new center
        minRms = 10000
        do i = 1, count
           call Distance(structures(1,1+(i-1)*npt),coor(1,1),rms)
           if (rms .lt. minRms) then
             minRms = rms
             minI = i
           end if
        end do

        ! update the j-th center
        do k = 1, npt
          centers(1,kk+k) = structures(1,k+(minI-1)*npt)
          centers(2,kk+k) = structures(2,k+(minI-1)*npt)
          centers(3,kk+k) = structures(3,k+(minI-1)*npt)
        end do

      end do ! j

      call removeSmallCenters(Csize,oldsize,count,expandNode,relevant,
     &                        oldRelevant,picked,Pt)

      end do  ! loop

      jump = 0
      do j = 1, relevant
        if ( picked(j) .ne. expandNode) then
      conform(picked(j)) = conform(picked(j))*P(expandNode,j-jump)*count
          if (P(expandNode,j-jump) .gt. 0.001) then
            write(6,'(a,i4,1x,i4,1x,f5.3)')
     &          'P: ',expandNode,picked(j),P(expandNode,j-jump)
          end if
        else
          jump = 1
        end if
      end do
      neighbors(expandnode) = relevant - 1

      return
      end
