      subroutine GetDave(davg ,d0)

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      
      double precision davg, d0(*)
      integer i,logN,log2

      davg = 0.d0
      if (.not.paral) then
        do i = 1,pseg+1
            davg  = davg + d0(i)
        enddo
        davg = davg / (pseg+1)
        return
      endif

C Parallel algorithm here:
C
      do i = 1, pseg+1
         if (.not.last) then
           if (i.le.pseg) davg  = davg + d0(i)
         else 
            davg  = davg + d0(i)
         endif
      enddo


C    compute log2(numprocs)
      logN=log2(numprocs)

      call Gather_Data(logN,davg,numprocs)

C   compute actuall average in processor 0
      if (first) davg = davg /(igrid-1)

C   redistribute computed average
      call Distribute_Data(logN,davg,numprocs)

      return
      end


C****************************************************************************
      function log2(a)
        implicit none

        integer a,log2

        integer tmpi

        log2=0
        tmpi=a
        do while (tmpi.gt.1)
           if (mod(tmpi,2).eq.1) tmpi=tmpi+1
           tmpi=tmpi/2
           log2 = log2 +1
        end do
        
        return
      end


C****************************************************************************
      subroutine Gather_Data(Nloops,Cdata,proc)

        implicit none

        integer Nloops,proc
        double precision Cdata 
      
        integer loop2,loop,i,rem,rem1,from,to
        
        loop2=1
        
        do loop = 1, Nloops
          loop2 = loop2 * 2

          rem = 0
          if (mod(proc,loop2).ne.0) rem =1

          do i = 0, proc/loop2 + rem -2
            from = (2*i +1)*loop2/2
            to = i * loop2
            call Transmit(from,to,Cdata,.false.)
          end do

          rem1 = 0
          if (mod(proc,loop2/2).ne.0) rem1 =1

          if ( (2*(proc/loop2 + rem)) .eq.
     &      (proc/(loop2/2) + rem1) ) then

            from = (2*(proc/loop2 + rem)-1)*loop2/2
            to = ((proc/loop2 + rem) -1) * loop2
            call Transmit(from,to,Cdata,.false.)

          endif
        end do 
        
        return
      end


C****************************************************************************
      subroutine Distribute_Data(Nloops,Cdata,proc)
        
        implicit none
        
        integer Nloops,proc
        double precision Cdata
        
        integer loop2,loop,i,rem,rem1,from,to
        
        loop2 = 2**(Nloops+1)
        
        do loop = Nloops, 1, -1

          loop2 = loop2 / 2

          rem = 0
          if (mod(proc,loop2).ne.0) rem =1

          do i = 0, proc/loop2 + rem -2
            to = (2*i +1)*loop2/2
            from = i * loop2
            call Transmit(from,to,Cdata,.true.)
          end do

          rem1 = 0
          if (mod(proc,loop2/2).ne.0) rem1 =1

          if ( (2*(proc/loop2 + rem)) .eq.
     &       (proc/(loop2/2) + rem1) ) then

            to = (2*(proc/loop2 + rem)-1)*loop2/2
            from = ((proc/loop2 + rem) -1) * loop2
            call Transmit(from,to,Cdata,.true.)

          endif
        end do 

        return
      end
