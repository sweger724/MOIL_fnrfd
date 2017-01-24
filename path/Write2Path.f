        subroutine Write2Path(ufile,filename,length,rr,npt,toOpen)
          
          implicit none

          integer length,npt,ufile
          character*80 filename
          double precision rr(3,*)
          logical toOpen

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/COORD.BLOCK'
          include 'COMMON/ENERGY.BLOCK'
          include 'COMMON/PATH.BLOCK'
          include 'COMMON/MASSFAC.BLOCK'
          include 'COMMON/PARALLEL.BLOCK'
c Local
          integer i,j,jloc,l,n, amount
          double precision rsave(3,maxpt*lgrid),etmp(lgrid)

           if (paral) then

               if (first) then
C Sending and immediately writing to save memory
C First save coordinates in first

                 if (toOpen) then
                    ufile=555
                    open(unit=ufile,status='unknown',
     &                   form='unformatted',file=filename(1:length))
                 endif
                 do i=2,pseg+1
                   do j = 1,npt
                      jloc = (i-1)*npt+j
                      do l = 1,3
                         rsave(l,jloc) = rr(l,jloc)
                      end do
                   end do
                 end do

                 do i = 1, pseg+1
                   call Write2File(ufile,npt,rr(1,(i-1)*npt+1),
     &                              e0(i),massfac,MASSWEIGHT)
                 end do

              endif  ! first

        do n=1, numprocs-1

           If (procID.eq.n) then
              Call Send_Double(rr(1,npt+1),3*npt*pseg,0,procID)
              Call Send_Double(e0(2),pseg,0,2*procID)
           endif

           If (first) then
              if (n.lt.mod(igrid-2,numprocs)) then
                  amount = (igrid-2)/numprocs +1
              else
                  amount = (igrid-2)/numprocs
              endif
              Call Recv_Double(rr(1,npt+1),3*npt*amount,n,n)
              Call Recv_Double(etmp(2),amount,n,2*n)

              do 51 i=2,amount+1

                call Write2File(ufile,npt,rr(1,(i-1)*npt+1)
     $                           ,etmp(i),massfac,MASSWEIGHT)

51            continue

           endif  ! first

        enddo !  enddo over numprocs
        
C Restore back my original coordinates
C and write the very last structure to uwipth

       if (last)
     >       Call Send_Double(r_final(1,1),3*npt,0,1123)
             Call Send_Double(e0(pseg+2),1,0,1124)
         if (first) then
             Call Recv_Double(r_final(1,1),3*npt,numprocs-1,1123)
             Call Recv_Double(etmp(1),1,numprocs-1,1124)

             call Write2File(ufile,npt,r_final(1,1)
     $                        ,etmp(1),massfac,MASSWEIGHT)

             if (toOpen) close(ufile)
             do i=2,pseg+1
               do j = 1,npt
                  jloc = (i-1)*npt+j
                  do l = 1,3
                     rr(l,jloc) = rsave(l,jloc)
                  end do
                end do
             end do
          endif
c this is the end of the new commuincations stuff
C If not paral

          else
            if (toOpen) then
               ufile=555
               open(unit=ufile,status='unknown',form='unformatted',
     1              file=filename(1:length))
            endif

            do 511 i=1,pseg+2
                call Write2File(ufile,npt,rr(1,(i-1)*npt+1),e0(i)
     1                           ,massfac,MASSWEIGHT)
 511        continue

            if (toOpen) close(ufile)

         end if

       return
       end 
