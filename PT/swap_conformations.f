        subroutine swap_conformations(step,accept)
        implicit none

         include 'COMMON/LENGTH.BLOCK'
         include 'COMMON/ENERGY.BLOCK'
         include 'COMMON/DYNA.BLOCK'
         include 'COMMON/CONNECT.BLOCK'
         include 'COMMON/COORD.BLOCK'
         include 'COMMON/VELOC.BLOCK'
         include 'COMMON/UNITS.BLOCK'
         include 'COMMON/PT.BLOCK'
         include 'COMMON/PARALLEL.BLOCK'
         
         double precision e_rep(2),t_rep(2) 
         double precision tmpC(3,maxpt),tmpV(3,maxpt)
         integer neighbor,irep, step, ji, ll, npt3 
         logical accept

         npt3 = 3 * npt
         accept = .false.

         if(random_index.gt.randomSize-4) 
     &         call prepare_random_numbers(step)

         random_index = random_index + 1
         irep=int(RANDOM_ARRAY(random_index)*(Nreplicas-1))
         !write(6,*)"RRR:",random_index,irep,replicaID
         random_index = random_index + 1  ! for deciding the acceptance/rejection 

         if (replicaID .ne. irep .and. replicaID.ne. (irep+1)) return


        if (replicaID.eq.irep)  then
          if (prll_on_off) call reduce_energies()
          e_rep(1)=e_total
          t_rep(1)=tempi(1)

          neighbor = procID + num_pes 
          !write(6,*)"comunicating",procID,neighbor
          call Recv_Double(e_rep(2),1,neighbor,1221)
          call Send_Double(e_rep(1),1,neighbor,1231)

          call Recv_Double(t_rep(2),1,neighbor,122)
          call Send_Double(t_rep(1),1,neighbor,123)
        endif 

        
        if(replicaID.eq.(irep+1)) then 
          if (prll_on_off) call reduce_energies()
          e_rep(2)=e_total
          t_rep(2)=tempi(1)
      
          neighbor = procID - num_pes 
          !write(6,*)"comunicating",procID,neighbor
          call Send_Double(e_rep(2),1,neighbor,1221)
          call Recv_Double(e_rep(1),1,neighbor,1231)

          call Send_Double(t_rep(2),1,neighbor,122)
          call Recv_Double(t_rep(1),1,neighbor,123)
        endif  

       ! make the decision about the exchange
       call ptempering(e_rep,t_rep,accept)
       
       ! exchange the coordinates and momenta
       !write(6,*)"AAA",accept
       if (accept) then
       
         do ji=1,npt
           do ll=1,3
             tmpC(ll,ji) = coor(ll,ji)
             tmpV(ll,ji) = velo(ll,ji)
           end do
         end do 
       
         if (replicaID .eq. irep) then

           neighbor = procID + num_pes
           call Recv_Double(coor(1,1),npt3,neighbor,124)
           call Recv_Double(velo(1,1),npt3,neighbor,125)
           
           call Send_Double(tmpC(1,1),npt3,neighbor,126)
           call Send_Double(tmpV(1,1),npt3,neighbor,127)
           
           !write(6,*)"RRR",sqrt(t_rep(1)/t_rep(2))
           do ji=1,npt
             do ll=1,3
               velo(ll,ji) = velo(ll,ji) * sqrt(t_rep(1)/t_rep(2))
             end do
           end do
           
         end if

         if (replicaID .eq. irep+1) then
           neighbor = procID - num_pes 
           call Send_Double(tmpC(1,1),npt3,neighbor,124)
           call Send_Double(tmpV(1,1),npt3,neighbor,125)

           call Recv_Double(coor(1,1),npt3,neighbor,126)
           call Recv_Double(velo(1,1),npt3,neighbor,127)

           !write(6,*)"RRR",sqrt(t_rep(1)/t_rep(2))
           do ji=1,npt
             do ll=1,3
               velo(ll,ji) = velo(ll,ji) * sqrt(t_rep(2)/t_rep(1))
             end do
           end do 
           
         end if
       
       end if  ! accept

c     now we report the acceptance ratio
      
      if ( replicaID.eq.irep .and. my_pe .eq. 0 ) then          
        sumswap1 = sumswap1 + 1
        if (accept) sum_accept1 = sum_accept1 + 1
        write(stdo,1222) step,t_rep(1),t_rep(2),sum_accept1/sumswap1 
      endif    

      if ( replicaID.eq.irep+1 .and. my_pe .eq. 0 ) then
        sumswap2 = sumswap2 + 1
        if (accept) sum_accept2 = sum_accept2 + 1
        write(stdo,1222) step,t_rep(1),t_rep(2),sum_accept2/sumswap2
      endif 

 1222 format('  Acc. ratio in steps=',i9,' temps = ',2x,f6.2,'<-->',
     >f6.2,2x,f4.2)
     
      end
