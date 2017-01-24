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
         
         double precision e_rep(2),t_rep(2),tmp1,tmp2 
         integer irep, step, ji, ll, npt3
         logical accept

         character*18 name
         integer namel

         character*80 comm_file 
         integer ucomm, flen
         logical exist


         name='swap_conformations'
         namel=18
         flen = 21
         npt3 = 3 * npt
         accept = .false.

         if(random_index.gt.randomSize-4) 
     &       call prepare_random_numbers(step)

         random_index = random_index + 1
         irep=int(RANDOM_ARRAY(random_index)*(Nreplicas-1))
         write(6,*)"RRR:",random_index,irep
         random_index = random_index + 1  ! for deciding the acceptance/rejection 

         if (replicaID .ne. irep .and. replicaID.ne. (irep+1)) return

         ! remove old files
         if (old_step .gt. 0 .and. my_pe .eq. 0) then
          write(comm_file,'(a7,i4.4,a1,i8.8,a2)') 'rm for_',replicaID,
     &                                    '_',old_step,'.*'
          call system(comm_file(1:22))
         end if


        if (replicaID.eq.irep)  then
          if (prll_on_off) call reduce_energies()
          e_rep(1)=e_total
          t_rep(1)=tempi(1)
          tmp1 = -1
          tmp2 = -1

          call getunit(ucomm)

          if (my_pe .eq. 0) then

          ! write my energy and temperature to a file
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID+1,'_'
     &                                    ,step,'.dat' 
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          write(ucomm,103,err=7)e_rep(1),t_rep(1)
          close(ucomm)
        
          ! check that the data are correctly written
          do while (abs(tmp1 - e_rep(1)) .gt .1.0  .or. 
     &              abs(tmp2 - t_rep(1)) .gt. 1.0) 
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            read(ucomm,103,err=27,end=27)tmp1,tmp2
27          continue            
            close(ucomm)
          end do
          
103       format(f20.8,1x,f10.5)

          ! let the processors in the neighboring replica know that
          ! the data is successfully written
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID+1,'_'
     &                                    ,step,'.don'
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          write(ucomm,103,err=7)0.0,0.0
          close(ucomm)

          end if

          ! wait until the neigh. replica writes data for me to disk
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,'_'
     &                                    ,step,'.don'
        
          exist = .false.
          write(6,*)"waiting for",replicaID+1
          do while (.not. exist)
            inquire(file=comm_file(1:flen),exist=exist)
          end do

          ! read the neighbor energy and temperature from the disk
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,'_'
     &                                    ,step,'.dat'
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          read(ucomm,103,err=8,end=8)e_rep(2),t_rep(2)
          close(ucomm)

          call freeunit(ucomm)
          
        endif 

        
        if(replicaID.eq.(irep+1)) then 
          if (prll_on_off) call reduce_energies()
          e_rep(2)=e_total
          t_rep(2)=tempi(1)
          tmp1 = -1
          tmp2 = -1
      
          call getunit(ucomm)

          if (my_pe .eq. 0) then

          ! write my energy and temperature to a file
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID-1,'_'
     &                                    ,step,'.dat'
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          write(ucomm,103,err=7)e_rep(2),t_rep(2)
          close(ucomm)

          ! check that the data are correctly written
          do while (abs(tmp1 - e_rep(2)) .gt .1.0  .or.
     &              abs(tmp2 - t_rep(2)) .gt. 1.0)
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            read(ucomm,103,err=17,end=17)tmp1,tmp2
17         continue
            close(ucomm)
          end do
          
          ! let the processors in the neighboring replica know that
          ! the data is successfully written
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID-1,'_'
     &                                    ,step,'.don'
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          write(ucomm,103,err=7)0.0,0.0
          close(ucomm)

          end if

          ! wait until the neigh. replica writes data for me to disk
          write(comm_file,'(a4,i4.4,a1,i8.8,a)') 'for_',replicaID,'_'
     &                                    ,step,'.don'

          exist = .false.
          write(6,*)"waiting for",replicaID-1
          do while (.not. exist)
            inquire(file=comm_file(1:flen),exist=exist)
          end do

          ! read the neighbor energy and temperature from the disk
          write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,'_'
     &                                    ,step,'.dat'
          open (unit=ucomm,status='unknown',file=comm_file(1:flen))
          read(ucomm,103,err=8,end=8)e_rep(1),t_rep(1)
          close(ucomm)

          call freeunit(ucomm) 

        endif  

       ! make the decision about the exchange
       call ptempering(e_rep,t_rep,accept)
       
       ! exchange the coordinates and momenta
       if (accept) then
       
         if (replicaID .eq. irep) then

           call getunit(ucomm)

           ! write data
           if (my_pe .eq. 0) then
         
            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID+1,
     &                                    '_',step,'.cor'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            call putdata(ucomm,'CHARM',coor)
            close(ucomm)

            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID+1,
     &                                    '_',step,'.vel'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            call putdata(ucomm,'CHARM',velo)
            close(ucomm)

            ! let the neigbor replica know that data is written
            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID+1,
     &                                    '_',step,'.fin'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            write(ucomm,103,err=7)0.0,0.0
            close(ucomm)
              
           end if

           ! read data
           ! first wait until the data are written by the neighbor to disk
           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.fin'
           exist = .false.
           write(6,*)"waiting for",replicaID+1
           do while (.not. exist)
             inquire(file=comm_file(1:flen),exist=exist)
           end do    
           
           ! read coordinates from the disk
           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.cor'
           
           open (unit=ucomm,status='unknown',file=comm_file(1:flen))
           call getdata(ucomm,'CHARM',coor)
           close(ucomm) 

           ! read velocities from the disk
           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.vel'
           open (unit=ucomm,status='unknown',file=comm_file(1:flen))
           call getdata(ucomm,'CHARM',velo)
           close(ucomm)
           
           ! scale velocities to my temperature
           do ji=1,npt
             do ll=1,3
               velo(ll,ji) = velo(ll,ji) * sqrt(t_rep(1)/t_rep(2))
             end do
           end do

           call freeunit(ucomm)
           
         end if

         if (replicaID .eq. irep+1) then
           
           call getunit(ucomm)

           ! write data
           if (my_pe .eq. 0) then

            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID-1,
     &                                    '_',step,'.cor'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            call putdata(ucomm,'CHARM',coor)
            close(ucomm)

            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID-1,
     &                                    '_',step,'.vel'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            call putdata(ucomm,'CHARM',velo)
            close(ucomm)

            write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID-1,
     &                                    '_',step,'.fin'
            open (unit=ucomm,status='unknown',file=comm_file(1:flen))
            write(ucomm,103,err=7)0.0,0.0
            close(ucomm)
              
           end if

           ! read data
           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.fin'
           exist = .false.
           write(6,*)"waiting for",replicaID-1
           do while (.not. exist)
             inquire(file=comm_file(1:flen),exist=exist)
           end do

           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.cor'
           
           open (unit=ucomm,status='unknown',file=comm_file(1:flen))
           call getdata(ucomm,'CHARM',coor)
           close(ucomm) 

           write(comm_file,'(a4,i4.4,a1,i8.8,a4)') 'for_',replicaID,
     &                                    '_',step,'.vel'
           open (unit=ucomm,status='unknown',file=comm_file(1:flen))
           call getdata(ucomm,'CHARM',velo)
           close(ucomm)

           call freeunit(ucomm)

           ! scale velocities to my temperature
           do ji=1,npt
             do ll=1,3
               velo(ll,ji) = velo(ll,ji) * sqrt(t_rep(2)/t_rep(1))
             end do
           end do 
           
         end if
       
       end if  ! accept

       old_step = step

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
    
      return


7       continue
        call alert(name,namel,'Error during write',18,1)
        return

8       continue
        call alert(name,namel,'Error during read',18,1)
        return

      end
