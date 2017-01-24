        subroutine prepare_random_numbers(step)

          implicit none

          include 'COMMON/LENGTH.BLOCK'
          include 'COMMON/PARALLEL.BLOCK'
          include 'COMMON/PT.BLOCK'
          include 'mpif.h'

          integer i, namel, step, flen, ucomm
          character *22 name
          character *80 comm_file
          logical exist

          namel = 22
          flen = 17
          

          ! draw sequence of random numbers
          
          ! it is important that all processors in replicaID=0
          ! draw randomSize random numbers, because random number 
          ! generators must be synchronized within a replica
          if (replicaID .eq. 0) call RANLUX(RANDOM_ARRAY,randomSize)
          
          ! broadcast this sequence to all processors/replicas
          if (replicaID .eq. 0 .and. my_pe .eq. 0 
     &         .and. Nreplicas.gt.1) then
            call getunit(ucomm)

            ! remove old files
            write(comm_file,'(a)')"rm -f ran_*.done random_numbers.dat"
            call system(comm_file(1:35))

            
            open (unit=ucomm,status='unknown',form='unformatted',
     &            file='random_numbers.dat')
            write(ucomm,err=8)(RANDOM_ARRAY(i),i=1,randomSize)
            write(6,*)"Generating random array:",RANDOM_ARRAY(1)
            close(ucomm)

            ! wait for the file system to update the file
            call system("sleep 10")

            write(comm_file,'(a4,i8.8,a5)') 'ran_',step,'.done'
            open (unit=ucomm,status='unknown',file=comm_file(1:17))
            write(ucomm,*)0.0,0.0
            close(ucomm)
            
            call freeunit(ucomm)
          end if
          
          if (replicaID .ge. 1) then

6           continue            
          
            call getunit(ucomm)
            write(comm_file,'(a4,i8.8,a5)') 'ran_',step,'.done'
            
             
             exist = .false.
             do while (.not. exist)
               inquire(file=comm_file(1:flen),exist=exist)
             end do
        
             open (unit=ucomm,form='unformatted',status='unknown',
     &             file='random_numbers.dat')
             read(ucomm,err=7)(RANDOM_ARRAY(i),i=1,randomSize)
             write(6,*)"Received random array:",RANDOM_ARRAY(1)
             close(ucomm)
             
             call freeunit(ucomm)
          end if
          
          ! initialize position in the sequence
          random_index = 0
          
          return

8       continue
        call alert(name,namel,'Error during write',18,1)
        return

7       continue
        call alert(name,namel,'Error during read',18,0)
        
        ! close the file 
        close(ucomm)
        call freeunit(ucomm)
        
        ! try again
        goto 6
        
        end
