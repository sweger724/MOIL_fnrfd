c -------------------------------------------------------------------------
c Gather the full list of shaked bonds in each processor

        subroutine gather_shake(nshak,ishak1,ishak2,dissq,shak_diff)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'mpif.h'
c        include '../include/Tcommlib_f.h'
c        include '../include/Trtlib_f.h'

        integer ishak1(*),ishak2(*),nshak
        double precision dissq(*),shak_diff(*)
        integer i,one,two,error1,error2
        integer nshaks(0:maxpe),displace(0:maxpe)        
        integer nshaks_tmp(0:maxpe)
        character*80 message

c@ RE   call tcommlib all gather(nshak,nshaks,1,4,error1)
c@ RE : Nov 21,2005 - adjusting MPI_ALLgather to new syntax
c
        one = 1
        two = 2
        write(*,*)' before MPI call #1 '
        write(*,*)' nshak, nshaks, MPI_INTEGER, MPI_COMM_WORLD, error1 '
        write(*,*) nshak, nshaks, MPI_INTEGER, MPI_COMM_WORLD, ierror 
        call MPI_allgather(nshak,one,MPI_INTEGER,nshaks(0),one,
     1          MPI_INTEGER,MPI_COMM_WORLD,ierror)
c@
        write(*,*)' nshaks ',(nshaks(i),i=0,num_pes-1)
        if (ierror.gt.0) then
                write(stdo,*)' Error in gather_shake '
                write(stdo,*)' error = ',ierror
                write(stdo,1)message
1               format(1x,a80)
        end if 
 
        nshak = 0
        displace(0)=0
        do 10 i=1,num_pes
           displace(i) = displace(i-1)+nshaks(i-1)
           nshak = nshak + nshaks(i-1)
10      continue

c@ RE        call tcommlibAllGatherV(ishak1,ishak1,4,displace,nshaks,error1)
c@ RE        call tcommlibAllGatherV(ishak2,ishak2,4,displace,nshaks,error2)
                write(*,*)' Just before all gatherV '
                write(*,*)' my_pe nshaks(my_pe) ishak1(1:nshaks(my_pe) '
            write(*,*) my_pe,nshaks(my_pe),(ishak1(i),i=1,nshaks(my_pe))
        call MPI_allgatherV(ishak1,nshaks(my_pe),MPI_INTEGER
     1   ,ishak1,nshaks_tmp,displace,MPI_INTEGER,MPI_COMM_WORLD,ierror)
        write(*,*)' after allgatherV #1 '
        if (ierror.ne.0) then
                write(*,*)' ierror = ',ierror
                stop
        end if
        call MPI_allgatherV(ishak2,nshaks(my_pe),MPI_INTEGER
     1   ,ishak2,nshaks_tmp,displace,MPI_INTEGER,MPI_COMM_WORLD,ierror)
c@
c@      write(*,*)' after allgatherv #1 ishak1 ',(ishak1(i),i=1,nshak)
c@      write(*,*)' after allgatherv #1 ishak2 ',(ishak2(i),i=1,nshak)
c@      stop
        if ((error1.gt.0).or.(error2.gt.0)) then
                if (error2.gt.0) error1 = error2
                write(stdo,*)' Error in gather_nshaks'
                write(stdo,*)' error = ',error1
                write(stdo,1)message
         end if

c@ RE        call tcommlibAllGatherV(dissq,dissq,8,displace,nshaks,error1)
c@ RE        call tcommlibAllGatherV(shak_diff,shak_diff,8,displace,
c@ RE     1                              nshaks,error2)
        call MPI_allgatherV(dissq,nshaks(my_pe),MPI_DOUBLE_PRECISION
     1   ,dissq,nshaks,displace,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)
        call MPI_allgatherV(shak_diff,nshaks(my_pe),MPI_DOUBLE_PRECISION
     1   ,shak_diff,nshaks,displace,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)

        if ((error1.gt.0).or.(error2.gt.0)) then
                if (error2.gt.0) error1 = error2
C@                call tcomm lib error message(error1,message)
                write(stdo,*)' Error1 in gether_nshaks'
                write(stdo,*)' error = ',error1
                write(stdo,1)message
         end if

           return
           end

c ************************************************************
        subroutine gather_shake_balance (shake_bal)

        integer shake_bal(0:*)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        integer error
        character*80 message

        call MPI_allgather(shake_bal(0),1,MPI_INTEGER,shake_bal,1,
     1          MPI_INTEGER,MPI_COMM_WORLD,error)

        if (error.gt.0) then
                write(stdo,*)' Error in gather_shake_balance'
                write(stdo,*)' error = ',error
                write(stdo,1)message
        end if
1       format(80a)
        

        return
        end

c *******************************************************************   
c get the shared particles from each processor
        subroutine get_shared_pt(shak_l,count,pts_num,
     *                    RBpts,p_dissq,p_shak_diff)


        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer shak_l(2,*)
        double precision RBpts(3,*)
        double precision p_dissq(*),p_shak_diff(*)
        double precision x,y,z,offset
        integer displace(0:maxpe),error,pts_num,i,j,temp,count
        integer npts2(0:maxpe),idx,pl1,pl2
        character*80 message

        pts_num = 0
        displace(0) = 0
        do 20 i=1,num_pes
           displace(i) = displace(i-1) + nptsh(i-1)
           pts_num = pts_num + nptsh(i-1)
           npts2(i-1) = 0
20      continue

c get the constrains for the rigid body itself (with no bonds
c to the light particles)
        temp = my_pe*4 + displace(my_pe)
        count = 0
        do 15 i=1,4
         do 25 j=i+1,4+nptsh(my_pe)
            x = RBpts(1,i) - RBpts(1,j)
            y = RBpts(2,i) - RBpts(2,j)
            z = RBpts(3,i) - RBpts(3,j)
            count = count + 1
            shak_l(1,count) = i+temp
            shak_l(2,count) = j+temp
            p_dissq(count) = (x*x + y*y + z*z)
            p_shak_diff(count) = p_dissq(count)*epshak
 25      continue
 15     continue

        count = 0
        do 75 i=1,num_pes
           temp = 4 + nptsh(i-1)
           npts2(i-1) = 6 + nptsh(i-1)*4
           displace(i) = displace(i-1) + npts2(i-1)
           count = count + npts2(i-1)     
 75     continue

c@        call tcommlibAllGatherV(p_dissq,p_dissq,8,displace,npts2,error)
c@        call tcommlibAllGatherV(p_shak_diff,p_shak_diff,8,
c@     1                                             displace,npts2,error)
c@      call tcommlibAllGatherV(shak_l,shak_l,2*4,displace,npts2,error)

        if (error.gt.0) then
c@                call tcomm lib error message(error,message)
                write(stdo,*)' Error in get_shared_pt '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if
        
c       do 99 i=1,count    
c          write(*,*)i,shak_l(1,i),shak_l(2,i),p_dissq(i)
c 99       continue

               
        return
        end


        subroutine gatherRB(pts,velocs,mass,my_pts,my_vel,
     *                                   pts_num,nptsh,displace)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        integer pts_num,nshpt(0:maxpe),nptsh(0:maxpe)
        integer displace(0:maxpe),error,dis,error1,i,j
        double precision pts(3,*),velocs(3,*),mass(*)
        double precision my_pts(3,*),my_vel(3,*)
        character*80 message

c get the rigid body points , and velocities from each processor        
c each processor has 4 points to represent the rigid body
c and additional points for the shared particles


        displace(0) = 0
        pts_num = 0
        do 10 i=1,num_pes
           nshpt(i-1) = nptsh(i-1)+4
           displace(i) = displace(i-1) + nshpt(i-1)
           pts_num = pts_num + nshpt(i-1)
10      continue

c@RE    call tcommliballGatherV(my_pts,pts,8*3,displace,nshpt,error)
c       call tcommliballGatherV(my_vel,velocs,8*3,displace,nshpt,error1)
c@RE    call tcommliballGatherV(mass,mass,8,displace,nshpt,error1)
        if ((error1.gt.0).or.(error.gt.0)) then
                if (error.gt.0) error1 = error
c@RE                call tcomm lib error message(error1,message)
                write(stdo,*)' Error in getRB'
                write(stdo,*)' error = ',error1
                write(stdo,1)message
1               format(1x,a80)
        end if

        do 50 i=1,pts_num
           velocs(1,i) = 0.d0
           velocs(2,i) = 0.d0
           velocs(3,i) = 0.d0
 50     continue
        
        return
        end


c@      subroutine gather_int (my_val,list)
c@
c@      include 'COMMON/LENGTH.BLOCK'
c@      include 'COMMON/PARALLEL.BLOCK'
c@      include 'COMMON/UNITS.BLOCK'
c@
c@      integer my_val,list(*),error
c@      character*80 message
c@
c@      call tcommlibAllGather(my_val,list,1,4,error)
c@        if (error.gt.0) then
c@               call tcomm lib error message(error,message)
c@               write(stdo,*)' Error in gather int '
c@               write(stdo,*)' error = ',error
c@               write(stdo,1)message
c@1              format(1x,a80)
c@        end if
c@
c@      return 
c@      end
C@      
C@      subroutine gather_int_pairs(my_array,array,displace,nptsh)
C@
C@      include 'COMMON/LENGTH.BLOCK'
C@      include 'COMMON/PARALLEL.BLOCK'
C@      include 'COMMON/UNITS.BLOCK'
C@
C@      integer array(2,*),my_array(2,*),i
C@      integer nptsh(0:maxpe),displace(0:maxpe),error
C@      character*80 message
C@
C@c     do 10 i=1,nptsh(my_pe)
C@c        write(*,*)'I have',my_array(1,i),my_array(2,i)
C@c 10     continue
C@      call tcommlibAllGatherV(my_array,array,4*2,displace,nptsh,error)
C@        if (error.gt.0) then
C@                call tcomm lib error message(error,message)
C@                write(stdo,*)' Error in gather double pairs '
C@                write(stdo,*)' error = ',error
C@                write(stdo,1)message
C@1               format(1x,a80)
C@        end if
C@      
C@      return
C@      end

C@      subroutine gather_double(array,displace,nptsh)
C@
C@      include 'COMMON/LENGTH.BLOCK'
C@      include 'COMMON/PARALLEL.BLOCK'
C@      include 'COMMON/UNITS.BLOCK'
C@
C@
C@      double precision array(2,*)
C@      integer nptsh(0:maxpe),displace(0:maxpe),error
C@      character*80 message
C@
C@      call tcommlibAllGatherV(array,array,8,displace,nptsh,error)
C@        if (error.gt.0) then
C@                call tcomm lib error message(error,message)
C@                write(stdo,*)' Error in gather double pairs '
C@                write(stdo,*)' error = ',error
C@                write(stdo,1)message
C@1               format(1x,a80)
C@        end if
C@      
C@      return
C@      end

        subroutine gather_double1(array,displace,nptsh)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/UNITS.BLOCK'


        double precision array(*)
        integer nptsh(0:maxpe),displace(0:maxpe),error
        character*80 message

C@      call tcommlibAllGatherV(array,array,8,displace,nptsh,error)
        call MPI_allgatherV(array,nptsh,MPI_DOUBLE_PRECISION
     1   ,array,nptsh,displace,MPI_DOUBLE_PREFCISION,MPI_COMM_WORLD)
        if (error.gt.0) then
c@                call tcomm lib error message(error,message)
                write(stdo,*)' Error in gather double pairs '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if

        return
        end

C@      subroutine gather_velo(array,displace,nptsh)
C@
C@      include 'COMMON/LENGTH.BLOCK'
C@      include 'COMMON/PARALLEL.BLOCK'
C@      include 'COMMON/UNITS.BLOCK'
C@
C@
C@      double precision array(3,*)
C@      integer nptsh(0:maxpe),displace(0:maxpe),error
C@      character*80 message
C@
C@c@    call tcommlibAllGatherV(array,array,8,displace,nptsh,error)
C@        if (error.gt.0) then
C@c@                call tcommlib error message(error,message)
C@                write(stdo,*)' Error in gather double pairs '
C@                write(stdo,*)' error = ',error
C@                write(stdo,1)message
C@1               format(1x,a80)
C@        end if
C@      return
C@      end
  
        subroutine gather_int1(array,displace,nptsh)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/UNITS.BLOCK'


        integer array(*),i
        integer nptsh(0:maxpe),displace(0:maxpe),error
        character*80 message

c@      call tcommlibAllGatherV(array,array,4,displace,nptsh,error)
        call MPI_AllGatherV(array,nptsh,MPI_INETEGR,array,nptsh,displace,
     1          MPI_INTEGER,MPI_COMM_WORLD,error)
        if (error.gt.0) then
c@               call tcomm lib error message(error,message)
                write(stdo,*)' Error in gather double pairs '
                write(stdo,*)' error = ',error
                write(stdo,1)message
1               format(1x,a80)
        end if

        return
        end


c Get integer array from the right neighbor, and send to the left.
        
        subroutine send_right_int (send,recv,size)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        integer send(*),recv(*),size,error
        character*80 message

c@      call TcommlibShiftRight(send,recv,size*4,error)

        if (my_pe.ne.num_pes-1) then
c this processor is not the last
         call MPI_send(send,size,MPI_INTEGER,my_pe+1,0,
     1          MPI_COMM_WORLD,error)
        else
c this is the last processor send back to the beginning

         call MPI_send(send,size,MPI_INTEGER,0,0,
     1          MPI_COMM_WORLD,error)
        end if

        if (my_pe.ne.0) then
         call MPI_recv(send,size,MPI_INTEGER,my_pe-1,0,
     1          MPI_COMM_WORLD,error)
        else
         call MPI_recv(send,size,MPI_INETGER,num_pes-1,0,
     1          MPI_COMM_WORLD,error)
        end if
        
        if (error.gt.0) then
c@                call tcommlib error message(error,message)
                write(*,*)' Error in get_right_int '
                write(*,*)' error = ',error
                write(*,1)message
1               format(1x,a80)
        end if

        return 
        end

        subroutine send_right_double (send,recv,size)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        double precision send(*),recv(*)
        integer size,error
        character*80 message

c@      call TcommlibShiftRight(send,recv,size*8,error)
c this processor is not the last

C SEND RIGHT
        if (my_pe.ne.num_pes-1) then
         call MPI_send(send,size,MPI_DOUBLE_PRECISION,my_pe+1,0,
     1          MPI_COMM_WORLD,error)
        else
c this is the last processor send back to the beginning

         call MPI_send(send,size,MPI_DOUBLE_PRECISION,0,0,
     1          MPI_COMM_WORLD,error)
        end if

C RECV RIGHT
        if (my_pe.ne.0) then
         call MPI_recv(send,size,MPI_DOUBLE_PRECISION,my_pe-1,0,
     1          MPI_COMM_WORLD,error)
        else
         call MPI_recv(send,size,MPI_DOUBLE_PRECISION,num_pes-1,0,
     1          MPI_COMM_WORLD,error)
        end if
        
        if (error.gt.0) then
c@                call tcomm lib error message(error,message)
                write(*,*)' Error in send_right_double '
                write(*,*)' error = ',error
                write(*,1)message
1               format(1x,a80)
        end if

        return 
        end

        subroutine send_left_int (send,recv,size)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        integer send(*),recv(*),size,error
        character*80 message

c@      call TcommlibShiftLeft(send,recv,size*4,error)


        if (my_pe.ne.0) then
c this processor is not the first
         call MPI_send(send,size,MPI_INTEGER,my_pe-1,0,
     1          MPI_COMM_WORLD,error)
        else
c this is the first processor send back to the end

         call MPI_send(send,size,MPI_INTEGER,num_pes-1,0,
     1          MPI_COMM_WORLD,error)
        end if

        if (my_pe.ne.num_pes-1) then
         call MPI_recv(send,size,MPI_INTEGER,my_pe+1,0,
     1          MPI_COMM_WORLD,error)
        else
         call MPI_recv(send,size,MPI_INETGER,0,0,
     1          MPI_COMM_WORLD,error)
        end if
        
        if (error.gt.0) then
c@                call tcomm lib error message(error,message)
                write(*,*)' Error in get_left_int '
                write(*,*)' error = ',error
                write(*,1)message
1               format(1x,a80)
        end if

        return 
        end

        subroutine send_left_double (send,recv,size)

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        double precision send(*),recv(*)
        integer size,error
        character*80 message

c@      call TcommlibShiftLeft(send,recv,size*8,error)

C SEND LEFT
        if (my_pe.ne.0) then
         call MPI_send(send,size,MPI_DOUBLE_PRECISION,my_pe-1,0,
     1          MPI_COMM_WORLD,error)
        else
c this is the last processor send back to the beginning

         call MPI_send(send,size,MPI_DOUBLE_PRECISION,num_pes-1,0,
     1          MPI_COMM_WORLD,error)
        end if

C RECV LEFT
        if (my_pe.ne.num_pes-1) then
         call MPI_recv(send,size,MPI_DOUBLE_PRECISION,my_pe+1,0,
     1          MPI_COMM_WORLD,error)
        else
         call MPI_recv(send,size,MPI_DOUBLE_PRECISION,0,0,
     1          MPI_COMM_WORLD,error)
        end if
        
        if (error.gt.0) then
c@                call tcomm lib error message(error,message)
                write(*,*)' Error in send_left_double '
                write(*,*)' error = ',error
                write(*,1)message
1               format(1x,a80)
        end if

        return 
        end


        subroutine gather_coord (crd,npt,pt_start,pt_end,disp,nptsh)
        
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        double precision crd(3,*)
        integer npt,pt_start,pt_end,disp(0:maxpe),nptsh(0:maxpe)

        integer i,k,error
        character*80 message

        k = pt_start - 1
        do 45 i=pt_start,pt_end
           crd(1,i-k)  = crd (1,i) 
           crd(2,i-k)  = crd (2,i) 
           crd(3,i-k)  = crd (3,i)
 45     continue

c       call tcommlibAllGatherV(crd,crd,8,disp,nptsh,error)
        if (error.gt.0) then
c@         call tcommlib error message(error,message)
           write(*,*)' Error in gather double pairs '
           write(*,*)' error = ',error
           write(*,1)message
 1         format(1x,a80)
        end if

        return
        end




