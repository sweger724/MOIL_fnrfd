      subroutine wpath_all()
c
c.v1.0 (last changed 24/2/97)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c
      integer ucrd
      logical fopen
c
      integer i,j,k,ist,ist1,ist2
c
      ucrd=-1
c
      if(fopen(uwpth)) then
         ucrd=uwpth
         write(stdo,*) ' Rewind and writing unit', ucrd
         rewind ucrd
      endif
c
c.....write only the structures actualy used. 
c.....(bc1) Two points are kept fixed at each extrema but only 
c...........one is used (and printed).
c.....(bc2)  it is not more used !!! (zva)
c
      k=0
      if (first.and.(ucrd.gt.0)) then
         do 100 i=0,pseg
            k = k + npt3
            write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1           (r(k+2+j),j=0,npt3-3,3),
     2           (r(k+3+j),j=0,npt3-3,3)
 100     continue
c     
         if (last) then
            k = k + npt3
            write(ucrd)onsager,(r(k+1+j),j=0,npt3-3,3),
     1           (r(k+2+j),j=0,npt3-3,3),
     2           (r(k+3+j),j=0,npt3-3,3)
            goto 999
         end if
      end if
c
c sending and receiving structures from other processors
c
      ist1=pseg+2
      ist2=pseg*proc_max+2
c
      do 200 ist=ist1,ist2
         call pgstruc(ist) 
         if(first.and.(ucrd.gt.0)) then
            write(ucrd) onsager,((coor(j,i),i=1,npt),j=1,3)
         endif
 200  continue
c     
 999  continue
c
      if(first.and.(ucrd.gt.0)) then
         write(stdo,*)'closing unit', ucrd
         call close_open_bin(ucrd)
      endif
c     
      return
      end
      
c
c  This part should go to the communications !!!!!
c
      subroutine pgstruc(ist)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
c
      include 'mpif.h'
c
c
      integer i,n,ist,ist1,ist2,l,ll
      integer rc
      integer status(MPI_STATUS_SIZE)
      character*80 err_msg
c
      ist1=pseg+2
      ist2=pseg*proc_max+2
      rc=0
c
      if(ist.lt.ist1.or.ist.gt.ist2) then
        write(stdo,*)
     >      ' Error: wrong structure',ist
        stop 100
      endif
c
      n=(ist-2)/pseg
      if(ist.eq.ist2) n=n-1
      i=ist-(n*pseg)+1   
c 
c      write(stdo,*) "   PATH: STRUCTURE =",ist
c      write(stdo,'(4x,2(a,i4))') "struc=",i," from Proc=",n
c      return
c       
      if((n.lt.1).or.(n.ge.proc_max)) then
        write(stdo,*)
     >      ' Error: out of target processor ! n=',n
        stop 101
      endif   
c       
      if((i.lt.3).or.(i.gt.(pseg+3))) then
        write(stdo,*)
     >      ' Error: out of target loc. structure ! =',i
        stop 102
      endif
c  
      if(first) then 
         write(stdo,'(4x,2(a,i4))') "Receive struc=",i," from Proc=",n
         Call MPI_Recv(coor, npt3, MPI_DOUBLE_PRECISION,
     >        n,ist,MPI_COMM_WORLD, status,rc )
         
         if (rc.ne.0) then
            write(stdo,*)' rc = ',rc
            write(stdo,*)   
     >           ' Error on receiv. struc i=',ist," from proc ID=",n
            call error message(rc,err_msg)
            write(stdo,*)err_msg(1:80)
            stop 103
         end if
      elseif (my_pe.eq.n) then
c     
         write(stdo,'(4x,2(a,i4))') "Sent struc=",i," from Proc=",n
         ll = -3 + npt3*(i-1)
         do 1000 l=1,npt
            ll = ll + 3
            coor(1,l) = r(ll+1)
            coor(2,l) = r(ll+2)
            coor(3,l) = r(ll+3)
 1000    continue  
c
         Call MPI_Send(coor, npt3, MPI_DOUBLE_PRECISION,
     >        0,ist,MPI_COMM_WORLD, rc )

         if (rc.ne.0) then
            write(stdo,*)' rc = ',rc
            write(stdo,*)
     >           ' Error on sending struc i=',ist," from proc ID=",n
            call error message(rc,err_msg)
            write(stdo,*)err_msg(1:80)
            stop 104
         end if
      endif
c
      return
      end
c



