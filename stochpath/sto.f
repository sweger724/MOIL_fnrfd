      program STO
c         
      include 'COMMON/LENGTH.BLOCK'  
      include 'COMMON/DYNA.BLOCK'   
      include 'COMMON/ACTPARA.BLOCK'  
      include 'COMMON/UNITS.BLOCK'
c     
c     
      integer ierr
      logical find,fopen
c     
c..   new local variables.
c.....name   - name of the subroutine that calls rline
c.....namel  - lenght of name
c.....dyna   - flag to direct the program to either simple minimization 
c     
      character*6 name
      integer namel
c     
      integer i
c     
c.....GENERAL INITIALIZATION
c     
      call init_sto
c  
c     
c.....And now we are ready to call the sto1 routine
c     
c      if (dyna) then
         write(stdo,'(/10x,a/)') 
     >              'Simulating annealing using dynamics !!!'
         call sto1_dyn
c      end if
c     
c     
      write(stdo,'(/10x,a)') "Summary CPU for STO :"
      call printcpu(stdo)
      write(stdo,'(/1x,5a10)') ("----------",i=1,5)
c     
ccccc      call printrtc(stdo)
c     
      write(stdo,'(/1x,5a10)') ("**********",i=1,5)
c
      if (fopen(stdo)) close(stdo)
      if (fopen(stdi)) close(stdi)
      if (fopen(uwcrd)) close(uwcrd)
      if (fopen(uwpth)) close(uwpth)
c     
      call comm_exit(ierr)
c 
c     
      stop
      end
