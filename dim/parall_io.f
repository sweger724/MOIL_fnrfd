      subroutine open_inout(stdi,stdo,index)
      
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      !include 'COMMON/UNITS.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      

      integer index,stdi,stdo
      character*80 filenamei,filename1
      integer len,len1
      
      len=8
      write(filenamei,'(a8)') 'inp_0000'
      open (unit=stdi,status='old',file=filenamei(1:len))
      
      len1=18
      write(filename1,'(a8,i4.4,a4)') 'dyn_out_',my_pe,'.log'
      open (unit=stdo,status='unknown',file=filename1(1:len1))
      
      return
      end
