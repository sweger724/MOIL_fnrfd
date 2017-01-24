      subroutine open_inout(stdi,stdo,index)
c     ======================================
c
c     (last changed 04/24/06) Peter Majek
c
c opens the input and output files of the parallel program
c
c declaration part:
c =================
c
      integer index,stdi,stdo
      character*80 filenamei,filename1
      integer len,len1
c
c
c execution part:
c ===============
c
c 
      len=8
      write(filenamei,'(a4,i4.4)') 'inp_',index
	write(*,*)' filenamei ', filenamei(1:8)
c     
      len1=18
c
      write(filename1,'(a8,i4.4,a4)') 'pth_out_',index,'.log'
c     
      open (unit=stdi,status='unknown',file=filenamei(1:len))
      rewind stdi
c
      open (unit=stdo,status='unknown',file=filename1(1:len1))
c
      return
      end

