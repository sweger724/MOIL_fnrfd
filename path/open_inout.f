      subroutine open_inout(stdi,stdo,index)
c
c opens the input and output files of the parallel program
c
c declaration part:
c =================
c
        implicit none
      integer index,stdi,stdo
      character*80 filenamei,filename1,filename2
      integer len,len1,len2
c
c execution part:
c ===============
c
      len=8
      write(filenamei,'(a4,i4.4)') 'inp_',index
        write(6,*)' filenamei ', filenamei(1:len)
     
      len1=18

      write(filename1,'(a8,i4.4,a4)') 'pth_out_',index,'.log'
     
      len2=18

      write(filename2,'(a8,i4.4,a4)') 'out_bin_',index,'.pth'
      
     
      open (unit=stdi,status='unknown',file=filenamei(1:len))
      rewind stdi

      open (unit=stdo,status='unknown',file=filename1(1:len1))

      return
      end

