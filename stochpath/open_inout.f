      subroutine open_inout(stdi,stdo,wpth,index)
c     ======================================
c
c.v1.0 (last changed 06/23/99) Zaloj
c
c opens the input and output files of the parallel program
c
c declaration part:
c =================
c
      integer index,stdi,stdo,wpth
      character*80 filenamei,filename1,filename2
      integer len,len1,len2
c
c
c execution part:
c ===============
c
c 
      len=10
c      write(filenamei,'(a4,i4.4)') 'inp_',index
      write(filenamei,'(a10)') 'allsto.inp'
c     
      len1=18
c
      write(filename1,'(a8,i4.4,a4)') 'sto_out_',index,'.log'
c     
      len2=18
c
      write(filename2,'(a8,i4.4,a4)') 'out_bin_',index,'.pth'
c     
      open (unit=stdi,status='old',file=filenamei(1:len))
c
      open (unit=stdo,status='unknown',file=filename1(1:len1))
c
      open (unit=wpth,status='unknown',form='unformatted',
     $        file=filename2(1:len2))
c
      return
      end

