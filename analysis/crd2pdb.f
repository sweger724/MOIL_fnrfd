       program CRD2PDB
c
       integer stdi,stdo,stderr
       character*80 line
       character*4 lbres,lbat
       character*3 lkeepH
       integer mnumat,natm,ires,imol,ichainInc
       character*1 lmol
       character*1 mollab(0:26) /' ','A','B','C','D','E','F','G',
     & 'H','I','J','K','L','M','N','O','P','Q','R','S',
     & 'T','U','V','W','X','Y','Z'/
       character*4 HydroArg(16) /'1HH1','1HH1','2HH1','1HH2','1HH2',
     &  '2HH1','2HH2','2HH2','HH11','1HH1','HH12','1HH2','HH21',
     &  '2HH1','HH22','2HH2'/
c
c......UNITS
       stdi=5
       stdo=6
       stderr=0
c......write(*,'(a40)')'File with more than one monomer (Y/N)?'
       imol=1
           ichainInc=0
       if ( iargc() .ge. 1 ) then
         call getarg(1, lkeepH)
         if (lkeepH(1:1).eq.'N' .or. lkeepH(1:1).eq.'n'
     &      .or. lkeepH(1:1).eq. '1' .or. lkeepH(1:1).eq. '0') then
           imol=0
         endif
       endif
       ltakres= 0
       lmol=mollab(imol)
       write(stdo,'(A33)')   'REMARK  99 File generared by MOIL'
       write(stdo,'(A33)')   'REMARK  99                       '
  50   read(stdi,'(a)', END= 500) line
       if ( line(1:1) .eq. '*' ) goto 50
       if ( line(10:10).ne.' ' ) goto 101
  100  read(stdi,'(a)', END= 500) line
  101  read(line(1:54),'(2i7,2(1X,A4),3F10.5)') natm,ires,lbres,lbat,
     &                x,y,z
        
       if (((lbres(1:4).eq.'NTER').or.(lbres(1:3).eq.'NTR'))
     &   .and.(ichainInc.eq.1)) then
          imol=imol+1
          lmol=mollab(imol)
          ltakres=ires-1
                  ichainInc=0
           endif

       if ((lbres(1:3).eq.'HEM').or.(lbres(1:3).eq.'TIP')) then
          write(stdo,120)natm,lbat,lbres,lmol, ires-ltakres,x,y,z
!HETATM   32 1HH1 HEM A   2      69.311  52.877  21.177  
  120     FORMAT ('HETATM',I5,2X,2A4,A1,I4,4X,3F8.3)
       elseif ((lbat(4:4).eq.'1').or.(lbat(4:4).eq.'2')) then
          if (lbres(1:3).eq.'ARG') then
             do i=0,8
                if (lbat(1:4).eq.HydroArg(2*i+1)) then
                   lbat(1:4)=HydroArg(2*i+2)
!                  write(*,*)' lbat HydA ',lbat(1:4),' ',HydroArg(2*i+1)
                   goto 130
                endif
             enddo
!               write(*,*)' lbat',lbat(1:4)
!             write(stderr,'(A25,/,A80)')
!     &            'Hydrogen not recognized',line(1:80)
             stop 6666
          endif
  130     continue
          write(stdo,364)natm,lbat,lbres,lmol, ires-ltakres,x,y,z
!ATOM     32 1HH1 ARG A   2      69.311  52.877  21.177  
  364     FORMAT ('ATOM',I7,1X,A4,1X,A4,A1,I4,4X,3F8.3)
       else
          write(stdo,365)natm,lbat,lbres,lmol, ires-ltakres,x,y,z
  365     FORMAT ('ATOM',I7,2X,2A4,A1,I4,4X,3F8.3)
       endif
!ATOM      2  N   LYS     2       0.000   0.000   0.000 -4.28 -0.33
!ATOM      1  N   SER A   1      76.336  51.991  27.698
          if ((lbres(1:4).eq.'CTER').or.lbres(1:3).eq.'CTR') then
           write(stdo,'(A3)')'TER'
                   ichainInc=1
       endif
       goto 100
  500  continue
       stop
       end



















