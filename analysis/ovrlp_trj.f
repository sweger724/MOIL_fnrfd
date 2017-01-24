         program ovrlp_trj

         include 'COMMON/LENGTH.BLOCK'
         include 'COMMON/COORD.BLOCK'
	 include 'COMMON/VELOC.BLOCK'
         include 'COMMON/CONNECT.BLOCK'
         include 'COMMON/LINE.BLOCK'
         include 'COMMON/DEBUG.BLOCK'
         include 'COMMON/FREEZ.BLOCK'
         include 'COMMON/CONVERT.BLOCK'
         include 'COMMON/UNITS.BLOCK'
           
c rms: rms of dynamic structrue respect to reference structure
c rms_ave: rms of dynamic structrue respect to averaged dynamic structure
         integer ucon,urdcrd,urcrd,uwcrd
         integer geti,of,i,j,ncoor,namel,n,k,kk
         integer rbin,wbin
         integer npick,npick2,level
         integer ipick(maxpt),jpick(maxpt)
         integer ipick2(maxpt),jpick2(maxpt)
         double precision dpot(3,maxpt)
c         double precision rms_ave
         double precision rms
         character*11 name
         logical find,fopen,fluct,rmsav,rmsre
         logical pickpt

         logical norew
         integer lpstr1
         COMMON /norw/ norew,lpstr1

         lpstr1=1

         stdi=5
         stdo=6
         rbin = 1
         wbin = 1 

         npt=0
         name='ovrlp_trj'
         namel=9
c  open junk file for rline
c
         jnkf=25
         open(unit=jnkf,status='scratch')
c    defalt parameters
         pickpt=.false.
          fluct=.false.
          rmsav=.false.
          rmsre=.false.
          norew=.true.
          n    =0
          ncoor=1
1         continue
          call rline(name,namel,stdi)
          if (find ('norw')) norew=.true.
          if (find ('file')) then
            if (find ('conn')) then
             ucon=of()
c  get connectivity
             call rconn(ucon)
            end if
            if (find ('rcrd')) then
             urcrd=of()
c   read referance structure
             call getcrd(urcrd,'CHARM')
             do 80 i=1,npt
             coor2(1,i)=coor(1,i)
             coor2(2,i)=coor(2,i)
             coor2(3,i)=coor(3,i)
80           continue
c@ Ron
c		write(*,*) ' coor2 '
c		do 81 i=1,npt
c			write(*,*)(coor2(j,i),j=1,3)
c81		continue
            end if
            if (find ('rdyc')) then
                 urdcrd=of()
                 rewind urdcrd
            end if
            if (find ('wcrd'))  uwcrd=of()
          else
            n = geti('#str',n)
            if (find('pic2')) then
             call pick(ipick2,npick2)
                j = 0
                do 7 i=1,npt
                        if (ipick2(i).ne.0) then
                                j = j + 1
                                jpick2(j) = i
                        end if
7               continue
                npick2 = j      
            end if
            if (find('pick')) then
             call pick(ipick,npick)
                j = 0
                do 6 i=1,npt
                        if (ipick(i).ne.0) then
                                j = j + 1
                                jpick(j) = i
                        end if
6               continue
                npick = j
                write(*,*)' after assignment npick ',npick
             pickpt=.true.
            end if
            if (find('action')) goto 5
          endif
          goto 1
5         continue
          if (.not. fopen(ucon)) then
            level=1
            call alert(name,namel,'ucon not opened',15,level)
          else if (.not. fopen(urdcrd)) then
            level=1
            call alert(name,namel,'urdcrd not opened',16,level)
          end if
c  initialize the vector nofreez and dpot
          inofrz=npt
          do 8 i=1,npt
           nofreez(i)=i
           dpot(1,i) = 0.d0
           dpot(2,i) = 0.d0
           dpot(3,i) = 0.d0
	   velo(1,i) = 0.d0
	   velo(2,i) = 0.d0
	   velo(3,i) = 0.d0
8         continue
c

        norew = .true.
        rewind urdcrd
         do 25 k=1,n
          if (.not.norew) rewind urdcrd
          call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
          call ovrlpck(coor2,coor,dpot,velo,jpick,npick,rms)
	  call wdyncrd(uwcrd,n,k,inofrz,nofreez,wbin)
25        continue
          stop
          end
