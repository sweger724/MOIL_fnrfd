      program rms_2dcd
c
c compare two dcd files of different length. All against all rmsd is computed
c the same connectivity file is assumed.
c
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/CCRD.BLOCK'
c
c
c
      logical dmin
      integer ipick(maxpt) ,npick 
      integer of,geti,nstru,logic
      integer namel,level
      integer ucon,uwrms,urcrd1,urcrd2,uallrms
      integer i,j,i1,i2,jmin
      integer nst1,nst2
      double precision rmsall,rms_min
      double precision mmm(maxpt)
      character*9 name
      logical find,fopen
      logical pickpt
      data ucon,uallrms,uwrms,urcrd1,urcrd2/5*99/

      norew = .false.
      dmin = .false.
      lpstr = 1

      stdi=5
      stdo=6
      totmon=0
      npt=0
      name='rms_2dcd'
      namel=8
      logic=1
c     open junk file for rline
c     
      jnkf=25
      open(unit=jnkf,status='scratch')
c
c     defalt parameters
c
c
      nst1=1
      nst2=1
c
c    line loop over the input data
c
      pickpt = .false.
 1    continue
      call rline(name,namel,stdi)
      if (find('file')) then
         if (find ('conn')) then 
c     open connectivity file
            ucon=of()
c     read connectivity
            call rconn(ucon)
c     open coordinate file #1
         else if(find ('rdc1')) then
            urcrd1=of()
         else if(find ('rdc2')) then
            urcrd2=of()
c     open output file for rms data
         else if(find('wrms')) then  
            uwrms=of()
         else if(find('wall')) then  
            uallrms=of()
         endif
      else 
         nst1=geti('st1s',nst1)
         nst2=geti('st2s',nst2)
         if (find('sele')) then
            call pick(ipick,npick)
            pickpt=.true.
         end if
         if (find ('dmin')) dmin = .true.
         if (find ('action')) goto  5
      end if
      go to 1
 5    continue
c
c---------------------------------------------
c
      if (.not.fopen(ucon)) then
         level=1
         call alert(name,namel,'ucon not opened',15,level)
      else if (.not. fopen(urcrd1)) then
         level=1
         call alert(name,namel,'urcrd1 not opened',17,level)
      else if (.not.fopen(urcrd2)) then
         level=1
         call alert(name,namel,'urcrd2 not opened',17,level)
      else if (.not.fopen(uwrms)) then
         level=1
         call alert(name,namel,'uwrms not opened',16,level)
      else if (.not.fopen(uallrms)) then
         level=1
         call alert(name,namel,'uallrms not opened',18,level)
      end if
c
c     use pick to store masses to zero
c      
      if(pickpt) then
         npick = 0
         do 100 i=1,npt
            if (ipick(i).eq.0) then
             mmm(i) = 0.d0
            else
             npick = npick + 1
             mmm(i) = ptms(i)
            end if
 100    continue
      else
        npick = npt
        do i=1,npick
         mmm(i) = 1.d0
        end do
      end if
c
c     initialize the vector nofreez  
c    
      inofrz=npt
      do 200 i=1,npt
         nofreez(i)=i
 200  continue
c      
c
c     
      write(*,*)'In Milestoning analysis, reporting first the '
      write(*,*)'index of dynamics then nearest anchor then rmsd'
      do i =1,nst1
       rewind urcrd1
       call rdyncrd(urcrd1,i,inofrz,nofreez,1)
       do i1 =1,npt
        do i2 = 1,3
         velo(i2,i1)=coor(i2,i1)
        end do
       end do
       rms_min = 9999.d0
       jmin = 0
       do j = 1,nst2
        rewind urcrd2
        call rdyncrd(urcrd2,j,inofrz,nofreez,1)
        call rmsd_weight(npt,velo,coor,rmsall,.false.,mmm)
        if (dmin) then
         if (rmsall.lt.rms_min) then
          rms_min = rmsall
          jmin = j
         end if
        else
         write(*,*)' i j rmsd ' ,i,j,rmsall
         write(uallrms,*)i,j,rmsall
        end if
       end do
       if (dmin) write(uwrms,*) i,jmin,rms_min
      end do
      
      stop
      end
