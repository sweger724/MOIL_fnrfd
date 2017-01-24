      program rms_p2p
c
c     --------------------------------------------------------------------------:
c
c     Compute rms of all structures included in one path file against all
c     structures in other path  file. On output generate a matrix of global rms
c     of all (i,j) pairs.
c
c     --------------------------------------------------------------------------:
c
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
      include 'COMMON/OVERLAP.BLOCK'
c
c
c
        parameter(rms_matrix_length=10000)
      integer ipick(maxpt) ,npick 
      integer of,geti,nstru,logic
      integer namel,level
      integer ucon,uwrms,upth1,upth2
      integer i,j,k,l,index,length1,length2
c      double precision rmsall(rms_matrix_length)
c @TFB: gfortran complains about above; instead use integer constant ->
c      double precision rmsall(10000)
      double precision rms,mmm(maxpt)
      character*7 name
      logical find,fopen
      logical pickpt
      data ucon,uwrms,upth1,upth2/4*99/

      stdi=5
      stdo=6
      name='rms_p2p'
      namel=7
        npick = 0
        length1 = 0
        length2 = 0
        pickpt = .false.
c     open junk file for rline
c     
      jnkf=25
      open(unit=jnkf,status='scratch')
c
c     defalt parameters
c
c    line loop over the input data
c
 1    continue
      call rline(name,namel,stdi)
      if (find('file')) then
         if (find ('conn')) then 
c     open connectivity file
            ucon=of()
c     read connectivity
            call rconn(ucon)
         else if (find ('rpt1')) then
            upth1=of()
            rewind upth1
         else if(find ('rpt2')) then
            upth2=of()
            rewind upth2
         endif
      else 
         if (find('sele')) then
            call pick(ipick,npick)
            pickpt=.true.
         end if
         length1 = geti('len1',length1)
         length2 = geti('len2',length2)
         if (find ('action')) goto  5
      end if
      go to 1
 5    continue
        level = 1
        if (length1.eq.0) 
     1    call alert(name,namel,'len path1 is zero',17,level)
        if (length2.eq.0) 
     1    call alert(name,namel,'len path2 is zero',17,level)
        if (npick.eq.0) 
     1    call alert(name,namel,'zero selection',14,level)
      if (.not.fopen(ucon)) then
         call alert(name,namel,'ucon not opened',15,level)
      else if (.not. fopen(upth1)) then
         call alert(name,namel,'urcrd1 not opened',17,level)
      else if (.not.fopen(upth2)) then
         call alert(name,namel,'urcrd2 not opened',17,level)
      end if
        level = 0
c
c     use pick to store masses to zero
c      
      if(pickpt) then
         npick = 0
         do 100 i=1,npt
            if (ipick(i) .eq. 1) then
               npick = npick + 1
               mmm(i) = ptms(i)
               rms_pick(npick) = i
            else
               mmm(i) = 0
            end if
 100    continue
      else
       npick = npt
      end if

        do 800 i=1,length1
        
                call rpath_seq(upth1,i)
        
                do k=1,npt
                  do l=1,3
                        coor2(l,k)=coor(l,k)
                  end do
                end do
                rewind upth2
                do 800 j=1,length2
                  call rpath_seq(upth2,j)
c
                  call rmsd_weight(npick,coor,coor2,rms,.false.,mmm)
                  write(*,*)' i j rms',i,j,rms
800     continue
      stop
      end
