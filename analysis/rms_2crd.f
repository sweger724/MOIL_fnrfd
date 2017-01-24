      program rms_2crd
      implicit none
c
c     --------------------------------------------------------------------------:
c
c     Compute rms of two structures (crd) format. The same connectivity file
c     We compute Kabsch matrix for a selected set of atom (ipick1) and compute
c     the rms for another set (ipick2)
c
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
      integer ipick(maxpt),ipick1(maxpt),ipick2(maxpt)
      integer npick1,npick2,npt2
      integer of,geti,nstru,logic
      integer namel,level
      integer ucon,uwrms,urcrd1,urcrd2
      integer i,j,ii,n, inpsty1, inpsty2, nprmono
      integer nst1s,nst1e,nst2s,nst2e
      double precision rmsall, tmp(3),mmm(maxpt), rms
      double precision mm2(maxpt)
      character*9 name
      logical find,fopen
      logical pickpt
      data ucon,uwrms,urcrd1,urcrd2/4*99/

      norew = .false.
      lpstr = 1

      stdi=5
      stdo=6
      totmon=0
      npt=0
      name='rms_2crd'
      namel=8
      logic=1
c     open junk file for rline
c     
      jnkf=25
      open(unit=jnkf,status='scratch')
c
c     defalt parameters
c
      inpsty1=0
      inpsty2=0
c
      nst1s=1
      nst1e=1
      nst2s=1
      nst2e=1
      nprmono=1
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
         else if (find ('rcc1')) then
            urcrd1=of()
         else if(find ('rcc2')) then
            urcrd2=of()
         endif
      else 
c
c the first select is for the group of atoms
         if (find('sel1')) then
            call pick(ipick,npick1)
            pickpt=.true.
            npick1 = 0
            do 2 i=1,npt
                if (ipick(i).ne.0) then
                        npick1 = npick1 + 1
                        ipick1(npick1) = i
                        rms_pick(npick1) = i
                        mmm(i) = ptms(i)
                else
                        mmm(i) = 0.d0
                end if
2           continue
         end if
         if (find('sel2')) then
           call pick(ipick,npick2)
           npick2 = 0
           do 3 i=1,npt
                if (ipick(i).ne.0) then
                 npick2 = npick2 + 1
                 ipick2(npick2) = i
                 mm2(i) = ptms(i)
                else
                 mm2(i)= 0.d0
                end if
3          continue
         end if
             
                
         end if
         if (find ('action')) goto  5
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
      end if
c
c     use pick to store masses to zero
c      
        call getcrd(urcrd1,'CHARM')
        
      do 250 j=1,npt
        do 240 i=1,3
           coor2(i,j)=coor(i,j)
 240    continue
 250  continue
c
        call getcrd(urcrd2,'CHARM')
        
        write(*,*)' npick1 npick2 coor coor2 '
        npt2 = npt
        call rmsd_weight2(npt,npt2,coor,coor2,rms,.false.,mmm,mm2)
       
      stop
      end
