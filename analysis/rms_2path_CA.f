      program rms_2path_CA
c
c     --------------------------------------------------------------------------:
c
c     This program is designed to analyse two sepatate path bynary files
c     (CHARMM/PATH) and compare the rms for multiple structures.
c
c     Author: V. Zaloj
c     Date: December 20, 1999
c
c     --------------------------------------------------------------------------
C This version only compute CA rms. -Alfredo
C
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
c
c
c
      integer ipick(maxpt) ,npick 
      integer of,geti,nstru,logic
      integer namel,level
      integer ucon,uwrms,urcrd1,urcrd2
      integer i,j,ii,n,nco
      integer nst1s,nst1e,nst2s,nst2e
      integer rmsind(maxmono), nprmono
      double precision rmsall, runi
      double precision velos(3,maxmono),coor2s(3,maxmono),
     & coors(3,maxmono)
      double precision mmm(maxpt)
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
      name='rms_2path'
      namel=9
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
c     open coordinate file #1
         else if (find ('rcc1').and.(inpsty1.eq.0)) then
            urcrd1=of()
            inpsty1=1
c         else if(find ('rdc1').and.(inpsty1.eq.0)) then
c            urcrd1=of()
c            inpsty1=2
         else if(find ('rpc1').and.(inpsty1.eq.0)) then
            urcrd1=of()
            inpsty1=3
c     open coordinate file #2
         else if (find ('rcc2').and.(inpsty2.eq.0)) then
            urcrd2=of()
            inpsty2=1
c         else if(find ('rdc2').and.(inpsty2.eq.0)) then
c            urcrd2=of()
c            inpsty2=2
         else if(find ('rpc2').and.(inpsty2.eq.0)) then
            urcrd2=of()
            inpsty2=3
c     open output file for rms data
         else if(find('wrms')) then  
            uwrms=of()
         endif
      else 
         nst1s=geti('st1s',nst1s)
         nst1e=geti('st1e',nst1e)
         nst2s=geti('st2s',nst2s)
         nst2e=geti('st2e',nst2e)
         nprmono=geti('#rmp',nprmono)
         if (find('frzp')) then
            call pick(ipick,npick)
            pickpt=.true.
         end if
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
      end if
c
c     use pick to store masses to zero
c      
      if(pickpt) then
         npick = 0
         do 100 i=1,npt
            if (ipick(i) .eq. 1) then
                npick = npick + 1
                rms_pick(npick) = i
                mmm(i) = ptms(i)
            else
                mmm(i) = 0.
            end if
 100    continue
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
      n=nst1s
      call getstrcrd(n,urcrd1,inpsty1,velo)
c        
C only write the coordinates for CA:
      nco=0
      do j=1,npt
        if (ptnm(j).eq.'CA')  then
          if (poimon(j).lt.(totmon-1)) then
            do k=poipt(poimon(j))+1,poipt(poimon(j)+1)
              if (ptnm(k).eq.'CA')  then
                nco=nco+1
                rms_pick(nco) = j
                mmm(k) = ptms(k)
                do i=1,3
                  velos(i,nco)=velo(i,k)-velo(i,j)
                enddo
              endif
                mmm(k) = 0
            enddo
          endif
        endif
      enddo
      do j=1,nco
        i = r
        runi=dsqrt(velos(1,i)**2 + velos(2,i)**2 + velos(3,i)**2 )
        velos(1,i)=velos(1,i)/runi
        velos(2,i)=velos(2,i)/runi
        velos(3,i)=velos(3,i)/runi
C        write(6,*) velos(1,i),velos(2,i),velos(3,i)
C        write(6,*) 'nco',nco
      enddo
      n=nst2s
      call getstrcrd(n,urcrd2,inpsty2,coor2)
c
      nco=0
      do j=1,npt
        if (ptnm(j).eq.'CA')  then
          if (poimon(j).lt.(totmon-1)) then
            
            do k=poipt(poimon(j))+1,poipt(poimon(j)+1)
              if (ptnm(k).eq.'CA')  then
                nco=nco+1
                do i=1,3
                  coor2s(i,nco)=coor2(i,k)-coor2(i,j)
                enddo
              endif
            enddo
            
          endif
          
        endif
      enddo
      do i=1,nco
        runi= dsqrt(coor2s(1,i)**2 + coor2s(2,i)**2 + coor2s(3,i)**2 )
        coor2s(1,i)=coor2s(1,i)/runi
        coor2s(2,i)=coor2s(2,i)/runi
        coor2s(3,i)=coor2s(3,i)/runi
      enddo
c     start loop over the structures
c
 210  continue
c
      do 250 j=1,nco
        do 240 i=1,3
           coors(i,j)=coor2s(i,j)
 240    continue
C       write(6,*)"Coor:",coors(1,j),coors(2,j),coors(3,j)
 250  continue
c
c     velo array will be changed by rotation
c  
      call rmsd_weight(nco,velos,coor2s,rmsall,.false.,mmm)
      write(uwrms,'(/10x,a,i4,a,i4)')
     >     "Structure ",nst1s," versus ",nst2s

      call wsortrms(uwrms,nprmono,rmsall,rmsarray,rmsind)

      write(uwrms,*) rmsall
      do 350 j=1,nco
        do 340 i=1,3
           coor2s(i,j)=coors(i,j)
 340    continue
 350  continue
c
      ii=0
      if(nst1e.gt.nst1s) then 
         n=1
         ii=ii+1
         call getstrcrd(n,urcrd1,inpsty1,velo)
         nst1s=nst1s+1
      nco=0
      do j=1,npt
      if (ptnm(j).eq.'CA')  then
      if (poimon(j).lt.(totmon-1)) then
      do k=poipt(poimon(j))+1,poipt(poimon(j)+1)
      if (ptnm(k).eq.'CA')  then
      nco=nco+1
      do i=1,3
      velos(i,nco)=velo(i,k)-velo(i,j)
      enddo
      endif
      enddo
      endif
      endif
      enddo
      do i=1,nco
      runi= dsqrt(velos(1,i)*velos(1,i)+velos(2,i)*velos(2,i)+
     & velos(3,i)*velos(3,i))
      velos(1,i)=velos(1,i)/runi
      velos(2,i)=velos(2,i)/runi
      velos(3,i)=velos(3,i)/runi
      enddo
      endif
c
      if(nst2e.gt.nst2s) then
         n=1
         ii=ii+1
         call getstrcrd(n,urcrd2,inpsty2,coor2)
      nco=0
      do j=1,npt
      if (ptnm(j).eq.'CA')  then
      if (poimon(j).lt.(totmon-1)) then
      do k=poipt(poimon(j))+1,poipt(poimon(j)+1)
      if (ptnm(k).eq.'CA')  then
      nco=nco+1
      do i=1,3
      coor2s(i,nco)=coor2(i,k)-coor2(i,j)
      enddo
      endif
      enddo
      endif
      endif
      enddo
      do i=1,nco
      runi= dsqrt(coor2s(1,i)*coor2s(1,i)+coor2s(2,i)*coor2s(2,i)+
     & coor2s(3,i)*coor2s(3,i))
      coor2s(1,i)=coor2s(1,i)/runi
      coor2s(2,i)=coor2s(2,i)/runi
      coor2s(3,i)=coor2s(3,i)/runi
      enddo
         nst2s=nst2s+1
      endif
c
c     repeat the loop for next structures
c
      if(ii.gt.0)  goto 210
c
      stop
      end
c
c
c
      subroutine getstrcrd(nstru,urcr,inpsty,scrd)
c
c subroutine to read coordinate of the structure
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
c
      integer nstru,urcr,inpsty
      double precision e0,scrd(3,*)

c local
      character*9 name
      integer namel,level,i,j,k
c      
      name  = 'getstrcrd'
      namel = 9
c
      if(nstru.le.0) return
c
      if(inpsty.ne.1.and.inpsty.ne.3) then
         level = 1
         write(stdo,'(10x,a,i5)') "inpsty=",inpsty
         call alert(name,namel,'Wrong inpsty!!!',15,level)
         return
      endif
c
c     CHARMM format
c
      if(inpsty.eq.1) then
         call getcrd(urcr,'CHARM')
         do 150 j=1,npt
            do 140 i=1,3
               scrd(i,j)=coor(i,j)
 140        continue
 150     continue
         return
      endif
c
c     PATH format
c
      if(inpsty.eq.3) then
         do 300 n=1,nstru
            read(urcr,err=999,end=999) e0,((scrd(j,i),i=1,npt),j=1,3)
            write(stdo,400) e0
 400        format(1x,' *** PATH FILE : ENERGY = ',e15.8)
 300     enddo 
         return
 999     continue
         level = 1
         call alert(name,namel,'Error while reading Path file',29,level)
         return
      end if
c
      return
      end

      Subroutine wsortrms(uwrms,nprmono,rmsall,rmsarr,rmsind)
      integer i,j,j1,j2,k,imax,imin,uwrms
      integer nprmono,nmax,rmsind(*)
      double precision rmsall,rmsmax,rmsmin,rmsarr(*)
       
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
c
      nmax=totmon
      do 10 i=1,nmax
         rmsind(i)=i
 10   continue
c
      kmax=(nmax-1)/2 
      kmax=min(nprmono,kmax)
c
      nprmono=min(totmon,max(1,nprmono))
c
      do 100 i=1,kmax
         imin=i
         imax=i
         rmsmax=rmsarr(i)
         rmsmin=rmsarr(i)
         j1=i+1
         j2=nmax-i+1
         do 50 j=j1,j2
            l=rmsind(j)
            if(rmsmax.lt.rmsarr(l)) then
                imax=j
                rmsmax=rmsarr(l)
            endif
            if(rmsmin.gt.rmsarr(l)) then
                imin=j
                rmsmin=rmsarr(l)
            endif
 50      continue
         k=rmsind(imax)
         rmsind(imax)=rmsind(i)
         rmsind(i)=k
         k=rmsind(imin)
         rmsind(imin)=rmsind(j2)
         rmsind(j2)=k              
 100  continue
c
C      write(uwrms,*) rmsall 
      write(uwrms,'(/10x,a,f12.6/)') "Total  R M S = ",rmsall 
c
      do 300 i=1,nprmono
         j=rmsind(i)
         write(uwrms,'(2x,2i8,f12.6)') i,j,rmsarr(j)
 300  continue

      return
      end
