         program fluctuation
         implicit none

         include 'COMMON/LENGTH.BLOCK'
         include 'COMMON/COORD.BLOCK'
         include 'COMMON/CONNECT.BLOCK'
         include 'COMMON/LINE.BLOCK'
         include 'COMMON/DEBUG.BLOCK'
         include 'COMMON/FREEZ.BLOCK'
         include 'COMMON/CONVERT.BLOCK'
         include 'COMMON/OVERLAP.BLOCK'
         include 'COMMON/UNITS.BLOCK'
           
c rms: rms of dynamic structrue respect to reference structure
c rms_ave: rms of dynamic structrue respect to averaged dynamic structure
         integer ucon,urdcrd,uwflu,urcrd,uwrms,uwrms_ave,uwcrd
         integer geti,of,nstep,i,j,ncoor,namel,n,k
         integer npick,level
         integer logic,ipick(maxpt)
         integer rbin
         double precision dt ,factor,t,resdmass
         double precision sqx(maxpt),sqy(maxpt),sqz(maxpt)
         double precision avesqx(maxpt),avesqy(maxpt),avesqz(maxpt)
         double precision tmpsqx(maxpt),tmpsqy(maxpt),tmpsqz(maxpt)
         double precision ave(3,maxpt)
         double precision tmpx(maxpt),tmpy(maxpt),tmpz(maxpt)
         double precision flucx(maxpt),flucy(maxpt),flucz(maxpt)
         double precision flucr(maxpt),getd,rms_ave,rms,fluc(maxmono)
         double precision mmm(maxpt)
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

         npt=0
         name='fluctuation'
         namel=11
c  open junk file for rline
c
         jnkf=25
         open(unit=jnkf,status='scratch')
c    defalt parameters
         pickpt=.false.
          fluct=.false.
          rmsav=.false.
          rmsre=.false.
          norew=.false.
          ncoor=1
1         continue
          call rline(name,namel,stdi)
          if (find ('norw')) norew=.true.
          if (find ('file')) then
            if (find ('conn')) then
             ucon=of()
c  get connectivity
             call rconn(ucon)
             do i=1,npt
              rms_pick(i) = i
              mmm(i) = 1
             end do
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
            end if
            if (find('wcrd')) uwcrd = of()
            if (find ('rdyc')) urdcrd=of()
            if (find('wrms')) then 
             uwrms=of()
             rmsre=.true.
            end if
            if (find('wave')) then
              uwrms_ave=of()
              rmsav=.true.
            end if
            if (find('wflu')) then
              uwflu=of()
              fluct=.true.
            end if
          else
            ncoor=geti('#crd',ncoor)
            nstep=geti('#ste',nstep)
            n=nstep/ncoor
            dt=(getd('step',(dt)))
            if (find('pick')) then
             call pick(ipick,npick)
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
          else if (.not. fopen(uwflu)) then
            if (fluct) then
             level=1
             call alert(name,namel,'uwflu not opened',16,level)
            end if
          end if
c  initialize the vector nofreez and coordinates
          inofrz=npt
          do 8 i=1,npt
           nofreez(i)=i
           tmpx(i)=0.d0
           tmpy(i)=0.d0
           tmpz(i)=0.d0
           tmpsqx(i)=0.d0
           tmpsqy(i)=0.d0
           tmpsqz(i)=0.d0
8         continue
c
          factor=1.0/dble(n)
c ============================================
        if (pickpt) then
         npick = 0
         do 1000 i=1,npt
          if(ipick(i) .ne. 0) then
           npick = npick + 1
           rms_pick(npick) = i
          else
           mmm(i) = 0.
          end if
1000     continue
        end if
c ============================================
        if (fluct .or. rmsav) then
         do 25 k=1,n
          call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
            logic=0
            call rmsd_weight(npt,coor2,coor,rms,.false.,mmm)
            do 20 i=1,npt
             tmpx(i)=tmpx(i)+coor(1,i)    
             tmpy(i)=tmpy(i)+coor(2,i)
             tmpz(i)=tmpz(i)+coor(3,i)
             sqx(i)=coor(1,i)*coor(1,i)
             sqy(i)=coor(2,i)*coor(2,i)
             sqz(i)=coor(3,i)*coor(3,i)
             tmpsqx(i)=tmpsqx(i)+sqx(i)
             tmpsqy(i)=tmpsqy(i)+sqy(i)
             tmpsqz(i)=tmpsqz(i)+sqz(i)
20          continue
            if (.not.norew) rewind urdcrd
25        continue
         do 35 i=1,npt
           ave(1,i)=factor*tmpx(i)
           ave(2,i)=factor*tmpy(i)
           ave(3,i)=factor*tmpz(i) 
           if (fluct) then
             avesqx(i)=factor*tmpsqx(i)    
             avesqy(i)=factor*tmpsqy(i)    
             avesqz(i)=factor*tmpsqz(i)    
             flucx(i)=avesqx(i)- ave(1,i)*ave(1,i)                      
             flucy(i)=avesqy(i)- ave(2,i)*ave(2,i)                      
             flucz(i)=avesqz(i)- ave(3,i)*ave(3,i)                      
             flucr(i)=flucx(i)+flucy(i)+flucz(i)
            end if
35        continue
         end if
         if (fluct) then
c calculate fluctuation for each residue
             do 50 i=1,totmon
               fluc(i)=0.d0
               resdmass=0.d0
               if (i .eq. 1) then
                do 40 j=1,poipt(1)
                 if (mmm(j).ne.0) then
                 fluc(i)=fluc(i)+flucr(j)*ptms(j)
                 resdmass=resdmass+ptms(j)
                 end if
40              continue
                else
               do 45 j=poipt(i-1)+1,poipt(i)
                if (mmm(j).ne.0) then
                fluc(i)=fluc(i)+flucr(j)*ptms(j)
                resdmass=resdmass+ptms(j)
                end if
45             continue
               end if
               if (resdmass .eq. 0.) go to 50
               fluc(i)=fluc(i)/resdmass 
               write(uwflu,55)i,fluc(i)
55             format(1x,i4,2x,f10.6)
50           continue
           end if                 
               
c   read dynamics coordinates
          rewind urdcrd 
          do 85 k=1,n
            if (.not.norew) rewind urdcrd
            call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
            if (rmsre) then
c   overlap dynamics structure and reference structure and calculate rms
             call rmsd_weight(npt,coor2,coor,rms,.false.,mmm) 
             t=k*ncoor*dt
             write(uwrms,82)t,rms
            end if

c   overlap dynamics structrue and average dynamic structure,calculate rms
            if (rmsav) then
             call rmsd_weight(npt,ave,coor,rms_ave,.false.,mmm) 
             t=k*ncoor*dt
             write(uwrms_ave,82)t,rms_ave
            end if
82          format(1x,f8.2,2x,f10.7)
85        continue
c         write(*,*) ' writing ',npick,' coordinates' 
c         npt = npick
c         do i=1,npt
c          do j=1,3
c          coor(j,i) = ave(j,i)
c          end do
c         end do
c	call putcrd(uwcrd,'CHARM')
          stop
          end
