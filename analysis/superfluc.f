             program superfluc

            implicit none
            include 'COMMON/LENGTH.BLOCK'
            include 'COMMON/COORD.BLOCK'
            include 'COMMON/CONNECT.BLOCK'
            include 'COMMON/LINE.BLOCK'
            include 'COMMON/DEBUG.BLOCK'
            include 'COMMON/FREEZ.BLOCK'
            include 'COMMON/CONVERT.BLOCK'
            include 'COMMON/UNITS.BLOCK'
            include 'COMMON/OVERLAP.BLOCK'
                
              
c rms: rms of dynamic structrue respect to reference structure
c rms_ave: rms of dynamic structrue respect to averaged dynamic structure
            integer ucon,urdcrd,uwflu,urcrd,uwrms,uwrms_ave
            integer geti,of,nstep,i,j,ncoor,namel,n,k,jump
            integer npick,level,iselect,ibig
            integer rbin
            integer logic,ipick(maxpt),jpick(maxpt)
            double precision factor,t,resdmass
            double precision sqx(maxpt),sqy(maxpt),sqz(maxpt)
            double precision avesqx(maxpt),avesqy(maxpt),avesqz(maxpt)
            double precision tmpsqx(maxpt),tmpsqy(maxpt),tmpsqz(maxpt)
            double precision ave(3,maxpt)
            double precision tmpx(maxpt),tmpy(maxpt),tmpz(maxpt)
            double precision flucx(maxpt),flucy(maxpt),flucz(maxpt)
            double precision flucr(maxpt),getd,rms_ave,rms,fluc(maxmono)
            character*11 name
            logical find,fopen,fluct,rmsav,rmsre
            logical pickpt,lrot
            double precision Tr,cosine, tmp_mass(maxpt)

                
         
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
c    defalut parameters
            lrot = .false.
            pickpt=.false.
             fluct=.false.
             rmsav=.false.
             rmsre=.false.
             norew=.false.
             ncoor=1
             n = 1
             jump = 1
c initialize ipick to default value
c       do 1 i=1,npt
c        ipick(i)  = 1
1        continue
c initialize jpick to default value
c        do 2 i=1,npt
c        jpick(i)  = 1
c 2      continue

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
               end if
               if (find ('rdyc')) urdcrd=of()
               if (find('wrms')) then 
                uwrms=of()
                rmsre=.true.
                lrot = .true.
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
                   n = geti('nstr',n)
                   jump = geti('jump',jump)
                   if (find('subs')) then
                           call pick(ipick,i)
                           iselect = 0
                           do 9 i=1,npt
                           if (ipick(i).eq.0) then
                            tmp_mass(i) = 0.d0
                           else
                            iselect = iselect + 1
                            tmp_mass(i) = ptms(i)
                           end if
9                       continue
                   end if

               end if
             if (find('action')) then 
               goto 5
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
c ============================================
           if (fluct .or. rmsav) then
            do 25 k=1,n
             rewind urdcrd
             write(*,*)' before a call to rdyncrd '
             call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
             write(*,*)' coor(1) ',(coor(j,1),j=1,3)
               logic=0
            call rmsd_weight(iselect,coor2,coor,rms,.false.,tmp_mass)
c
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
                    fluc(i)=fluc(i)+flucr(j)*ptms(j)
                    resdmass=resdmass+ptms(j)
40              continue
                   else
                  do 45 j=poipt(i-1)+1,poipt(i)
                   fluc(i)=fluc(i)+flucr(j)*ptms(j)
                   resdmass=resdmass+ptms(j)
45             continue
                  end if
                  if (resdmass .eq. 0.) go to 50
                  fluc(i)=fluc(i)/resdmass 
                  write(uwflu,55)i,fluc(i)
55             format(1x,i4,2x,f10.6)
50           continue
              end if                 
                  
c   read dynamics coordinates
c             rewind urdcrd 
             do 85 k=1,n,jump
               if (.not.norew) rewind urdcrd
               call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
              if (rmsre)then
c   overlap dynamics structure and reference structure and calculate rms

             call rmsd_weight(iselect,coor2,coor,rms,.false.,tmp_mass)

                Tr = rotat(1,1)+rotat(2,2)+rotat(3,3)
                cosine = 0.5d0 * (Tr - 1.d0)

              write(stdo,*)'The trace of the rotation matrix is ',Tr
              write(stdo,*)'The cosine of the rotation angle is ',cosine
              t=k

              write(uwrms,82)t,rms,cosine

               end if
        
c   overlap dynamics structrue and average dynamic structure,calculate rms
               if (rmsav) then
           call rmsd_weight(iselect,ave,coor,rms_ave,.false.,tmp_mass) 
                
                write(uwrms_ave,82)t,rms_ave
82          format(1x,f7.2,2(2x,f5.3))

83        format(6(1x,f6.4),/)

84        format(3(1x,f6.4),/)

                endif
                
85        continue

            
             stop
             end








