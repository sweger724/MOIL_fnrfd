           program rms_residue

           implicit none
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/DEBUG.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/FREEZ.BLOCK'
           include 'COMMON/CONVERT.BLOCK'
           include 'COMMON/CCRD.BLOCK'
c           include 'COMMON/OVRLAP.BLOCK'
           integer ipick(maxpt) ,npick 
           integer of,geti,nstru,logic
           integer namel,i,level,j,k
           integer urcrd,ucon,urdcrd,uwrms,urcor
           integer rbin
           double precision rms(maxmono)
           double precision rms_ave(maxmono)
           double precision coor_ave(3,maxpt)
           double precision mmm(maxpt)
           character*3 name
           logical find,fopen,rdyncor
           logical pickpt,mean
           data ucon,urcrd,urdcrd/3*99/
           data rms_ave/maxmono*0.d0/

           norew = .false.
           lpstr = 1
           call vinit(coor_ave, 0.d0, 3*maxpt)
           stdi=5
           stdo=6
           rbin = 1

           totmon=0
           npt=0
           name='rms'
           namel=3
           logic=1
*  open junk file for rline
*
            jnkf=25
            open(unit=jnkf,status='scratch')
* defalt parameters
            nstru=1
            pickpt=.false.
            rdyncor=.false.
            mean = .false.
1           continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('conn')) then  
                ucon=of()
* get connectivity 
                call rconn(ucon)
               end if
               if (find ('rcrd')) then
                urcrd=of()
* read reference coordinate
                 call getcrd(urcrd,'CHARM')
                 do 10 i=1,npt
                   coor2(1,i)=coor(1,i)
                   coor2(2,i)=coor(2,i)
                   coor2(3,i)=coor(3,i)
10               continue 
               end if
               if (find ('rdyc')) then
               urdcrd=of()
               rdyncor=.true.
               end if
               if (find ('rcor')) urcor=of()
               if(find('wrms'))  uwrms=of()
              else 
               nstru=geti('#str',nstru)
               if (find('mean')) mean = .true.
               if (find('pick')) then
                call pick(ipick,npick)
                pickpt=.true.
               end if
               if (find ('action')) goto  5
             end if
             go to 1
5            continue

          if (.not. fopen(ucon)) then
             level=1
             call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd)) then
             level=1
             call alert(name,namel,'urcrd not opened',16,level)
            else if (.not. fopen(uwrms)) then
             level=1
             call alert(name,namel,'uwrms not opened',16,level)
            else if (rdyncor .and. (.not. fopen(urdcrd))) then
             level=1
             call alert(name,namel,'urdcrd not opened',16, level)
           end if

             
           if(pickpt) then
            npick = 0
            do 1000 i=1,npt
             if (ipick(i) .eq. 0) then
              mmm(i) = 0.d0
             else
              npick = npick + 1
              mmm(i) = ptms(i)
              end if
1000       continue
            else
             npick = npt
            end if
* initialize the vector nofreez      
            inofrz=npt
            do 6 i=1,npt
             nofreez(i)=i
6           continue
   
         do 7 i=1,nstru
                rewind urdcrd
                call rdyncrd(urdcrd,i,inofrz,nofreez,rbin)
                call rmsd_weight(npick,coor2,coor,rms,.false.,mmm)
                call vadd(coor_ave,coor,3*npt)
7        continue
         call scalar_mu_vec1(1.d0/nstru, coor_ave, 3*npt)

                
         do 8 i=1,nstru
                rewind urdcrd
                call rdyncrd(urdcrd,i,inofrz,nofreez,rbin)
                call rmsd_weight(npick,coor_ave,coor,rms,.false.,mmm)
                call vadd(rms_ave,rms,totmon)
8        continue
         call scalar_mu_vec1(1.d0/nstru, rms_ave, totmon)
             do 50 i=1,totmon
              if (rms_ave(i) .eq. 0.) go to 50
              write(uwrms,11) i,rms_ave(i)
11            format(1x,i7,2x,f10.7)
50           continue
             stop
             end
             
