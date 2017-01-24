           program superrms
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/DEBUG.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/FREEZ.BLOCK'
           include 'COMMON/CONVERT.BLOCK'
           include 'COMMON/CCRD.BLOCK'
           include 'COMMON/OVERLAP.BLOCK'
c          integer ipick(maxpt) ,npick 
           integer ipick(maxpt),jpick(maxpt)
           integer of,geti,nstru,logic
           integer namel,i,level,npick,iselect
           integer urcrd,ucon,urdcrd,uwrms,urcor
           integer rbin
           double precision  mmm(maxpt)
           double precision rms
           character*3 name
           logical find,fopen,rdyncor
           logical pickpt
           data ucon,urcrd,urdcrd/3*99/
c          double precision TR,cosine
           double precision Tr,cosine,values(3),vec(3)
          
           logical CAonly, self_rms 
           integer k,loop2

           CAonly= .false.
           norew = .true.
           self_rms = .false.
           lpstr = 1

           stdi=5
           stdo=6
           rbin = 1

           totmon=0
           iselect = 0
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
1           continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('conn')) then  
                ucon=of()
* get connectivity 
                call rconn(ucon)
               end if
               if (find ('rcrd')) urcrd=of()
               if (find ('rdyc')) then
                  urdcrd=of()
                  rdyncor=.true.
               end if
               if (find ('rcor')) urcor=of()
               if(find('wrms'))  uwrms=of()
              else 
               if (find('CAon')) CAonly=.true.
               nstru=geti('#str',nstru)
              if (find('pick')) then
               call pick(ipick,npick)
                iselect = 0
                do 9 i=1,npt
                if (ipick(i).ne.0) then
                        iselect = iselect + 1
                        mmm(i) = ptms(i)
                        rms_pick(iselect) = i
                else
                        mmm(i) = 0
                end if
9           continue
            end if
               if (find('self')) then
                 self_rms = .true.
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
            else if ((.not.rdyncor) .and. (.not. fopen(urcor))) then
             level=1
             call alert(name,namel,'urcor not opened',16, level)
           end if

           if (self_rms) rdyncor = .false.
             
* initialize the vector nofreez      
            inofrz=npt
            do 8 i=1,npt
             nofreez(i)=i
8           continue
   
           if(CAonly) then
             call getcrd(urcrd,'CHARM')
             rewind urcrd
             iselect=0
             do 105 k=1,npt
               if (ptnm(k).eq.'CA  ')  then
                 iselect = iselect + 1
                 rms_pick(iselect) = k
                 mmm(k) = ptms(k)
               else
                 mmm(k) = 0
               end if

C    CG models might store secondary structure information 
C    in the more(k) array
C    coil = 0, alpha helix = 1, beta sheet = 2, 
C    beta turn = 3, pi helix = 4
C              if (more(k).eq. 0.00)  then
C                ptms(k)=0.d0
C              end if
105          continue
             write(*,*)' iselect = ',iselect
           else if (iselect .eq. 0) then
            iselect = npt
            do i = 1,npt
             mmm(i) = ptms(i)
            end do
           end if
        
        if (rdyncor)  then
           call getcrd(urcrd,'CHARM')
           do i=1,npt
             coor2(1,i)=coor(1,i)
             coor2(2,i)=coor(2,i)
             coor2(3,i)=coor(3,i)
           end do

        
           do 124 loop=1, nstru
              call rdyncrd(urdcrd,loop,inofrz,nofreez,rbin)
              call rmsd_weight(npt,coor,coor2,rms,.false.,mmm)

              write(uwrms,11) rms
11            format(f10.7,2x,f10.7)
124         continue

          else if (self_rms) then 
            do loop = 1, nstru
              rewind urdcrd
              norew = .false.
              call rdyncrd(urdcrd,loop,inofrz,nofreez,rbin)
             
              do i=1,npt
                coor2(1,i)=coor(1,i)
                coor2(2,i)=coor(2,i)
                coor2(3,i)=coor(3,i)
              end do

              rewind urdcrd
              norew = .true.
              do loop2=1, nstru
                call rdyncrd(urdcrd,loop2,inofrz,nofreez,rbin)
                call rmsd_weight(npt,coor,coor2,rms,.false.,mmm)

                write(uwrms,11) rms
              end do
            end do

          else
            call getcrd(urcrd,'CHARM')
            do 13 i=1,npt
                coor2(1,i)=coor(1,i)
                coor2(2,i)=coor(2,i)
                coor2(3,i)=coor(3,i)
13          continue
            call getcrd(urcor,'CHARM')
            
            call rmsd_weight(npt,coor,coor2,rms,.false.,mmm)
            
            write(uwrms,11) rms
          endif
                
         stop
         end

