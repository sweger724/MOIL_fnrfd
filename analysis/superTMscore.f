           program superTMscore
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
           include 'COMMON/OVERLAP.BLOCK'
c          integer ipick(maxpt) ,npick 
           integer ipick(maxpt)
           integer of,geti,nstru,logic
           integer namel,i,level,npick,iselect
           integer urcrd,ucon,urdcrd,uwrms,urcor
           integer rbin
c          double precision rms(maxmono)
           double precision rms
           character*3 name
           logical find,fopen,rdyncor
           logical pickpt
           data ucon,urcrd,urdcrd/3*99/
c          double precision TR,cosine
           double precision Tr,cosine,values(3),vec(3)
          
           logical CAonly 
           integer j,k, loop

           double precision x1(maxpt),y1(maxpt),z1(maxpt)
           double precision x2(maxpt),y2(maxpt),z2(maxpt)
           integer n1(maxpt), n2(maxpt)
           double precision TM, Lcomm,npt2

           CAonly= .false.
           norew = .false.
           lpstr = 1

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
                        rms_pick(iselect) = i
                end if
9           continue
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

             
* initialize the vector nofreez      
            inofrz=npt
            do 8 i=1,npt
             nofreez(i)=i
8           continue
   
           if(CAonly) then
             iselect=0
             do 105 k=1,npt
               if (ptnm(k).ne.'CA  ')  then
                   ptms(k) = 0
                   rms_pick(iselect) = k
               else
                   iselect = iselect + 1
               end if
105          continue
           end if
        
        if (rdyncor)  then
            call getcrd(urcrd,'CHARM')
           do 124 loop=1, nstru
              rewind urcrd
                do 10 i=1,npt
                   coor2(1,i)=coor(1,i)
                   coor2(2,i)=coor(2,i)
                   coor2(3,i)=coor(3,i)
10              continue
              call rdyncrd(urdcrd,loop,inofrz,nofreez,rbin)
              rewind urdcrd
              
              do i=1, npt
                x1(i) = coor(1,i)
                y1(i) = coor(2,i)
                z1(i) = coor(3,i)
                x2(i) = coor2(1,i)
                y2(i) = coor2(2,i)
                z2(i) = coor2(3,i)
                n1(i) = i
                n2(i) = i
C               write(6,*)"x1:",x1(i),y1(i),z1(i)
              enddo
              
              call TMscore(npt,x1,y1,z1,n1,npt,x2,y2,z2,n2,TM,rms,Lcomm)
C             call rmsd_weight(npt,coor,coor2,rms,.false.,ptms)

              write(uwrms,11) TM
11            format(f10.7,2x,f10.7)
124         continue
          else
            call getcrd(urcrd,'CHARM')
            do 13 i=1,npt
                coor2(1,i)=coor(1,i)
                coor2(2,i)=coor(2,i)
                coor2(3,i)=coor(3,i)
13          continue
            call getcrd(urcor,'CHARM')
            
            call rmsd_weight(iselect,coor,coor2,rms,.false.,ptms)
            
            write(uwrms,11) rms
          endif
                
         stop
         end

