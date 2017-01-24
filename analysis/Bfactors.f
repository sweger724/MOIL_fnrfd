           program Bfactors

           implicit none
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/CCRD.BLOCK'
           
           integer of,geti,nstru,lstart,lend,level
           integer namel,i,loop,k,j,l
           double precision Bx(maxpt),By(maxpt),Bz(maxpt)
           integer urcrd1,urcrd2,ucon,upath,uwpth,beg,end
           character*8 name
           logical find,fopen
           data ucon,urcrd1,urcrd2,upath,uwpth/5*99/
          
           double precision rx,ry,rz,r1,r2,e0
           
           stdi=5
           stdo=6
           totmon=0
           npt=0
           name='Bfactors'
           namel=8
           nstru=1

*  open junk file for rline
*
           jnkf=25
           open(unit=jnkf,status='scratch')

1         continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('conn')) then  
                  ucon=of()
                  call rconn(ucon)
               end if
               if (find ('rcor')) urcrd1=of()
               if (find ('rpth')) upath=of()
            else 
               nstru=geti('#str',nstru)
               if (find ('action')) goto  5
            end if
            go to 1
5         continue

          if (.not. fopen(ucon)) then
             level=1
             call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd1)) then
             level=1
             call alert(name,namel,'urcrd1 not opened',16,level)
            else if (.not. fopen(upath)) then
             level=1
             call alert(name,namel,'upath not opened',16, level)
           end if

           call getcrd(urcrd1,'CHARM')
           do 10 i=1,npt
             do 10 l=1,3
                coor2(l,i)=0.d0
10         continue

#          beg=70
#          end=140

          beg=351
          end=700

          do loop=1, nstru
               k=(loop-1)*npt
               read(upath,err=999,end=999)e0,((coor(j,i),i=1,npt),j=1,3)
                if (nstru.gt.0) then
                   write(stdo,100)loop,e0
                end if
               if (loop .lt. end .and. loop .ge. beg) then
                 do i=1,npt
                     coor2(1,i) = coor2(1,i) + coor(1,i)
                     coor2(2,i) = coor2(2,i) + coor(2,i)
                     coor2(3,i) = coor2(3,i) + coor(3,i)
                 end do
               end if
           end do

           do i =1, npt
             coor2(1,i) = coor2(1,i)/(end-beg+1)
             coor2(2,i) = coor2(2,i)/(end-beg+1)
             coor2(3,i) = coor2(3,i)/(end-beg+1)
           end do

           rewind upath

           do 20 loop=1, nstru
                k=(loop-1)*npt
               read(upath,err=999,end=999)e0,((coor(j,i),i=1,npt),j=1,3)
                if (nstru.gt.0) then
                   write(stdo,100)loop,e0
100   format(1x,' *** READING PATH FILE ',i5,' ENERGY = ',e15.8)
                end if
               if (loop .lt. end .and. loop .ge. beg) then
                 do i=1,npt
                     Bx(i) = Bx(i) + (coor(1,i) - coor2(1,i))**2
                     By(i) = By(i) + (coor(2,i) - coor2(2,i))**2
                     Bz(i) = Bz(i) + (coor(3,i) - coor2(3,i))**2
                 end do
               end if
20         continue

           do i=1,npt
             Bx(i) = (Bx(i) + By(i) + Bz(i))/(end-beg+1)
             write(6,*) "bfactor: ",moname(poimon(i))," ",ptnm(i),Bx(i)
           end do

           write(6,*) 'DONE'
        stop

999   continue
      level = 1
      call alert(name,namel,'Error while reading Path.',29,level)
      end
