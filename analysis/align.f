           program align

C
C  This program takes three crd files as an input, 
C  First two structures share same connectivity file, third structure is an extension 
C  of the 2nd structure.
C 
C  This program finds a transformation of the second structure to the first one
C  and then applies this transformation to the third structure.
C
C EXAMPLE: Structures 1,2 are in C-alpha representation of a protein in two different
C configurations. The third structure is the all-atom representation of the 2nd structure. You
C might want to get a quess for all-atom description of the 1st structure. You can use this program and
C aply the rigid body transformation that takes structure 2 to structure 1 on the all-atomistic structure 3.
C
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
c          integer ipick(maxpt) ,npick 
           integer ipick(maxpt),jpick(maxpt)
           integer of,geti,nstru,logic
           integer namel,i,level,npick,iselect
           integer urcrd1,urcrd2,urcrd3,ucon1,ucon2,uwcrd
           double precision rmsall,Cntr1(3),Cntr2(3),Cntr3(3)
           double precision coor3(3,maxpt),tmp_mass(maxpt)
           character*3 name
           logical find,fopen
           logical pickpt
           data urcrd1,urcrd2,urcrd3,ucon1,ucon2,uwcrd/6*99/
c          double precision TR,cosine
           double precision Tr,cosine,values(3),vec(3),rotat(3,3)
          
           logical CAonly 
           integer l,j,k,loop,npt3
           double precision tmp(3), rms

           CAonly= .false.
           norew = .false.

           stdi=5
           stdo=6
           totmon=0
           npt=0
           name='rms'
           namel=3
           logic=1
*  open junk file for rline
            jnkf=25
            open(unit=jnkf,status='scratch')

1       continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('con1')) then
                ucon1=of()
                call rconn(ucon1)
               end if
               if (find ('con2')) ucon2=of()

               if (find ('rcr1')) urcrd1=of()
               if (find ('rcr2')) urcrd2=of()
               if (find ('rcr3')) urcrd3=of()
               if (find('wcrd')) uwcrd = of()
              else 
               if (find('CAon')) CAonly=.true.
               if (find ('action')) goto  5
             end if
             go to 1

5            continue

          if (.not. fopen(ucon1)) then
             level=1
             call alert(name,namel,'ucon1 not opened',15,level)
            else if (.not. fopen(ucon2)) then
            level=1
             call alert(name,namel,'ucon2 not opened',15,level)
            else if (.not. fopen(urcrd1)) then
             level=1
             call alert(name,namel,'urcrd1 not opened',16,level)
            else if (.not. fopen(urcrd2)) then
             level=1
             call alert(name,namel,'urcrd2 not opened',16,level)
            else if (.not. fopen(urcrd3)) then
             level=1
             call alert(name,namel,'urcrd3 not opened',16,level)
            else if (.not. fopen(uwcrd)) then
             level=1
             call alert(name,namel,'uwcrd not opened',16,level)
           end if

           if(CAonly) then
             j=0
             do 105 k=1,npt
               if (ptnm(k).ne.'CA  ')  then
                   ptms(k)=0.0d0
               end if
105          continue
           end if
        
            call getcrd(urcrd2,'CHARM')
            do i=1,npt
              do l = 1,3
                coor2(l,i) = coor(l,i)
              end do
            end do

            call Find_Center(npt,coor2,Cntr2)

            call rconn(ucon2)
            npt3=npt
            call getcrd(urcrd3,'CHARM')
            do i=1,npt3
              do l = 1,3
                coor3(l,i) = coor(l,i)
              end do
            end do

            call Find_Center(npt3,coor3,Cntr3) 
            
            call rconn(ucon1)
            call getcrd(urcrd1,'CHARM')

            call Find_Center(npt,coor,Cntr1)

            do l =1,3
              do i=1,npt        
                 coor2(l,i) = coor2(l,i) - Cntr2(l)
                 coor(l,i)  = coor(l,i)  - Cntr1(l)
              end do
              do i=1,npt3
                 coor3(l,i) = coor3(l,i) - Cntr3(l)
              end do
            end do
           
            call rmsd_weight(npt,coor(1,1),coor2(1,1),rms,.false.,
     1        ptms)
        
            do k=1,npt3
              tmp(1)=rotat(1,1)*coor3(1,k)+rotat(1,2)*coor3(2,k)
     1              +rotat(1,3)*coor3(3,k)
              tmp(2)=rotat(2,1)*coor3(1,k)+rotat(2,2)*coor3(2,k)
     1              +rotat(2,3)*coor3(3,k)
              tmp(3)=rotat(3,1)*coor3(1,k)+rotat(3,2)*coor3(2,k)
     1              +rotat(3,3)*coor3(3,k)

              coor3(1,k) = tmp(1)
              coor3(2,k) = tmp(2)
              coor3(3,k) = tmp(3)
            end do
 
            do i=1,npt3
              do l = 1,3
                coor(l,i) = coor3(l,i) + Cntr1(l)
              end do
            end do 
         
            call rconn(ucon2)
            call putcrd(uwcrd,'CHARM') 
                
         stop
         end
