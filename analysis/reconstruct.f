           program reconstruct
        
C       This program takes in all atomistic representation
C       of two different conformations (CRD files) of a same molecule 
c       (single wcon file). It also takes a coarse grained (CA particles only)
c       representation of a trajectory between these two endpoints.

c       On output it returns all-atomistic representation of the trajectory
c       by using an algorithm described in 

c       Peter Majek, Harel Weinstein, Ron Elber: Pathways of conformational 
c       transitions in proteins in Coarse-Graining of Condensed Phase and 
c       Biomolecular Systems, ed. Gregory Voth (2008), in section VI.


           implicit none
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/CCRD.BLOCK'
           
           integer of,geti,nstru,lstart,lend
           integer namel,i,loop,k,j,l, extra
           integer urcrd1,urcrd2,ucon,upath,uwpth
           character*11 name
           logical find,fopen
           data ucon,urcrd1,urcrd2,upath,uwpth/5*99/
          
           integer CAindex(maxmono),level, monDecrease,NN
           double precision CA(3,maxmono), CA2(3,maxmono)
           double precision e0,r(3,maxpt*lgrid) 
C       if reverse=true then the ordering of frames is 1,2,...,n-1,n,n-1,...,2
           logical reverse

           stdi=5
           stdo=6
           totmon=0
           npt=0
           name='reconstruct'
           namel=11
           nstru=1

           reverse=.false.
           
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
               if (find ('rcr1')) urcrd1=of()
               if (find ('rcr2')) urcrd2=of()
               if (find ('rpth')) upath=of()
               if (find ('wpth')) uwpth=of()
            else 
               nstru=geti('#str',nstru)
               if (find('join')) reverse=.true.
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
            else if (.not. fopen(urcrd2)) then
             level=1
             call alert(name,namel,'urcrd2 not opened',16,level)
            else if (.not. fopen(uwpth)) then
             level=1
             call alert(name,namel,'uwpth not opened',16,level)
            else if (.not. fopen(upath)) then
             level=1
             call alert(name,namel,'upath not opened',16, level)
           end if

             
             j=0
             monDecrease = 0
             do 105 k=1,npt
               if (ptnm(k).eq.'CA  ')  then
                   j=j+1
                  CAindex(j) = k
               end if
               if (ptnm(k).eq.'OX2 ')  then
                  j=j+1
                  CAindex(j) = k
               end if
               if (ptnm(k).eq.'HX2 ')  then
                  j=j+1
                  CAindex(j) = k
                  monDecrease = monDecrease + 1
                end if
105          continue

           do 20 loop=1, nstru
                k=(loop-1)*npt
                write(6,*)"CAs: ",totmon-monDecrease,j
                NN = totmon - monDecrease
                read(upath,err=999,end=999)e0,((coor(j,i),i=1,NN),j=1,3)
                if (nstru.gt.0) then
                   write(stdo,100)loop,e0
100   format(1x,' *** READING PATH FILE ',i5,' ENERGY = ',e15.8)
                end if
                
                extra = 0
                do i=1,totmon
                  if(ptnm(CAindex(i)).ne.'HX2 ') then
                    do l=1,3
                        r(l,CAindex(i)+k)=coor(l,i-extra)
                    end do
C                   write(6,*)"R(CAindex(i)+k) = ",CAindex(i),r(1,CAindex(i)+k)
                  else
                    extra = extra + 1
                  endif
                end do
20         continue

           if(reverse) then
             do loop=1,nstru-2
                k = (nstru+loop-1)*npt
                j = (nstru-loop-1)*npt
                do i=1,totmon
                  do l =1,3
                    r(l,CAindex(i)+k)=r(l,CAindex(i)+j)
                  end do 
                end do
             end do
           endif
           write(6,*) 'path read'
           call getcrd(urcrd2,'CHARM')
           do 10 i=1,npt
             do 10 l=1,3
               coor2(l,i)=coor(l,i)
10         continue        
           write(6,*) 'coor2 read'
           call getcrd(urcrd1,'CHARM')
           write(6,*) 'coor1 read'
           if (reverse) then
             lstart=1
             lend=2*nstru-2
           else
             lstart=1
             lend=nstru
           endif

           call Align(r(1,1), coor(1,1))
           call Align(r(1,(nstru-1)*npt+1), coor2(1,1))


           do 124 loop= 1, nstru
                k=(loop-1)*npt
                call Add_atoms5(loop,nstru,r(1,k+1),CAindex(1))
124        continue
           
           if (reverse) then
             do loop= nstru+1, 2*nstru-2
                k=(loop-1)*npt
                call Add_atoms5(2*nstru-loop,nstru,r(1,k+1),CAindex(1))
             end do
           endif

C          Write result to a file
           do 52 loop=lstart,lend
                   k = (loop-1)*npt
                   write(uwpth) 0.d0,
     $                  (r(1,j),j=k+1,k+npt),
     $                  (r(2,j),j=k+1,k+npt),
     $                  (r(3,j),j=k+1,k+npt)
52         continue
        stop

999   continue
      level = 1
      call alert(name,namel,'Error while reading Path.',29,level)
      end

C************************************************************************
C************************************************************************

        subroutine Add_atoms5(index,igrid,result,CAindex)

        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/CCRD.BLOCK'

        integer index,igrid, CAindex(maxmono)
        double precision result(3,maxpt), sc1, sc2
        double precision rel1(3), rel2(3)
        double precision CenterR(3),Center1(3),Center2(3)
        double precision A3(3,5),B3(3,5),R3(3,5)
        double precision RotMatrix1(3,3),RotMatrix2(3,3)

        integer i,j,k,l,jjj

          sc2=(index-1.0)/(igrid-1.0)
          sc1=1.d0-sc2

        do j = 1, totmon
           jjj = j
           if(ptnm(CAindex(j)).eq.'HX2 ') jjj=j+3
           if(ptnm(CAindex(j-1)).eq.'HX2 ') jjj=j+2
           if(ptnm(CAindex(j-2)).eq.'HX2 ') jjj=j+1

           if(ptnm(CAindex(j)).eq.'OX2 ') jjj=j-2
           if(ptnm(CAindex(j+1)).eq.'OX2 ') jjj=j-1

           do l=1,3
             do k =1,5
               A3(l,k) = coor(l,CAindex(jjj-3+k))
             end do
           enddo
           call Find_Center(5,A3,Center1)

           do l=1,3
             do k =1,5
               B3(l,k) = coor2(l,CAindex(jjj-3+k))
             end do
           enddo
           call Find_Center(5,B3,Center2)

           do l=1,3
             do k =1,5
               R3(l,k) = result(l,CAindex(jjj-3+k))
C              write(6,*)"CAindex(jjj-2+k): ",CAindex(jjj-2+k)
C     $                 ,result(l,CAindex(jjj-2+k))
             end do
           enddo
           call Find_Center(5,R3,CenterR)

           do l=1,3
            do k = 1,5
              A3(l,k) = A3(l,k) - Center1(l)
              B3(l,k) = B3(l,k) - Center2(l)
              R3(l,k) = R3(l,k) - CenterR(l)
            end do
           enddo
           call Find_Rotation(5,R3,A3,RotMatrix1(1,1))
           call Find_Rotation(5,R3,B3,RotMatrix2(1,1))

           do i =poipt(j-1)+1, poipt(j)
              if(ptnm(i).ne.'CA  ' .and. ptnm(i).ne.'OX2 ')then
                do l=1, 3
                  rel1(l)=coor(l,i) - Center1(l)
                  rel2(l)=coor2(l,i)- Center2(l)
                end do
                call RotateVector(rel1,RotMatrix1)
                call RotateVector(rel2,RotMatrix2)
                do l=1,3
                  result(l,i)=CenterR(l) + sc1*rel1(l) + sc2*rel2(l)
                end do
              endif
            enddo

          end do

        return
        end

C************************************************************************
C************************************************************************

        subroutine Add_atoms(index,igrid,result,CAindex)
        
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/CCRD.BLOCK'

        integer index,igrid, CAindex(maxmono)
        double precision result(3,maxpt), sc1, sc2
        double precision rel1(3), rel2(3)
        double precision CenterR(3),Center1(3),Center2(3)
        double precision A3(3,3),B3(3,3),R3(3,3)
        double precision RotMatrix1(3,3),RotMatrix2(3,3)

        integer i,j,k,l,jjj

          sc2=(index-1.0)/(igrid-1.0)
          sc1=1.d0-sc2
          
        do j = 1, totmon
           jjj = j
           if(ptnm(CAindex(j)).eq.'HX2 ') jjj=j+2
           if(ptnm(CAindex(j)).eq.'OX2 ') jjj=j-1
           if(ptnm(CAindex(j-1)).eq.'HX2 ') jjj=j+1
           
           do l=1,3
             do k =1,3
               A3(l,k) = coor(l,CAindex(jjj-2+k))
             end do
           enddo  
           call Find_Center(3,A3,Center1)
           
           do l=1,3
             do k =1,3
               B3(l,k) = coor2(l,CAindex(jjj-2+k))
             end do 
           enddo
           call Find_Center(3,B3,Center2)
           
           do l=1,3
             do k =1,3
               R3(l,k) = result(l,CAindex(jjj-2+k))
C              write(6,*)"CAindex(jjj-2+k): ",CAindex(jjj-2+k)
C     $                 ,result(l,CAindex(jjj-2+k))
             end do
           enddo
           call Find_Center(3,R3,CenterR)
           
           do l=1,3
            do k = 1,3
              A3(l,k) = A3(l,k) - Center1(l)
              B3(l,k) = B3(l,k) - Center2(l)
              R3(l,k) = R3(l,k) - CenterR(l)
            end do  
           enddo
           call Find_Rotation(3,R3,A3,RotMatrix1(1,1))
           call Find_Rotation(3,R3,B3,RotMatrix2(1,1))
            
           do i =poipt(j-1)+1, poipt(j)
              if(ptnm(i).ne.'CA  ' .and. ptnm(i).ne.'OX2 ')then
                do l=1, 3
                  rel1(l)=coor(l,i) - Center1(l)
                  rel2(l)=coor2(l,i)- Center2(l)
                end do
                call RotateVector(rel1,RotMatrix1)
                call RotateVector(rel2,RotMatrix2)
                do l=1,3
                  result(l,i)=CenterR(l) + sc1*rel1(l) + sc2*rel2(l)
                end do
              endif
            enddo 
            
          end do

        return
        end

C************************************************************************
C************************************************************************
        subroutine Align(ref,target)
C   ref - is whole atom structure but only CA, CTRM position are correct
C   target - is all atom structure, same as source but rotated/translated to different position
C   the subroutine allignes CAs from ref with those in target
C   gets the rotation and translation of the allignment and applies
C   it to the target as whole

        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
          
          double precision ref(3,*), target(3,*), tmp(3), rms
          double precision refCA(3,maxmono), tarCA(3,maxmono)
	  double precision mmm(maxmono)
          integer nCA,k,i
         
          nCA=0
          do k=1,npt
               if (ptnm(k).eq.'CA  '.or.ptnm(k).eq.'OX2 ' )  then
                  nCA = nCA + 1
                  rms_pick(nCA)=k
                  mmm(nCA) = ptms(k)
                  refCA(1,nCA) = ref(1,k)
                  refCA(2,nCA) = ref(2,k)
                  refCA(3,nCA) = ref(3,k)

                  tarCA(1,nCA) = target(1,k)
                  tarCA(2,nCA) = target(2,k)
                  tarCA(3,nCA) = target(3,k)
               end if
          enddo 
          write(6,*) "Align nCA= ",nCA 
          
          call rmsd_weight(nCA,ref(1,1),target(1,1),rms,.false.,mmm)

          do k = 1,npt
             target(1,k)=target(1,k)-Bcm(1)
             target(2,k)=target(2,k)-Bcm(2)
             target(3,k)=target(3,k)-Bcm(3)
          end do

          do k=1,npt
            tmp(1)=rotat(1,1)*target(1,k)+rotat(1,2)*target(2,k)
     1            +rotat(1,3)*target(3,k)
            tmp(2)=rotat(2,1)*target(1,k)+rotat(2,2)*target(2,k)
     1            +rotat(2,3)*target(3,k)
            tmp(3)=rotat(3,1)*target(1,k)+rotat(3,2)*target(2,k)
     1            +rotat(3,3)*target(3,k)

            target(1,k) = tmp(1) + Acm(1)
            target(2,k) = tmp(2) + Acm(2)
            target(3,k) = tmp(3) + Acm(3)
          end do


        end
