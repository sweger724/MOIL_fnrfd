           program clstrs
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
           integer maxclus,maxlist
           parameter (maxclus=1000,maxlist=1000000)

           integer icluster(maxlist)
           integer iseed(0:maxclus)
           integer ipick(maxpt)
           integer irep(maxclus)
           integer counter(maxclus)

           integer i,j,k,l,imin
           integer junk,nlist,nclus
           integer of,geti,nstru
           integer namel,level,npick,iselect
           integer ulist,ucon,ucrd,uwclu

           real randvec(maxclus)
           real rmscluster(maxlist)
          
           double precision getd
           double precision irms(maxclus)
           double precision average(3,maxpt,maxclus)
           double precision clusters(3,maxpt,maxclus)
           double precision  mmm(maxpt)
           double precision rms,rmsmin,rcut
           character*7 name
           character*1 tmp
           logical find,fopen,kmin,distance
           data ucon,ulist,uwclu/3*99/

           norew = .true.
           iseed(0) =1

           stdi=5
           stdo=6

           totmon=0
           npt=0
           name='cluster'
           namel=7
c  open junk file for rline
c
            jnkf=25
            open(unit=jnkf,status='scratch')
c default parameters
            rcut = 3.d0
            nclus = 500
            kmin = .false.
            distance = .false.
            level = 0
1           continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('conn')) then  
                ucon=of()
c get connectivity 
                call rconn(ucon)
               end if
               if (find ('list')) ulist=of()
               if(find('wclu'))  uwclu=of()
              else 
               if (find('kmin')) then
                  kmin = .true.
                  nclus = geti('#cls',nclus)
               else if (find('dist')) then
                  distance = .true.
                  rcut = getd('rcut',rcut)
               end if
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
               if (find ('action')) goto  5
             end if
             go to 1
5            continue
             if (iselect.eq.0) then
              level = 1
              call alert(name,namel,'no selection',12,level)
             else if (.not.kmin .and. .not.distance) then
              level = 1
              call alert(name,namel,'no clus method',14,level)
             else if (.not. fopen(ulist)) then
              level =1
              call alert(name,namel,'list not open',13,level)
             else if (.not. fopen(ucon)) then
              level =1
              call alert(name,namel,'wcon not open',13,level)
             else if (.not. fopen(uwclu)) then
              level =1
              call alert(name,namel,'wclu not open',13,level)
             end if
             
c check hoe many structure do we have
            rewind ulist
            nstru = 0
            do i=1,maxlist
              read(ulist,7,end=6)tmp
              nstru = nstru + 1
            end do
6           continue

            if (kmin) then
c generate seeds
             call ranlux(randvec,maxclus)
             do i=1,maxclus
              iseed(i) = int((randvec(i)+0.5))*nstru
             end do

             do i=1,maxclus
              rewind ulist
              do j=1,iseed(i)-1
                read(ulist,7)tmp
7               format(a1)
              end do
              call rline(name,namel,ulist)
              ucrd = of()
              call getcrd(ucrd,'pdb')
              close (ucrd)
              do k=1,npt
               do l=1,3
                clusters(l,k,i) = coor(l,k)
               end do
              end do
             end do

c create clusters

c zero the average coordinate

            do i=1,maxclus
             do l=1,npt
              do k=1,3
               average(k,l,i) = 0
              end do
             end do
            end do

            rmsmin = 999.d0
            imin   = -1
            rewind ulist
            do i=1,nlist
             call rline(name,namel,ulist)
             ucrd = of()
             call getcrd(ucrd,'pdb')
             close (ucrd)
             do j=1,maxclus
c Note: coor is moving
              if (iseed(j).ne.i) then
               call rmsd_weight(npt,clusters(1,1,j)
     1          ,coor,rms,.false.,mmm)
               if (rms.lt.rmsmin) then
                imin = iseed(j)
                rmsmin = rms
                counter(j) = counter(j) + 1
               end if
              end if
             end do
             icluster(i)=imin
             rmscluster(i) = rmsmin
c average is done after optimal overlap
             do j=1,3
              do k=1,npt
               average(j,k,imin) = average(j,k,imin) + 
     1         + mmm(k)*coor(j,k)
              end do
             end do
            end do
c compute average structures for current cluster
            do i=1,maxclus
             do j=1,npt
              do k=1,3
               average(k,j,i) = average(k,j,i)/counter(i)
              end do
             end do
            end do
c find the structure closest to the average
            do j=1,maxclus
             rewind ulist
             rmsmin = 999.
             imin   = -1
             do i=1,nlist
              call rline(name,namel,ulist)                 
              ucrd = of()
              call getcrd(ucrd,'pdb')
              close (ucrd)
              if (iseed(j).ne.i) then
               call rmsd_weight(npt,average(1,1,j),
     1           coor,rms,.false.,mmm)
               if (rms.lt.rmsmin) then
                rmsmin = rms
                irep(j) = i
                irms(j) = rmsmin
               end if
              end if
             end do
            end do

            do j =1,maxclus
             rewind ulist
             do i=1,irep(j)-1
              read(ulist,7)tmp
             end do
             call rline(name,namel,ulist)
             ucrd = of()
             call getcrd(ucrd,'pdb')
             close (ucrd)
             do k=1,npt
              do l=1,3
               clusters(l,k,j) = coor(l,k)
              end do
             end do
            end do
c write some output
            do j=1,maxclus
              write(*,*)irep(j)
            end do
c stop for now. Still need to assess distnce between clusters
c and cluster sizes
            end if
            stop
         end

