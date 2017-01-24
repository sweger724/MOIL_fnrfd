        subroutine rchain(urcrd,styl,interpolate)
c
c a subrotuine to read initial chain coordinates. Three styles
c are supported: 
c (i)   DYNAmics (binary single precision format for compatability
c       with QUANTA)
c (ii)  PATH (double precision format recommend for computations.
c (iii) INIT (initialize a chain using linear interpolation, urcrd
c               contains the names of normal coordinate files for
c               reactants and products).
c (iv)  INTR (interpolation. Missing structures are interpolated
c             (linearly) between existing structure. Here urcrd IS NOT a
c             coordinate file but rather a set of directions
c             written in a free format (the crd file must be PATH)
c
c             [number of old structures to be read]
c             [indices of the old structures in the new set- 1 4 5.. ]
c             [name of old coordinate file]
c
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        integer urcrd,interpolate
        integer rbin
        logical find
        integer of, stru_per_proc, add_extra
        character*4 styl
        logical skip
c local
        integer i,j,k,l,i1,i2,namel,jprev,jint
        integer u1,u2,countsave
        integer nstru,istru(lgrid)
        integer nofreez(maxpt),skipcount
        double precision div,rms, add(3,maxpt)
        integer istart_eff,kstru_eff,istart_file,j1,iend_eff,n

        character*6 name

        rbin = 1

        norew = .true.
        lpstr = 1

        name  = 'rchain'
        namel = 6

        npt3 = 3*npt
        npt2 = 2*npt

        if (styl.eq.'DYNA') then
C ADD
           if (paral) then
              write(6,*) 'STOP, DYNA OPTION IS NOT PARALLELIZED'
              stop
            endif
         
            i1 = npt
            rewind urcrd
            
            do 1 i=1,igrid
               j = (i-1)*npt
               call rdyncrd(urcrd,i,i1,nofreez,rbin)
               do 1 k=1,npt
                 l = j + k
                 r(1,l)  = coor(1,k)
                 r(2,l)  = coor(2,k)
                 r(3,l)  = coor(3,k)
1           continue

          else if (styl.eq.'PATH') then
             stru_per_proc = (igrid-2)/numprocs
             rewind urcrd
             
C            Read very first structure
             if (first) then
               write(stdo,*)' Reading path format i = ',1
               call rpath_seq(urcrd,1)
               do  k=1,npt
                 do l = 1,3
                    r_initial(l,k) = coor(l,k)
                 end do
               end do 
               call rmsd_weight(npt,coor,r_initial,rms,.true.,ptms)
             endif

             if(paral) then
               if(first) then
                 istart_eff=2
                 if (last) then
                   iend_eff=igrid-1
                 else
                   iend_eff=stru_per_proc+1
                   if (mod(igrid-2,numprocs).ne.0) then
                      iend_eff = iend_eff +1
                   endif
                 endif

C WE ARE GOING TO TAKE OUT INTERPOLATION FROM this parallel code
C Use the serial version if you need to interpolate structures
               do  i=istart_eff,iend_eff
                   write(stdo,*)' Reading path format i = ',i

                   call rpath_seq(urcrd,i)
                   write(stdo,*)' Storing path format i = ',i
                   call rmsd_weight(npt,r_initial,coor,rms,.true.,ptms)
                   j = (i-1)*npt
                   do  k=1,npt
                      l = j + k
                      r(1,l)  = coor(1,k)
                      r(2,l)  = coor(2,k)
                      r(3,l)  = coor(3,k)
                   end do
               end do
c
c reading on master processor and sending to other processors
c
               if(.not.last) then
                 do n=2,numprocs
                   istart_eff=iend_eff+1
                   iend_eff=istart_eff+stru_per_proc-1
C  If not igrid-2/proc then first mod(igrid-2,numprocs) must receive one extra structure
                   if (n.le.mod(igrid-2,numprocs))
     >                           iend_eff = iend_eff +1
C     The last processor must also receive the final structure
                   if (n.eq.numprocs) iend_eff = iend_eff + 1
                        
                   do i=istart_eff,iend_eff
                     write(stdo,*)' Reading path format i = ',i
                     call rpath_seq(urcrd,i)
                     write(stdo,*)' Sending path format i = ',i,
     >                            ' to PROC_ID = ',(n-1)
                    call rmsd_weight(npt,r_initial,coor,rms,.true.,ptms)
                     call send_struc(i,(n-1))
                   end do
                 end do
                 
               endif  ! if not last

C   store final structure in the first processor
               do  k=1,npt
                 r_final(1,k)  = coor(1,k)
                 r_final(2,k)  = coor(2,k)
                 r_final(3,k)  = coor(3,k)
               end do
C         (not first)  
          else
c
c    Receiving section (only for paralel code !!!)
c
            if(mod(igrid-2,numprocs).lt.procID) then
              add_extra = mod(igrid-2,numprocs)
            else
              add_extra = procID
            endif
            istart_eff=procID*stru_per_proc + 2 + add_extra
            iend_eff=istart_eff+pseg-1
        

            n=procID
C
            do i=istart_eff,iend_eff
              write(stdo,*)' Receiving path format i = ',i,
     >                     ' at PROC_ID = ',n
              call recv_struc(i,n)
              j = (i-istart_eff+1)*npt
              do  k=1,npt
                    l = j + k
                    r(1,l)  = coor(1,k)
                    r(2,l)  = coor(2,k)
                    r(3,l)  = coor(3,k)
              end do
            enddo
                
C     The last processor must also receive the final structure
            if (last) then
              write(stdo,*)' Receiving path format i = ',igrid,
     >                     ' at PROC_ID = ',n
              call recv_struc(igrid,n)
              do  k=1,npt
                    r_final(1,k)  = coor(1,k)
                    r_final(2,k)  = coor(2,k)
                    r_final(3,k)  = coor(3,k)
              end do 
            end if
c
c     Closing the PATH reading and storing section !!!!
c
           endif
C If serial:
        else
           write (6,*) 'interpolate mode ',interpolate
           if(interpolate.eq.0.or.interpolate.eq.1) then
              do  k=1,npt
                r(1,k)  = r_initial(1,k)
                r(2,k)  = r_initial(2,k)
                r(3,k)  = r_initial(3,k)
              end do
              do  i=2,igrid-1

                 call rpath_seq(urcrd,i)
                 call rmsd_weight(npt,r_initial,coor,rms,.true.,ptms)
                 if (interpolate.eq.0) then
                    j = (i-1)*npt

C       If we are interpolating up then leave a space btwn
C       each structure

                 else if (interpolate.eq.1) then
                    j = (i-1)*2*npt
                 end if

                 do  k=1,npt
                    l = j + k
                    r(1,l)  = coor(1,k)
                    r(2,l)  = coor(2,k)
                    r(3,l)  = coor(3,k)
                 end do
c       if we are adding intermediate structures
c       and we have already read in the first structure
c       then interpolate btw this one and the one before this

                 if(interpolate.eq.1.and.i.ne.1) then
c       this structure,one before , and the one we will interpolate
                    j = (i-1)*2*npt
                    jprev = (i-1)*2*npt-2*npt
                    jint = (i-1)*2*npt-npt
                    do k = 1,npt
                       r(1,jint+k)=(r(1,j+k)+r(1,jprev+k))/2.0d0
                       r(2,jint+k) =(r(2,j+k)+r(2,jprev+k))/2.0d0
                       r(3,jint+k) =(r(3,j+k)+r(3,jprev+k))/2.0d0
                    end do
                 end if
              end do

c interpolate mode is 2 (go down)

           else

              skip=.false.
              countsave = 0
              skipcount = 10000000
              skpno=skpno-1
              do i = 1,igrid

                 if (skipcount.ge.skpno) then
                    write (6,*) 'reading in structure ',i
                    call rpath_seq(urcrd,i)
                    countsave = countsave+1
                    j = (countsave-1)*npt
                    do k = 1,npt
                       l = j + k
                       r(1,l)  = coor(1,k)
                       r(2,l)  = coor(2,k)
                       r(3,l)  = coor(3,k)

                    end do
                    skipcount = 0
                    skip=.true.
                 else
                    call rpath_seq(urcrd,i)
                    write (6,*) 'skipping structure ',i
                    skip=.false.
                    skipcount = skipcount+1
                 end if
              end do

              igrid = countsave
              write (6,*) 'I have removed intermediate structures'
              write (6,*) 'and changed igrid to ',igrid
           end if

C            Read the very last structure
             write(stdo,*)' Reading path format i = ',igrid
             call rpath_seq(urcrd,igrid)
             call rmsd_weight(npt,r_initial,coor,rms,.true.,ptms)
             
             do  k=1,npt
               do l = 1,3
                 r_final(l,k)  = coor(l,k)
               end do
             end do

             if (interpolate.eq.1) then
                 jprev = (2*igrid-3)*npt-npt
                 jint  = (2*igrid-2)*npt-npt
                 do k = 1,npt
                       r(1,jint+k) =(r_final(1,k)+r(1,jprev+k))/2.0d0
                       r(2,jint+k) =(r_final(2,k)+r(2,jprev+k))/2.0d0
                       r(3,jint+k) =(r_final(3,k)+r(3,jprev+k))/2.0d0
                 end do
                 igrid = (igrid-1)*2 +1
                 write (6,*) 'resetting our grid size to ',igrid
             end if

         end if

       write(6,*) 'I finish to read the initial path'
         

        else if (styl.eq.'INIT') then
         stru_per_proc = (igrid-2)/numprocs
        
         if (first) then
          if (urcrd.ne.stdi) rewind urcrd
c getting unit number of reactants (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
            u1 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
          end if
c getting unit number of products (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
            u2 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
          end if

          call getcrd(u1,'CHARM')
          do  k=1,npt
            do l = 1,3
              r_initial(l,k)  = coor(l,k)
            end do
          end do
C       center r_initial to 0,0,0        
          call rmsd_weight(npt,coor,r_initial(1,1),rms,.true.,ptms)

          call getcrd(u2,'CHARM')
          do  k=1,npt
            do l = 1,3
              r_final(l,k)  = coor(l,k)
            end do
          end do
C       center r_final to 0,0,0        
          call rmsd_weight(npt,coor,r_final(1,1),rms,.true.,ptms)

c overlapping products with respect to reactants
      call rmsd_weight(npt,r_initial(1,1),r_final(1,1),rms,.false.,ptms)

          do  k=1,npt
            do l = 1,3
              r(l,k)  = r_initial(l,k)
            end do
          end do 

c calculating step size.
          div = dble(1.d0/(igrid-1))

          do i=1,npt
            add(1,i) = (r_final(1,i)-r_initial(1,i))*div
            add(2,i) = (r_final(2,i)-r_initial(2,i))*div
            add(3,i) = (r_final(3,i)-r_initial(3,i))*div
          end do



        end if ! first

        
         if (paral) then
C           write(6,*) 'STOP, INIT (cini) OPTION IS NOT PARALLELIZED'
C           stop

           if (first) then
C              prepare mine structures
             do i = 2, pseg+1
               write(stdo,*)' Preparing structure i = ',i
               call PrepareStructure(r(1,(i-2)*npt+1),add(1,1))
               do j = 1,npt
                 r(1,(i-1)*npt+j) = coor(1,j)
                 r(2,(i-1)*npt+j) = coor(2,j)
                 r(3,(i-1)*npt+j) = coor(3,j)
               end do
             end do

C          send structures to other processors
             iend_eff = pseg+1
             if(.not.last) then
               do n=2,numprocs
                 istart_eff=iend_eff+1
                 iend_eff=istart_eff+stru_per_proc-1
C  If not igrid-2/numprocs then first mod(igrid-2,numprocs) must receive one extra structure
                 if (n.le.mod(igrid-2,numprocs)) iend_eff = iend_eff +1
C     The last processor must also receive the final structure
                 if (n.eq.numprocs) iend_eff = iend_eff + 1

                 do i=istart_eff,iend_eff
                   write(stdo,*)' Preparing structure i = ',i
                   call PrepareStructure(coor(1,1),add(1,1))
                   if(i.eq.igrid) then
                     do j = 1,npt
                         coor(1,j) = r_final(1,j)
                         coor(2,j) = r_final(2,j)
                         coor(3,j) = r_final(3,j)
                     end do
                   endif
                   write(stdo,*)' Sending path format i = ',i,
     >                          ' to PROC_ID = ',(n-1)
                   call send_struc(i,(n-1))
                 end do
               end do
             endif  ! if not last


           else ! if not first
c
c    Receiving section (only for paralel code !!!)
c
             if(mod(igrid-2,numprocs).lt.procID) then
               add_extra = mod(igrid-2,numprocs)
             else
               add_extra = procID
             endif
             istart_eff=procID*stru_per_proc + 2 + add_extra
             iend_eff=istart_eff+pseg-1

             n=procID
             do i=istart_eff,iend_eff
               write(stdo,*)' Receiving path format i = ',i,
     >                     ' at PROC_ID = ',n
               call recv_struc(i,n)
               j = (i-istart_eff+1)*npt
               do  k=1,npt
                 l = j + k
                 r(1,l)  = coor(1,k)
                 r(2,l)  = coor(2,k)
                 r(3,l)  = coor(3,k)
               end do
             enddo

C     The last processor must also receive the final structure
            if (last) then
              write(stdo,*)' Receiving path format i = ',igrid,
     >                     ' at PROC_ID = ',n
              call recv_struc(igrid,n)
              do  k=1,npt
                    r_final(1,k)  = coor(1,k)
                    r_final(2,k)  = coor(2,k)
                    r_final(3,k)  = coor(3,k)
              end do
            end if

           end if ! if first

        
C IF Serial:
         else
          
c calculate intermediate structures
c
          do i=2,igrid
            k= (i-1)*npt
            call PrepareStructure(r(1,k-npt+1),add(1,1))
            do j = 1,npt
C              write(6,*)"CCC:",coor(1,j)
              r(1,k+j) = coor(1,j)
              r(2,k+j) = coor(2,j)
              r(3,k+j) = coor(3,j)
            end do
          end do

          do j = 1,npt
              r_final(1,j) = coor(1,j)
              r_final(2,j) = coor(2,j)
              r_final(3,j) = coor(3,j)
          end do

        endif

        else if (styl.eq.'INTR') then
           do 7 i=1,lgrid
              istru(i) = 0
 7         continue
         if (urcrd.ne.stdi)  rewind urcrd
         read(urcrd,*,err=999)nstru
         read(urcrd,*,err=999)(istru(i),i=1,nstru)
         if (istru(1).ne.1) then
          level = 1
          call alert(name,namel,'Missing reactants',17,level)
         else if (istru(nstru).ne.igrid) then
          level = 1
          call alert(name,namel,'Missing products',16,level)
         end if
         call rline(name,namel,urcrd)
         if (find('file')) then
          u1 = of()
         else
          level = 1
          call alert(name,namel,'Missing file name',17,level)
         end if
         do 8 i=1,nstru
          call rpath(u1,i)
          j = (istru(i)-1)*npt
          do 8 k=1,npt
           l = j + k
           r(1,l)  = coor(1,k)
           r(2,l)  = coor(2,k)
           r(3,l)  = coor(3,k)
8        continue
         do 11 i=1,nstru-1
          j = (istru(i)-1)*npt  + 1
          k = (istru(i+1)-1)*npt + 1
          l = istru(i)*npt+1
c Remove overlap because of hydrogen (AVI)
c         call rmsd_weight(npt,r(1,j),r(1,k),rms,.false.,ptms)
          do 9 i1=0,npt-1
           r(1,l+i1) = (r(1,j+i1)+r(1,k+i1))*0.5d0
           r(2,l+i1) = (r(2,j+i1)+r(2,k+i1))*0.5d0
           r(3,l+i1) = (r(3,j+i1)+r(3,k+i1))*0.5d0
9         continue
11        continue
         end if
         return
999      continue
         level = 1
         call alert(name,namel,'Error reading nstru data',24,level)
         return
         end


C*******************************************************************
        subroutine PrepareStructure(old,add)

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/PATH.BLOCK'

        integer i
        double precision add(3,*),old(3,*)

C       prepare mine structures
        do i = 1, npt
          coor(1,i)  = old(1,i) + add(1,i)
          coor(2,i)  = old(2,i) + add(2,i)
          coor(3,i)  = old(3,i) + add(3,i)
        end do

        return
        end
