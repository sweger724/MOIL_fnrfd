      subroutine path_2(nstep,temper, constant_tmp,
     1     pointr,ipick,nselec,d0,grdcmx,grdcmy,grdcmz,
     2     grdlx,grdly,grdlz,udata,ucrd,
     3     debug,scalar,sigmav,nlist,cgrd2,
     5     dtopt,clo,npri)    
c    
      implicit none
      double precision tmp1, tmp2, tmp3
      double precision temper, cgrd2
      integer nstep,nselec,nlist,ncalls
      integer ucrd,udata,finished,npri
      logical debug
c     
c     npri   -  print some useful(?) data each NPRI steps
c     nstep  -  number of minimization steps
c     nselec -  number of selected particles, on which the chain constraints
c     are imposed.
c     ucrd   -  unit number of file on which the coordinates (path format)
c     are written.
c     udata  -  write info on the run on unit UDATA
c     nwcrd  -  write coordinates on ucrd each NWCRD steps
c     debug  - if .true. print a LOT of debugging info
c     
      double precision grdcmx(*),grdcmy(*),grdcmz(*)
      double precision grdlx(*),grdly(*),grdlz(*)
      double precision d0(*)
      double precision scalar(*),sigmav(*)
      double precision estred

      integer n, ji, jf, l
c     
c     dmass - double precision mass vector (m(i=1,npt)(j),j=1,3*igrid)
c     divms - double precision 1/mass vector(1/m(i=1,npt)(j),j=1,3*igrid)
c     
c     d0         -  vector of distances between i,i+1 & i,i+2 pairs
c     grdcm[x-z] -  gradient of center of mass constraints
c     grdl[x-z]  -  gradient of infitesimal rotation constraints
c     r - coordinates (rx(i=1,npt),ry(i=1,npt),rz(i=1,npt)(j=1,igrid))
c     dsall- forces      (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c     
      integer pointr(*),ipick(*)
c     
c     pointr - a pointer to the selected particles
c     (which are subject to chain const.
c     


C     
C     next common blocks are used in the energy calculation.
C     list, debugging, coordinate transfer.. etc
C     
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/MASSFAC.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/PATH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      
c     
c     local
      integer i,j,k,ndegf,imin
      integer ncall
      double precision sener,dtopt
      logical constant_tmp
      double precision step,factor,ss(2),tsa1(2*lgrid-2)
      double precision clo,dclo,run_temp, allmass
      
      integer amount

      ncall  = 0
      npt3  = npt*3
      
      ndegf= (pseg)*npt3
        
      dclo=clo/(nstep+1)
      write(6,'(10x,a,i6,1x,i3)') '# of degrees of freedom= ',ndegf,pseg
C  compute energy of the first stucture
      if (first) then
        do i =1, npt
          do l =1,3
             coor(l,i)=r_initial(l,i)
          end do
        end do
        
        if (eenmyes) call enm_lists(r_initial(1,1),r_final(1,1))
        if (eCGyes)  call CGinit()
        
        if (esymyes) call squeeze()
        call nbondm()
        if (esymyes) call syminit()
        call eforce()
        call wener(stdo)

        e0(1) = e_total
      endif


C  compute energy of the last stucture
      if (last) then
        do i =1, npt
          do l =1,3
             coor(l,i)=r_final(l,i)
          end do
        end do
        if (esymyes) call squeeze()
        call nbondm()
        if (esymyes) call syminit()
        call eforce()
        call wener(stdo)

        e0(pseg+2) = e_total
      endif        

      call initPath()

      allmass = 0.d0
      if (MASSWEIGHT.eq.1) then
         do i = 1,npt
            massfac(i) = 1.0d0/dsqrt(ptms(i))
            allmass = allmass + ptms(i)
         end do
         call mass_weight_coordinates(pseg,r(1,npt+1))
         if (first)  then
             call mass_weight_coordinates(1,r_initial)
             call mass_weight_coordinates(1,r_final)
         else if (last) then
            call mass_weight_coordinates(1,r_final)
         end if
      else
         do i = 1,npt
            allmass = allmass + 1.d0
            massfac(i) = 1.0d0
         end do
      end if


        if (paral) then
C  initial and final structures must be comunicated to each processor
        if(first) then
         do i=1,npt
           coor(1,i) = r_initial(1,i)
           coor(2,i) = r_initial(2,i)
           coor(3,i) = r_initial(3,i)
         enddo

         do i=2,numprocs
           call send_struc(1,(i-1))
         end do
        else
          call recv_struc(1,procID)
          do i=1,npt
            r_initial(1,i) = coor(1,i)
            r_initial(2,i) = coor(2,i)
            r_initial(3,i) = coor(3,i)
          enddo
!  initiate the CG force fields (they need coordinates for initialization)
          if (eenmyes) call enm_lists()
          if (eCGyes)  call CGinit()
          
        endif

        if(first) then
         do i=1,npt
           coor(1,i) = r_final(1,i)
           coor(2,i) = r_final(2,i)
           coor(3,i) = r_final(3,i)
         enddo

         do i=2,numprocs
           call send_struc(igrid,(i-1))
         end do
        else
          call recv_struc(igrid,procID)
          do i=1,npt
            r_final(1,i) = coor(1,i)
            r_final(2,i) = coor(2,i)
            r_final(3,i) = coor(3,i)
          enddo
        endif

      endif

C     Here I have to update fixed end structures:      
      call Communicate_Positions(npt,r,.true.)

            smlp=nlist
            ncalls=nlist
c     
c     Start minimization loop
c     
            write(6,*) 'nstep, smlp',nstep,smlp
            estred = 0.01d0

            do  10  imin = 1,nstep,smlp
               
               finished=0
      
               run_temp = temper *(1.d0 -  1.d0 * (imin-1.d0)/nstep)
C               call path_powel(ndegf,smlp,cgrd2,estred,ncalls,
C     1              sener,pointr,ipick,nselec,
C     2              d0,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,
C     3               scalar,sigmav,finished) 

               call path_anneal(ncalls,run_temp,constant_tmp,
     1              pointr,ipick,nselec,d0,grdcmx,grdcmy,grdcmz,
     2              grdlx,grdly,grdlz,udata,ucrd,
     3              debug,scalar,sigmav,
     4              nlist,sener,dtopt,clo,dclo,npri)

C               call path_test(ncalls,run_temp,constant_tmp,
C     1              pointr,ipick,nselec,d0,grdcmx,grdcmy,grdcmz,
C     2              grdlx,grdly,grdlz,udata,ucrd,
C     3              debug,scalar,sigmav,
C     4              nlist,sener,dtopt,clo,dclo,npri) 

               
     
               write(udata,100) imin+ncalls-1, sener
C               write (udata,*) 'chain energy ',sener
C               write (udata,*) 'clo = ',clo
               
               if (paral) then
                  ss(1)=sener
                  if(.not.first) then
                     call Send_Double(ss,1,0,procID)
                  endif
                  if (first) then
                     do n=1,numprocs-1
                        call Recv_Double(ss,1,n,n)
                        sener=sener+ss(1)
                     enddo
                  endif
               endif
               write (udata,'(a,f13.3)') 'TOTAL ACTION: ',
     &                                   dsqrt(sener/(npt3*(igrid-2)))
 100     format(1x,'chain energy after ',i6,' energy calls = ',f16.3)
               write(udata,*)' Stru #    dist(i,i+1)  e '
 
               do i=1,pseg+1
                  tmp1 = d0(i)/dsqrt(allmass)
                  write(udata,101)i,tmp1,e0(i)
 101              format(1x,i6,1x,2(f14.3,1x))
               end do
               
               tmp1 = 0.d0
               do 4 j=2,pseg+1
                  k = (j-1)*npt
c     write (6,*) 'Structure no ',j

                  do 4 i=1,npt
                     tmp1 = tmp1 + dsall(1,i+k)**2
     1                    + dsall(2,i+k)**2 + dsall(3,i+k)**2
 4             continue
 
                  if (paral) then
                     tmp3=tmp1

C     Send all the gradients to main processor to check convergence!
                     if (.not.first) then
                        call Send_Double(tmp1, 1, 0,procID)
                     endif
C     Receive in first
                     if (first) then
                        do n=1,numprocs-1
                           call Recv_Double(tmp1,1,n,n)
                           tmp3=tmp1+tmp3
                        enddo
                        tmp3 = dsqrt(tmp3/dfloat(npt3*(igrid-2)))
                        write (6,*) 'step no ',imin
                        write(udata,102)tmp3


 102                    format(//,1x,' Current gradient ',f10.5,//)
                        if (tmp3.lt.cgrd2) then
                           write(6,*) 'GRadient converged!'
                           finished=1
                        endif
C     Let other processes know that convergence was achieved or not
                        do n=1,numprocs-1
                           call Send_Integer(finished,1,n,imin)
                        enddo
                     endif
                     if(.not.first) then
                        call Recv_Integer(finished,1,0,imin)
                     endif
                  else
                     tmp1 = dsqrt(tmp1/dfloat(ndegf))
                     write (6,*) 'step no ',imin
                     write(udata,102)tmp1

                     if (tmp1.lt.cgrd2) then
                        write(6,*) 'GRadient converged!'
                        finished=1
                     endif
                  endif
                  if (finished.eq.1) then
                     goto 200
                  end if

 10            continue

 200           continue

C                   FINAL WRITINGS
            write(udata,*)' Writing coordinate at step',imin
            call Write2Path(ucrd,"aaa",3,r(1,1),npt,.false.)


C     Sending and writing final distances between structures to FIRST. 
               if (paral) then
                  if (.not.first) then
                     do i=1,pseg+1
                        tsa1(i)= d0(i)/dsqrt(allmass)
                     enddo
                     call Send_Double(tsa1,pseg+1,0,procID)
                  else ! first
                     write(udata,*) 'Processor 0'
                     
                     do i=1,pseg+1
                           tmp1 = d0(i)/dsqrt(allmass)
                           write(udata,1010)i,tmp1
                     enddo
                     do n=1,numprocs-1
                        write(udata,*) 'Processor ',n
                        if (n.lt.mod(igrid-2,numprocs)) then
                          amount = (igrid-2)/numprocs +1 
                        else
                          amount = (igrid-2)/numprocs
                        endif
                      call Recv_Double(tsa1,amount+1,n,n)
                        do i=1,amount+1
                           write(udata,1010)i,tsa1(i)
                        enddo
                     enddo
                  endif ! first
               endif ! paral

 1010          format(1x,i6,1x,f14.5,1x)
               write (6,*) 'last line of path_2 '

            return
            end
