      subroutine dyna_mfep(proc,temp,step,step2,pmlst,tmlst
     1,nmlst,nmom
     1,nens ,mxstps,nsvel,pointr,nselec,dmass,divms,grdcmx,grdcmy 
     1     ,grdcmz,grdlx,grdly,grdlz,fact1,fact2,udata,urcrd 
     2     ,uwfpt,uwfpp,uwcrd,uwpep,uene,nene,nwcrd,ntest,ptemp,nequ 
     2     ,newv,irand,debug,scalar,sigma,nlist,massw,shakm,shakb 
     3     ,shakl,nwpep,prp,prn,plp,pln,pmp,pmn,unitgrad,orth,tfp 
     3     ,udot,udot1,udst,umom,ucrb,umlst,upfr,upvl
     4     ,num_mlsts,gamma,smooth,nwav,frep,fixend,repul) 

c     deb-- for mshk of TIP3--shakm and tolmshk added------
c     deb tfac was changed to step2 as tfac was not being used.
C     TF: tolmshk eliminated in favor of tolcons, for compatibility
C     with COMMON/MSHAKE.BLOCK. 8-Feb-2005.
      
c     temp   -  assigned temperature
c     step   -  time integration step (in AMKA units)
c     step2  -  step*step/2.d0
c     pmlst  -  previous milestone number
c     tmlst  -  current milestone number
c     nmlst  -  next milestone number
c     nmom   -  number of steps between distance calculations
c     nens   -  number of configs in urcrd
c     mxstps -  max number of integration steps, not counting
C     equilibration steps. Set to -1 for no max.
c     nsvel  -  number of steps between velocity scaling
c     pointr -  pointer to particles selected for constraint
c     nselec -  number of selected particles, on which the chain
C     constraints are imposed.
c     dmass  -  double precision mass vector (m(i=1,npt)
c     divms  -  double precision 1/mass vector(1/m(i=1,npt)
c     grdcm[x-z] -  gradient of center of mass constraints
c     grdl[x-z]  -  gradient of infitesimal rotation constraints
c     fact[1-2]  - constant vectors useful for integration (see fp.f)
c     udata  -  write diagnostic info to this unit
c     urcrd  -  initial configs path file unit num
c     uwcrd  -  trajectory (path) file unit num
c     uwpep  -  peptide trajectory (path) file unit num
c     upfrc  -  file unit for projected force
c     upvl   -  file unit for projected velocity
c     uene   -  energy file unit num
c     nwcrd  -  write coordinates on uwcrd each nwcrd steps
c     ntest  -  test the constraints each ntest steps
c     ptemp  -  temperature at which to pre-equilibrate
c     nequ   -  number of pre-equilibration steps
c     newv   -  assign new velocities each NEWV steps
c     irand  -  random number generator seed
c     debug  -  if .true. print a LOT of debugging info
c     scalar / sigma - constraint arrays (see crbm & comc)
c     nlist  -  re-build non-bonded list every nlist steps
c     massw  -  boolean: use mass-weighting
c     shakm, shakb, shakl: SHAKE booleans. Use MSHAKE, but not SHAKE.
c     pmp, pmn, plp, pln, prp prn: see fp.f

c     LOCAL
c     coor - coordinates (coor(i,j) i=x,y,z j=1,...,npt )
c     velo - velocities (velo(i,j) i=x,y,z j=1,...,npt)
c     dpot -  forces  (dpot(i,j) i=x,y,z j=1,..,npt)

      implicit none
      integer nomlst
      parameter(nomlst=0)
      integer umlst,upfr,upvl,num_mlsts
      double precision temp,old_tmpr,tmpr,step,step2,ptemp
      integer nmom,nsvel,nselec,iens,nens,mxstps
      integer pmlst,tmlst,nmlst,atmlst,nwpep,natm,patm,nene
      integer udata,urcrd,uwcrd,uwpep,uwfpt,uene,nwcrd,ntest
      integer nequ,newv,irand,nlist,uwfpp,umom,udot,udot1,udst
      logical tfp,orth,debug,massw,rvrs,change
      double precision dmass(*),divms(*)
      double precision grdcmx(3,*),grdcmy(3,*),grdcmz(3,*)
      double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
      double precision fact1(*),fact2(*)
      double precision scalar(*),sigma(*)
      double precision cm(3)
      double precision nranf, nsysf, projf
      integer pointr(*), l,nwav
      double precision factor
      data factor/1.98768d-3/
      logical ucrb,fixend,repul
     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/CCRD.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK' 
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/OVERLAP.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'mpif.h'



      character*20 name
      integer namel,level

c     local
      integer i,j,k,m,npt3,smlp,idyna,istp,sflag,ip
      integer ii, ndegf, iat1, iat2,nselec3
      logical shakm,shakl,shakb
      double precision tmp
      double precision hot_umbr,vfac,sigmav(7),gtemp
      double precision pi
      double precision tmp_rand(6)
      double precision dotmp
      double precision pepcrd(3,maxpt),coor_2(3,maxpt)
      
      
c     plp, pln (prp, prn (pmp, pmn)) : left (right (mid)) fp plane
c     point and normal
      double precision plp(3,maxpt), pln(3,maxpt)
      double precision prp(3,maxpt), prn(3,maxpt)
      double precision pmp(3,maxpt), pmn(3,maxpt),pmp_old(3,maxpt)
      double precision unitgrad(3,maxpt),dpmp(3,maxpt)
      double precision rms,dst_array(maxpt),ll,ndpmp
      double precision sener

      
c     fpt: in units of time step
      double precision fpt, frac
      
c     for file read
      double precision e0
      
c     track configurations
      double precision dst, coorstart(3,maxpt),dsmax,dst2,dst3
      double precision dummy(3,maxpt), tmp_mass(maxpt),dummy2(3,maxpt)
      double precision unitgradnxt(3,maxpt),avg_oldnxt(3,nselec)
      integer npick3

c     functions
      double precision vscalar, dist2_fp,dist2_pp

      double precision dp(6)

c     averge coordinate 
      double precision avg_pep(3,nselec),avg_old(3,nselec)
      double precision av(3,nselec),flc(3,nselec)
      integer ierr
      integer status(MPI_STATUS_SIZE)
      integer proc
      integer p1  ,itr,itrall,frep,ifrep,iupdate
c     string constants 
      double precision gamma,smooth,dstmin
       double precision xx(3,nselec*200),pmpt(3,maxpt),xx1(3,nselec*200)
      double precision ptx(3,maxpt)
c     for trace
      double precision norm1, norm2,norm_dpmp
      double precision dpold(3*maxpt)
      double precision  lambda, rho ,dmass_tot
      double precision dotl,dotr,dotm ,d3

      name = 'dyna_mfep'
      namel = 7
      frac = 0
       itr=0 
       itrall=0 
       ifrep=0 
       iupdate=0 
       norm_dpmp=0.d0 
c     ------------------------------------------------------

      if ( ucrb .and. torscstr ) then
         level = 1
         call alert(name, namel,'ucrb w/ tors constraints',40,level)
      end if


      if ( (mxstps .lt. 0) .and. (orth .or. (.not. tfp)) ) then
         level = 1
         call alert(name, namel,'perpetual simulation',40,level)
      end if


      
      if ( orth ) then
         if ( (.not. ucrb) .and. (.not. torscstr) .and. (.not. rgcst)
     1        ) then
            level = 0
            call alert(name, namel, 'orth without constraint',40,level)
         end if
         if ( shakb .or. shakl ) then
            level = 1
            call alert(name, namel,'Shake with orth',40,level)
         end if
      else
         if ( (.not. tfp) .and. (nens .gt. 1) ) then
            level = 1
            call alert(name, namel, '> 1 free sim traj',40,level)
         end if
      end if




c     ------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      npt3  = npt*3
      nselec3  = nselec*3
      smlp  = 1

c     deb---------Substract shake degrees of freedom-------------

      if (shakm) then
         ndegf = 3*inofrz - nshakm 
      else
         ndegf = 3*inofrz
      endif


c     deb--------------------------------------------------------
c     Data for velocities constraints.
c     sigmav for the velocities is zero (in contrast to the
c     coordinates)
      
      do 1 i=1,7
         sigmav(i)=0.d0
 1    continue
      
      if (nlist .ne. 0) smlp = nlist
      if (mxstps .gt. 0) mxstps = mxstps + nequ

         read(urcrd) e0,((coorstart(k,i),i=1,npt),k=1,3) 
       
 
c      calculate the 1/total mass for the picked particles 
            dmass_tot=0.d0 
                  do j = 1,nselec
                     i = pointr(j)
                     dmass_tot=ptms(i)+dmass_tot
                end do
              
        
        
      
c     BEGIN ENSEMBLE LOOP
      do 7900 iens = 1,nens


         atmlst = tmlst
         dotm =  0.d0

 
         call vdcopy ( coorstart, coor, npt3 )


c       initializing the average coordinate for each milestone 
       do i=1,nselec
        do k=1,3
          avg_pep(k,i) = 0.0d0 
          pepcrd(k,i) = 0.0d0
          pmp_old(k,i)=pmp(k,i)
        enddo          
       enddo 


c        BEGIN OUTER INTEGRATION LOOP
c        initialize integration loop
c        do internal update each inbfrq for the non-bonded AND the images
c        (if applicable). A single non-bonded list is generated according
c        to the first structure.
         idyna = 1

 110     if ( nlist .ne. 0 ) then
            if (esymyes) call squeeze()
            call nbondm()
            if (esymyes) call syminit()
         end if


 
c        write(6,*) 'I am here 1'
        call eforce()
c        write(6,*) 'I am here 1 after'


c        BEGIN INNER INTEGRATION LOOP
         do 100 istp = idyna,idyna+smlp-1


c           check for termination by maximum steps
            if ((mxstps .gt. -1) .and. (istp .gt. mxstps)) then
               write(*,*) ' DONE: istp, mxstps:', istp, mxstps
               go to 7800
            end if

          if ( orth ) write(upfr,*) projf(dpot,unitgrad,nselec,pointr)

c           set temperature
            if ( istp .le. nequ ) then
               gtemp = ptemp
            else
               gtemp = temp
            end if


            dotmp = dotm
c            end if

c           calculate the average conformation of the chain at the plane 
               if(istp.ge.nequ.and.ifrep.eq.0) then 
               itr=itr+1 

               do k = 1,3
                  do j = 1,nselec
                     i = pointr(j)
                  avg_pep(k,j)=avg_pep(k,j) +coor(k,i)


                end do
               end do
             
    
                       if(itr.eq.nwav) then
                       itrall=itrall+1

                        do k= 1,3 
                        do j=1,nselec
                         
                         avg_pep(k,j)= avg_pep(k,j) /nwav
                         pepcrd(k,j)=pepcrd(k,j)+avg_pep(k,j)
                         pmpt(k,j)=pepcrd(k,j)/itrall
                       
                        enddo 
                        enddo 

 
                         
c                        if(mod(itrall,50).eq.0) then 
c                        itrall = 0 
c                        do k= 1,3 
c                        do j=1,nselec
                         
c                         pepcrd(k,j)=0 
c                       
c                        enddo 
c                        enddo 
c                        endif 
                       
                        
c                 do j=1,nselec
c                    i= pointr(j)
c                 write(80+procID,80) istp,itr,i,j,pepcrd(1,j)
c     &,pepcrd(2,j),pepcrd(3,j), coor(1,i),coor(2,i),coor(3,i),
c     $pepcrd(1,j)-coor(1,i)
c                 enddo 
 
            
                        endif 
          
             endif 

             

80        format(4i7,12f12.5) 
             
123       format(i6,3f12.5) 


           if(procID.gt.0.and.procID.lt.proc-1.and.istp.gt.nequ.and.
     >     ifrep.eq.0.and.itr.eq.nwav) then 
             ifrep=1
            itr=1

            do k = 1,3
              do j = 1,nselec
                  avg_pep(k,j)=0. 

                end do
               end do

           do i=1,nselec
           do j=1,3 
           dpmp(j,i)=(pmpt(j,i)-pmp(j,i))/frep

           enddo 
           enddo 
            endif 
            
           if(procID.gt.0.and.procID.lt.proc-1.and.istp.gt.nequ.and.
     >        ifrep.eq.1) then 
           iupdate=iupdate +1 

c           step=step/50. 
c           step2=step*step/2.d0
           do i=1,nselec
           do j=1,3 
           pmp(j,i)=pmp(j,i) + dpmp(j,i)*gamma    
           enddo 
           enddo 
           endif 
                 
           if(iupdate.eq.frep) then 
           ifrep=0
           iupdate=0 

c           step=step*50. 
c           step2=step*step/2.d0

           endif  


124      format(9f12.5)

134   format(a2,i4,9f13.8)

c    now we are smoothing the string 
c   ***************************************



       if(procID . eq. (proc-1) ) then

      do i=1,nselec
      do k=1,3
      prp(k,i)=pmp(k,i)
      enddo
      enddo

         Call MPI_Recv(plp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                        proc-2,1,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.(proc-2)) then
          Call MPI_Send(pmp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                         (proc-1),1,MPI_COMM_WORLD, ierr )
       endif


c     first milestone does not have prev 
       if(procID . eq. 0 ) then

      do i=1,nselec
      do k=1,3
      plp(k,i)=pmp(k,i)
      enddo
      enddo

         Call MPI_Recv(prp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                        1,23,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.1) then
          Call MPI_Send(pmp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                         0,23,MPI_COMM_WORLD, ierr )
       endif



       do i=1,proc-2
       if(procID.eq.i) then
         Call MPI_Recv(plp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                        i-1,13,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i-1)) then
          Call MPI_Send(pmp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                         i,13,MPI_COMM_WORLD, ierr )
       endif


        if(procID.eq.i) then
         Call MPI_Recv(prp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                        i+1,24,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i+1)) then
          Call MPI_Send(pmp, 3*nselec, MPI_DOUBLE_PRECISION,
     >                       i,24,MPI_COMM_WORLD, ierr )
       endif
       enddo
       
      if(procID.gt.0.and.procID.lt.proc-1) then 
      if(ifrep.ne.0.and.istp.gt.nequ) then 
      do i=1,nselec 
      do k=1,3
      pmp(k,i)=(1-smooth)*pmp(k,i)+smooth*0.5d0*(plp(k,i)+ prp(k,i))
      enddo 
      enddo 
      endif 
       endif 


c   *****************************************
c     Now reparametrize the string to satisfy equal distance between points
         if(procID.eq.0) then
         do ip=1,nselec
          do j=1,3  
         xx(j,ip) = pmp(j,ip) 
         enddo      
         enddo            
         endif


          do 111 i=1,proc-1

          if(procID.eq.i) then
          Call MPI_Send(pmp, nselec*3, MPI_DOUBLE_PRECISION,
     >                         0,5,MPI_COMM_WORLD, ierr )
         endif
         if(procID.eq.0) then
         Call MPI_Recv(ptx, nselec*3, MPI_DOUBLE_PRECISION,
     >                   i,5,MPI_COMM_WORLD, status,ierr )

         k = i*nselec
         do ip=1,nselec
          do j=1,3
         xx(j,ip+k) = ptx(j,ip)
         enddo
        enddo
         endif


111     continue


c      reparametrization of the string by linear interpolation  
       if(procID.eq.0.) then 
  
        if(ifrep.ne.0.and.istp.gt.nequ)  then 
        call  reparam(xx,nselec,proc,xx1,ll,istp)


c        if(repul)  call repulsion(xx1,nselec,proc,rho,lambda,sener)


         
        else
         do ip=1,nselec*proc
         do j=1,3 
         xx1(j,ip)=xx(j,ip) 
       enddo 
         enddo 
         endif 


 
         endif
   
133    format(2i8,6f12.5) 


           do 112 i=1,proc-1
          if(procID.eq.0) then
           k = i*nselec
           do ip=1,nselec
            do j=1,3
           ptx(j,ip)=xx1(j,ip+k) 
            
      
            enddo 
            enddo 

          Call MPI_Send(ptx, nselec3, MPI_DOUBLE_PRECISION,
     >                         i,666,MPI_COMM_WORLD,ierr )
          elseif (procID.eq.i) then


          Call MPI_Recv(pmp, nselec3, MPI_DOUBLE_PRECISION,
     >                   0,666,MPI_COMM_WORLD, status,ierr )

         endif


 112     continue


C ******************************************
c End of REPARAMETRIZATION 
c

c   ***************************************
c    now move the string to its new position 

       if(procID . eq. (proc-1) ) then

      do i=1,nselec
      do k=1,3
      prp(k,i)=pmp(k,i)
      enddo
      enddo

         Call MPI_Recv(plp, nselec3, MPI_DOUBLE_PRECISION,
     >                        proc-2,1,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.(proc-2)) then
          Call MPI_Send(pmp, nselec3, MPI_DOUBLE_PRECISION,
     >                         (proc-1),1,MPI_COMM_WORLD, ierr )
       endif


c     first milestone does not have prev 
       if(procID . eq. 0 ) then

      do i=1,nselec
      do k=1,3
      plp(k,i)=pmp(k,i)
      enddo
      enddo

         Call MPI_Recv(prp,nselec3, MPI_DOUBLE_PRECISION,
     >                        1,23,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.1) then
          Call MPI_Send(pmp,nselec3, MPI_DOUBLE_PRECISION,
     >                         0,23,MPI_COMM_WORLD, ierr )
       endif



       do i=1,proc-2
       if(procID.eq.i) then
         Call MPI_Recv(plp, nselec3, MPI_DOUBLE_PRECISION,
     >                        i-1,13,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.(i-1)) then
          Call MPI_Send(pmp, nselec3, MPI_DOUBLE_PRECISION,
     >                         i,13,MPI_COMM_WORLD, ierr )
       endif


        if(procID.eq.i) then
         Call MPI_Recv(prp, nselec3, MPI_DOUBLE_PRECISION,
     >                        i+1,24,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i+1)) then
          Call MPI_Send(pmp, nselec3, MPI_DOUBLE_PRECISION,
     >                       i,24,MPI_COMM_WORLD, ierr )
       endif
       enddo
 


c      if(ifrep.ne.0.and.nequ.lt.istp) then 
 
      call mvplane(plp,prp,pmn,nselec)
      call normalize ( pmn, nselec3, unitgrad, norm1 )


c     print out some useful information
c       if ((istp .eq. 1) .or. (istp / nwcrd * nwcrd) .eq. istp) then

       do i=1,proc-2
       if(procID.eq.i) then
         Call MPI_Recv(pln, nselec3, MPI_DOUBLE_PRECISION,
     >                        i-1,13,MPI_COMM_WORLD, status,ierr )

      endif

       if(procID.eq.(i-1)) then
          Call MPI_Send(pmn, nselec3, MPI_DOUBLE_PRECISION,
     >                         i,13,MPI_COMM_WORLD, ierr )
       endif


        if(procID.eq.i) then
         Call MPI_Recv(prn, nselec3, MPI_DOUBLE_PRECISION,
     >                        i+1,24,MPI_COMM_WORLD, status,ierr )
      endif

       if(procID.eq.(i+1)) then
          Call MPI_Send(pmn, nselec3, MPI_DOUBLE_PRECISION,
     >                       i,24,MPI_COMM_WORLD, ierr )
       endif
       enddo

c           dot products
           
            call vecmin_fp ( coor, pmp, dummy, nselec, pointr )
            dotm = vscalar ( pmn, dummy, nselec3 )
               call vecmin_fp ( coor, plp, dummy, nselec, pointr )
               dotl = vscalar ( pln, dummy, nselec3 )
               call vecmin_fp ( coor, prp, dummy, nselec, pointr )
               dotr = vscalar ( prn, dummy, nselec3 )



         if ( procID.gt.0.and.procID.lt.proc-1) then 
c        
c          do i=2,proc-2
c          if(procID.eq.i) then
c         Call MPI_Recv(d3, 1, MPI_DOUBLE_PRECISION,
c     >                        i-1,13,MPI_COMM_WORLD, status,ierr )
c          endif 
     

c          if(procID.eq.(i-1)) then
c                 Call MPI_Send(dotr, 1, MPI_DOUBLE_PRECISION,
c     >                         i,13,MPI_COMM_WORLD, ierr )
c          endif        
            
 
c      -----------------------------------------------
c      If conformations are closing to a neighbor 
c      inverse the velocity 
c      ---------------------------------------------

         if (dotr.gt.1..or.dotl.lt.1) then
                do k = 1,3
                  do j = 1,nselec
                     i = pointr(j)
                  velo(k,i)= - velo(k,i)


                end do
               end do

c        if (dotr.gt.0.or.dotl.lt.0) then 
c            do k=1,npt 
c             do j=1,3 
c             coor(j,k)=coorstart(j,k) 
c            enddo 
c           enddo 
         endif  
       
      endif


               write(stdo,1000)istp,dotl,dotm,dotr
1000           format(1x,'ISTP = ',i6,/,1x,
     1     ' distances from  planes l m r ',3(e10.3,1x))

        

c              if(mod(istp,nwcrd).eq.0) then
c                  do j = 1,nselec     
c                  write(90+procID,1231) (plp(k,j),k=1,3)
c     >,(pmp(k,j),k=1,3),(prp(k,j),k=1,3),(pmn(k,j),k=1,3)
c                  end do
c              endif 
c          endif 
c      calculate the angle between plane slopes

c       if(ifrep.ne.0.and.nequ.lt.istp) then 


c       do  i=1,proc-1
c          if(procID.eq.i-1) then
c          Call MPI_Send(pmn, nselec3, MPI_DOUBLE_PRECISION,
c     >                         i,666,MPI_COMM_WORLD,ierr )
c          elseif (procID.eq.i) then
c          Call MPI_Recv(dummy2, nselec3, MPI_DOUBLE_PRECISION,
c     >                   i-1,666,MPI_COMM_WORLD, status,ierr )

c     calculate the dot product 
c       dst=0.d0 
c        do ip=1,nselec
c        do j=1,3
c        dst=dst + pmn(j,ip)*dummy2(j,ip) 
c        enddo 
c        enddo     


        

c        endif 
c       enddo 
     
c      endif 

c      calculate the norm of dpmp 
              if(mod(istp,nwcrd).eq.0) then
           norm_dpmp=0.d0 
           do i=1,nselec
           do j=1,3 
           norm_dpmp=norm_dpmp+(pmp_old(j,i)-pmp(j,i))**2.
           enddo 
           enddo 

           norm_dpmp=norm_dpmp           

c                  do j = 1,nselec     
c                  write(90+procID,1231) (pmp(k,j),k=1,3)
c     >,(pmp_old(k,j),k=1,3)
c                  end do
c              endif 
          endif 
1231    format(12(f12.5))
     
           


c     get center of mass constraints and orthonormalize all
c      constraints. (dpold and velo are dummies)
         call comc_fxd(pmp,pmn,divms,npt,pointr,nselec,grdcmx,
     1     grdcmy,grdcmz,grdlx,grdly,grdlz,dpold(1),
     2     dpold(1+npt),dummy,sigma,ucrb,debug)



          dst2=0.d0
          dst=0.d0  

          if(mod(istp,nwcrd).eq.0) then 
            
                  do k = 1,3
                  do j = 1,nselec
                     dummy(k,j) = pmp(k,j)
                  end do
               end do

                  do k = 1,3
                  do j = 1,nselec
                     i=pointr(j)
                     dummy2(k,j)=coorstart(k,i)
                  end do
               end do

c                 do j=1,npt
c                   tmp_mass(j)=ptms(j)
c                    ptms(j)=0.0d0
c                 enddo                      
                 
c                  do k = 1, nselec
c                    ptms(k) = tmp_mass(pointr(k))
c                  end do
            

c                  call rmsd_weight(nselec,dummy,dummy2,rms,.false.)

c                  do k = 1,npt
c                    ptms(k) = tmp_mass(k)
c                  end do

c      calculate the distance traveled by pmp   
                dst3=0.d0 
                 do k=1,3
                 do j=1,nselec
                  
                  dst3=  dst3+ (dummy2(k,j)-dummy(k,j))**2
                  enddo 
                  enddo 
                dst3=sqrt(dst3)/nselec   
                                 


c      calculate the distance between i i-1  
                dst=0.d0 
                 do k=1,3
                 do j=1,nselec
                  dst=  dst+ (plp(k,j)-pmp(k,j))**2
                  enddo 
                  enddo 
                dst=sqrt(dst)/nselec   
c      calculate the distance between i-1 i+1  
                dst2=0.d0 
                 do k=1,3
                 do  j=1,nselec 
                 dst2=  dst2+ (plp(k,j)-prp(k,j))**2
                  enddo 
                  enddo 
                dst2=sqrt(dst2)/nselec   


c                  write(*,*) 'rot matrix wrt initial point:'
c                  write(stdo,*) ''
c                  write(stdo,*)' rotation matrix '
c                  do k = 1,3
c                    write(stdo,3001)(rotat(k,j),j=1,3)
c3001                format(1x,3(f10.4,1x))
c                  end do
      
c              if(mod(istp,nwcrd).eq.0) then
c                  do j = 1,nselec     
c                  write(80+procID,1231) (plp(k,j),k=1,3)
c     >,(pmp(k,j),k=1,3),(prp(k,j),k=1,3),(pmn(k,j),k=1,3)
c                  end do
c             endif 
c1231    format(12f12.5) 
c    ______________________________
          
           
         if(procID.gt.0) then 
         write(stdo,122) istp,dst3,dst,dst2
         endif 
c   add up all dst3 and write the result to pth_0001.log file for further measure of convergence 

           ll=0.d0  
           do i=1,proc-1  

          if(procID.eq.i) then        
          Call MPI_Send(dst3, 1, MPI_DOUBLE_PRECISION,
     >                         0,13,MPI_COMM_WORLD, ierr )
            endif
           if(procID.eq.0) then 
         Call MPI_Recv(dst, 1, MPI_DOUBLE_PRECISION,
     >                        i,13,MPI_COMM_WORLD, status,ierr )
           endif      
           ll = ll+dst

          enddo 

c    add up all norm_dpmp and write the result to pth_0001.log file for measure of convergence     

           ndpmp=0.d0
           do i=1,proc-1

          if(procID.eq.i) then        
          Call MPI_Send(norm_dpmp, 1, MPI_DOUBLE_PRECISION,
     >                         0,13,MPI_COMM_WORLD, ierr )
            endif
           if(procID.eq.0) then 
         Call MPI_Recv(dst, 1, MPI_DOUBLE_PRECISION,
     >                        i,13,MPI_COMM_WORLD, status,ierr )
           endif      
           ndpmp = ndpmp+dst

          enddo
         if(procID.eq.0) then 
         write(stdo,122) istp,ll,ndpmp
         endif




         endif 

122     format(i14,8x,f15.10,3f20.16)  


c           write peptide traj
c           if ( ( nwpep .eq. 1 ) .or. ( ( nwpep .gt. 0 ) .and. ( (
c     1           istp / nwpep * nwpep ) .eq. ( istp - 1 ) ) .and. (
c     2           istp .gt. nequ ) ) ) then
c               do k = 1,3
c                  do j = 1,nselec
c                     i = pointr(j)
c                     pepcrd(k,j) = coor(k,i)
           
                      
c                  end do
c               end do
c               write(uwpep) e_total, ((pepcrd(k,i),i=1,nselec),k=1,3)

c            end if



c           write the FIRST trajectory to file. So when orth
c           is .true. this if block writes the entire milestone
c           ensemble.
            if ( iens .eq. 1 ) then

            if (mod(istp,nwcrd).eq.0) then

                 
                  do k = 1,3
                  do j = 1,nselec
                     pmp_old(k,j)=pmp(k,j)
                  end do
                  end do 



                  do i=1,npt 
                    do j=1,3 
                     dummy2(j,i)=coor(j,i)
                    enddo 
                  enddo 

     

                  do k = 1,3
                  do j = 1,nselec
                     i=pointr(j)
                     dummy2(k,i)=pmp(k,j)
c                    dummy2(k,i)=avg_pep(k,j)/itr 
                  end do
               end do

                  write(uwpep) e_total,(dummy2(1,i),i=1,npt),
     1                 (dummy2(2,i),i=1,npt),(dummy2(3,i),i=1,npt)
                endif 

                


              if (mod(istp,nwcrd).eq.0) then
                 
                  write(uwcrd) e_total,(coor(1,i),i=1,npt),
     1                 (coor(2,i),i=1,npt),(coor(3,i),i=1,npt)
             endif 
c                 write distance travelled

         if ( istp.eq.1) write(umom,1221) 

1221        format(1x,"Distance traveled from init conf :")
  
          if ( umom .gt. 0 .and. ((istp/nmom * nmom).eq.istp) ) then
                     dst = sqrt ( dist2_fp ( coorstart, coor, nselec,
     1                    pointr ) )
c        write the reaction progress 

c           do i=1,proc-2  

c          if(procID.eq.(i+1)) then        
c          Call MPI_Send(coor, npt3, MPI_DOUBLE_PRECISION,
c     >                         i,13,MPI_COMM_WORLD, ierr )
 
 
c           endif
c           if(procID.eq.i) then 
c         Call MPI_Recv(coor_2, npt3, MPI_DOUBLE_PRECISION,
c     >                        i+1,13,MPI_COMM_WORLD, status,ierr )
c           endif      
c          enddo 

c
c                  dst2=0.d0     
c                  do k = 1,3
c                  do j = 1,nselec
c                     i=pointr(j)
c                     dst2=pmn(k,j)*(coor_2(k,i)-coor(k,i))+dst2
c                  end do
c               end do
                  
c             if(dst2.lt.0.and.istp.gt.nequ) then 
c             repul=.true. 
c             else 
c             repul=.false. 
c             endif   

                 write(umom,1234) istp,dst/nselec
                  end if
        
1234         format(i8,2x,f12.7) 

c                 verify non-drift of center of mass
c                  call com_fp ( coor, dmass, cm, nselec, pointr )
c                  write(*,1002)  cm(1), cm(2), cm(3)
c1002              format(1x,' com ',2(f10.4))



c                 verify rotational non-drift
c                  do k = 1, npt
c                    do l = 1,3
c                      dummy(l,k) = coor(l,k)
c                    end do
c                    tmp_mass(k) = ptms(k)
c                    ptms(k) = 0.d0
c                  end do
                  
               end if


c     ---------------------------------------------------
c     If the flow of control reaches this point, then the
c     trajectory has not yet terminated. So proceed to
c     the Verlet step.
c     ----------------------------------------------------


c           sample the velocities at the right temperature
            if ( (istp.eq.1) .or. ((newv.ne.0).and.(istp/newv*newv.eq.
     1           istp)) ) then
               call velinit ( gtemp,1,0 )
            end if

c        do i = 1,npt
c        write(6,*) i,velo(1,i),velo(2,i),velo(3,i)
c        enddo

c           correct newly assigned velos to satisfy constraints
            call crbm ( velo,scalar,sigmav,grdcmx,grdcmy
     1           ,grdcmz,grdlx,grdly,grdlz,divms,pmn,npt,nselec
     1           ,pointr,istp,ntest,ucrb,debug,udata )

c        write(6,*) 'VELO'
c        do i = 1,npt
c        write(6,*) i,velo(1,i),velo(2,i),velo(3,i)
c        enddo


       if(istp.eq.1)   then 
        write(stdo,1220) 
1220     format(1x,"        Time              init mlst      (i-1)-i   
     & (i-1)-(i+1)")
        write(stdo,*)"      -----------       
     & ------------   ----------  -----------" 
       endif 


c     Verlet stage 1
c     calculate a "free" (unconstrained) step
c     Make a step for unfrozen particles (velocity is
c     used to store the step) dividing finally by dt
c     will recover the magnitude required for velocity
c     calculations
            do 10 ii = 1,inofrz
               i = nofreez(ii)
               do 1010 k = 1,3
                  velo(k,i) = velo(k,i)*step - fact1(i)*dpot(k,i)
 1010          continue
 10         continue

c           call matrix shake if (shakm)
            if (shakm) call mshakpt(step,step2)

c            if(procID.gt.0.and.procID.lt.proc-1) then 
c           Verlet stage 2
            tmp = 1.d0/step
            do 115 ii = 1,inofrz
               i = nofreez(ii)
               do 1115 k = 1,3
                  coor(k,i) = coor(k,i) + velo(k,i) 
                  velo(k,i) = velo(k,i) * tmp
 1115          continue
 115        continue
           

          

c           correct new position to satisfy constraints
c           call Correct Rigid Body Motion (CRBM)
c           AFTER THIS, BONDS WILL BE SLIGHTLY DISTORTED
c           Tony

      
            call crbm(coor,scalar,sigma,grdcmx,grdcmy,grdcmz
     1           ,grdlx,grdly,grdlz,divms,pmn,npt,nselec,pointr
     1           ,istp,ntest,ucrb,debug,udata)

c           compute new force
c        write(6,*) 'I am here 2'
c        do i = 1,npt
c        write(6,*) i,coor(1,i),coor(2,i),coor(3,i)
c        enddo
            call eforce()
c        write(6,*) 'I am here 2 after'
            



c           Verlet stage 3
c           calculate "free" velocities
            do 20 ii = 1,inofrz
               i = nofreez(ii)
               do 120 k = 1,3
                  velo(k,i) = velo(k,i) - fact2(i)*dpot(k,i)
 120           continue
 20         continue

       

c             endif     

c           deb for mshk of TIP3 ----velocity correction--------
            if (shakm) call mshakvl(step)


c           VELOCITY SCALING
c           Two options:
c           if the run is normal MD, do the following:
c           if the temperature is more than 20 degrees different
c           from the assigned temperature: scale the velocities
            tmpr = hot_umbr ( velo, dmass, npt, ndegf, massw )
            if ( (nsvel.ne.0) .and. (istp/nsvel*nsvel.eq.istp)) then
               vfac = dsqrt ( gtemp / tmpr )

c              scale the velocities
               do 989 ii = 1,inofrz
                  i = nofreez(ii)
                  do 1989 k = 1,3
                     velo(k,i) = velo(k,i) * vfac
 1989             continue
 989           continue

               old_tmpr = tmpr
               tmpr = hot_umbr ( velo, dmass, npt, ndegf, massw )
c               write(uene,1001)old_tmpr,tmpr,vfac
1001           format(1x,' old temp, new temp, velo scal fac '
     1          ,3(1x,f9.3))
  

            end if
            


c           energy facts
            if ( ( nene .eq. 1 ) .or. ( istp / nene * nene ) .eq. (
     1           istp  ) ) then
               write(uene,1240) istp 
1240            format(1x,'At dynamics step',i8)
c               write(uene,1005)  iens, istp
c1005           format(1x,'-----IENS, istp: ',2i6)
               write(uene,25) tmpr
 25            format(1x,'current temperature is ' , f9.3)
               write(uene,26) e_total + factor * tmpr * ndegf / 2
 26            format(1x,'total energy (kin + pot) = ', f12.3)
               call wener(uene)
               call flush(uene)
            end if

            
 100     continue 
c        END INNER INTEGRATION LOOP
          idyna = idyna + smlp

         
          goto 110

c        END OUTER INTEGRATION LOOP



c     ----------------------------------------------------
c     If the flow of control reaches this point, then the
c     trajectory is done, either by mxstps or by convergence
c     Now record some data.
 7800    continue


         write (*,1006) iens
1006     format(1x,'TRAJECTORY DONE:',1x,i6)

         
         if ( istp .eq. 1 ) then
            fpt = sflag * 1.d-3
         else if (istp .gt. mxstps) then
            fpt = mxstps 
         else if ( frac .lt. 0 ) then
            fpt = sflag * ( istp - 1 )
         else
            fpt = sflag * ( istp - 2 + frac )
         end if

c     write mlst number and final distance
         write(*,1007) tmlst
1007     format(1x,'tmlst:',i6)
         if ( udst .gt. 0 ) write(udst,1008) dst
1008     format(1x,' dst:',f10.4)




 7900 continue
c     END ENSEMBLE LOOP
      
      
      return
      end
