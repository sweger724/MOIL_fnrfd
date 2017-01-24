      subroutine dyna_fp(temp,step,step2,pmlst,tmlst,nmlst,npri,nens 
     1     ,mxstps,nsvel,pointr,nselec,dmass,divms,grdcmx,grdcmy 
     1     ,grdcmz,grdlx,grdly,grdlz,fact1,fact2,udata,urcrd 
     2     ,uwfpt,uwfpp,uwcrd,uwpep,uene,nene,nwcrd,ntest,ptemp,nequ 
     2     ,newv,irand,debug,scalar,sigma,nlist,massw,shakm,shakb 
     3     ,shakl,nwpep,prp,prn,plp,pln,pmp,pmn,unitgrad,orth,tfp 
     3     ,udot,udot1,udst,umom,ucrb,umlst,upfr,upvl
     4     ,num_mlsts) 

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
c     npri   -  number of steps between diagnostic output
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
c     uwfpt  -  first passage time file unit num
c     uwfpp  -  first passage configs path file unit num
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
      integer npri,nsvel,nselec,iens,nens,mxstps
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
      integer pointr(*), l
      double precision factor
      data factor/1.98768d-3/
      logical ucrb
      
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


      character*20 name
      integer namel,level

c     local
      integer i,j,k,npt3,smlp,idyna,istp,sflag
      integer ii, ndegf, iat1, iat2,nselec3
      logical shakm,shakl,shakb
      double precision tmp
      double precision hot_umbr,vfac,sigmav(7),gtemp
      double precision pi
      double precision tmp_rand(6)

      double precision pepcrd(3,maxpt)
      
c     dotl (dotr, dotm): left (right, mid) first passage dot product
      double precision dotl, dotr, dotm, dotlp, dotrp, dotmp
      
c     plp, pln (prp, prn (pmp, pmn)) : left (right (mid)) fp plane
c     point and normal
      double precision plp(3,maxpt), pln(3,maxpt)
      double precision prp(3,maxpt), prn(3,maxpt)
      double precision pmp(3,maxpt), pmn(3,maxpt)
      double precision unitgrad(3,maxpt)
      double precision rms
      
c     fpt: in units of time step
      double precision fpt, frac
      
c     for file read
      double precision e0
      
c     track configurations
      double precision dst, coorstart(3,maxpt)
      double precision dummy(3,maxpt), tmp_mass(maxpt)
      integer npick3

c     functions
      double precision vscalar, dist2_fp

      double precision dp(6)

      name = 'dyna_fp'
      namel = 7
      frac = 0

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


      if ( tfp ) then
         if ( nequ .gt. 0 ) then
            level = 1
            call alert(name,namel,'neqstep not zero',40,level)
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

      write(*,*) 'ndegf: ', ndegf

c     deb--------------------------------------------------------

c     Data for velocities constraints.
c     sigmav for the velocities is zero (in contrast to the
c     coordinates)
      
      do 1 i=1,7
         sigmav(i)=0.d0
 1    continue
      
      if (nlist .ne. 0) smlp = nlist
      if (mxstps .gt. 0) mxstps = mxstps + nequ


c     print out first few lines of relevant milestones
c      write(*,*) 'pmp:'
c      write(*,*) ((pmp(k,i),k=1,3),i=1,12)
c      write(*,*) 'pmn:'
c      write(*,*) ((pmn(k,i),k=1,3),i=1,12)
c      if ( tfp ) then
         if ( pmlst .ne. nomlst ) then
            call getmlst(umlst,plp,pln,pmlst,num_mlsts,nselec,pointr)
c            write(*,*) 'plp:'
c            write(*,*) ((plp(k,i),k=1,3),i=1,12)
c            write(*,*) 'pln:'
c            write(*,*) ((pln(k,i),k=1,3),i=1,12)
         end if
         if ( nmlst .ne. nomlst ) then
            call getmlst(umlst,prp,prn,nmlst,num_mlsts,nselec,pointr)
c            write(*,*) 'prp:'
c            write(*,*) ((prp(k,i),k=1,3),i=1,12)
c            write(*,*) 'prn:'
c            write(*,*) ((prn(k,i),k=1,3),i=1,12)
         end if
c      end if
      

      
c     BEGIN ENSEMBLE LOOP
      do 7900 iens = 1,nens
         write (*,*) '----------------------------------------------'
         write (*,*) 'IENS: ', iens

         atmlst = tmlst
         dotm =  0.d0

c         if ( tfp ) then 
            natm = mod ( atmlst + num_mlsts, num_mlsts ) + 1
            patm = mod ( atmlst + num_mlsts - 2, num_mlsts ) + 1
            call getmlst(umlst,prp,prn,natm,num_mlsts,nselec,pointr)
            call getmlst(umlst,plp,pln,patm,num_mlsts,nselec,pointr)
            dotl = +1.d0
            dotr = -1.d0
            fpt = 0.d0
c         end if

c        read next structure
         read(urcrd) e0, ((coorstart(j,i),i=1,npt),j=1,3)
         call vdcopy ( coorstart, coor, npt3 )


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


         call eforce()
         

c        BEGIN INNER INTEGRATION LOOP
         do 100 istp = idyna,idyna+smlp-1

            if ( orth ) write(upfr,*) projf(dpot,unitgrad,nselec,pointr)

c           check for termination by maximum steps
            if ((mxstps .gt. -1) .and. (istp .gt. mxstps)) then
               write(*,*) ' DONE: istp, mxstps:', istp, mxstps
               go to 7800
            end if


c           set temperature
            if ( istp .le. nequ ) then
               gtemp = ptemp
            else
               gtemp = temp
            end if


            dotmp = dotm
            if ( tfp ) then
               dotlp = dotl
               dotrp = dotr
            end if


c           dot products
            npick3 = 3 * nselec
            call vecmin_fp ( coor, pmp, dummy, nselec, pointr )
            dotm = vscalar ( pmn, dummy, npick3 )
c            if ( tfp ) then 
               call vecmin_fp ( coor, plp, dummy, nselec, pointr )
               dotl = vscalar ( pln, dummy, npick3 )
               call vecmin_fp ( coor, prp, dummy, nselec, pointr )
               dotr = vscalar ( prn, dummy, npick3 )
c            end if
            
            
            
c     print out some useful information
            if (.not. tfp  ) then
                 write(*,999)dotl,dotr 
999             format(1x,' OEQ distances, left and right ',2(1x,e10.4))
            end if
            if ((istp .eq. 1) .or. (istp / npri * npri) .eq. istp) then
               write(*,*) ''
               write(*,1000)istp,dotl,dotm,dotr
1000           format(1x,'ISTP = ',i6,/,1x,
     1          ' distance from planes l m r ',3(e10.3,1x))
            end if


c           write peptide traj
            if ( ( nwpep .eq. 1 ) .or. ( ( nwpep .gt. 0 ) .and. ( (
     1           istp / nwpep * nwpep ) .eq. ( istp - 1 ) ) .and. (
     2           istp .gt. nequ ) ) ) then
               do k = 1,3
                  do j = 1,nselec
                     i = pointr(j)
                     pepcrd(k,j) = coor(k,i)
                  end do
               end do
               write(uwpep) e_total, ((pepcrd(k,i),i=1,nselec),k=1,3)
            end if



c           write the FIRST trajectory to file. So when orth
c           is .true. this if block writes the entire milestone
c           ensemble.
            if ( iens .eq. 1 ) then

               if ((mod(istp,nwcrd) .eq. 0) .and. (istp .gt. nequ)) then
C Initialize udot1 as 99 in fp.f and make changes here
                  if ( udot1 .ne. 99 ) write(udot1,*) dotl, dotm, dotr

c                 write structure
                  write(uwcrd) e_total,(coor(1,i),i=1,npt),
     1                 (coor(2,i),i=1,npt),(coor(3,i),i=1,npt)

c                 write distance travelled
                  if ( umom .gt. 0 ) then
                     dst = sqrt ( dist2_fp ( coorstart, coor, nselec,
     1                    pointr ) )
                     write(umom,*) dst
                  end if

c                 verify non-drift of center of mass
                  call com_fp ( coor, dmass, cm, nselec, pointr )
                  write(*,1002)  cm(1), cm(2), cm(3)
1002              format(1x,' com ',2(f10.4))

c                 verify rotational non-drift
                  
                  do k = 1, npt
                    do l = 1,3
                      dummy(l,k) = coor(l,k)
                    end do
                    tmp_mass(k) = ptms(k)
                    ptms(k) = 0.d0
                  end do
                  
                  do k = 1, nselec
                    ptms(pointr(k)) = tmp_mass(pointr(k))
!                    write(6,*) "KKK",k,pointr(k),tmp_mass(pointr(k))
                  end do
                  
                  call rmsd_weight(npt,coorstart,dummy,rms,.false.,ptms)

                  do k = 1,npt
                    ptms(k) = tmp_mass(k)
                  end do

                  write(*,*) 'rot matrix wrt initial point:'
                  write(stdo,*) ''
                  write(stdo,*)' rotation matrix '
                  do k = 1,3
                    write(stdo,3001)(rotat(k,j),j=1,3)
3001                format(1x,3(f10.4,1x))
                  end do
                  write(stdo,*) ''
     
               end if
            end if


c           Now check for first passage
            if ( tfp ) then

               change = .false.

c              RIGHT (NEXT) first passage
               sflag = 1
               if ( dotr .ge. 0 ) then
                  if ( natm .eq. nmlst ) then
                     frac = -dotrp / (dotr - dotrp)
                     if ( frac .lt. 0 ) then
                        write(*,1003) dotrp, dotr
1003                    format(1x,'dotrp, dotr:',2(1x,f10.4))
                        call alert(name, namel,'negative frac!',40,0)
                     end if
                     goto 7800
                  end if
                  change = .not. ( (atmlst .eq. num_mlsts) .and.
     1                             (nmlst .eq. nomlst) )
                  if ( change ) atmlst = natm
                  goto 777
               end if


c              LEFT (PREV) first passage
               sflag = -1
               if ( dotl .le. 0 ) then
                  if ( patm .eq. pmlst ) then
                     frac = dotlp / (dotlp - dotl)
                     if ( frac .lt. 0 ) then
                        write(*,1004)  dotlp, dotl
1004                    format(1x,'dotlp, dotl:',2(1x,f10.4))
                        call alert(name, namel,'negative frac!',40,0)
                     end if
                     goto 7800
                  end if
                  change = .not. ( (atmlst .eq. 1) .and.
     1                             (pmlst .eq. nomlst) )
                  if ( change ) atmlst = patm
               end if

               
 777           if ( change ) then
                  natm = mod ( atmlst + num_mlsts, num_mlsts ) + 1
                  call getmlst(umlst,prp,prn,natm,num_mlsts,
     1                 nselec,pointr)
                  patm = mod ( atmlst + num_mlsts-2, num_mlsts ) + 1
                  call getmlst(umlst,plp,pln,patm,num_mlsts,
     1                 nselec,pointr)
                  write(*,*) 'NOW AT MILESTONE', atmlst
               end if

            end if


c     ---------------------------------------------------
c     If the flow of control reaches this point, then the
c     trajectory has not yet terminated. So proceed to
c     the Verlet step.
c     ----------------------------------------------------


c           sample the velocities at the right temperature
            if ( (istp.eq.1) .or. ((newv.ne.0).and.(istp/newv*newv.eq.
     1           istp)) ) then
            call RLUXGO(223,irand,0,0)   
            call velinit (gtemp,1,0 )
            end if


c           correct newly assigned velos to satisfy constraints
            call crbm ( velo,scalar,sigmav,grdcmx,grdcmy
     1           ,grdcmz,grdlx,grdly,grdlz,divms,pmn,npt,nselec
     1           ,pointr,istp,ntest,ucrb,debug,udata )


            if ( tfp ) write(upvl,*) projf(velo,unitgrad,nselec,pointr)


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
            call eforce()

c           Verlet stage 3
c           calculate "free" velocities
            do 20 ii = 1,inofrz
               i = nofreez(ii)
               do 120 k = 1,3
                  velo(k,i) = velo(k,i) - fact2(i)*dpot(k,i)
 120           continue
 20         continue
            


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
               write(*,1001)old_tmpr,tmpr,vfac
1001           format(1x,' old temp, new temp, velo scal fac '
     1          ,3(1x,f9.3))
            end if
            


c           energy facts
            if ( ( nene .eq. 1 ) .or. ( istp / nene * nene ) .eq. (
     1           istp - 1 ) ) then

               write(uene,1005)  iens, istp
1005           format(1x,'-----IENS, istp: ',2i6)
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
c     trajectory is done, either by mxstps or by first
c     passage. Now record some data.
 7800    continue

         write (*,1006) iens
1006     format(1x,'TRAJECTORY DONE:',1x,i6)

         
         if ( istp .eq. 1 ) then
            fpt = sflag * 1.d-3
         else if ((istp .gt. mxstps) .and. (mxstps.gt.0)) then
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

c     write dot products
C Initialize udot as 99 in fp.f and make changes here
         if ( ( udot .ne. 99 ) .and. ( orth .or. tfp ) ) write(udot,
     1        1009) dotl, dotm, dotr
1009     format(1x,3(1x,f10.4))



c     first passage data
         if ( tfp ) then
            write(*,1011)  fpt
1011        format(1x,'fpt:',f11.3)
            write(uwfpt,*) fpt
            write(uwfpp) e_total, (coor(1,i),i=1,npt),
     1           (coor(2,i),i=1,npt), (coor(3,i),i=1,npt)
         end if


 7900 continue
c     END ENSEMBLE LOOP
      
      
      return
      end
