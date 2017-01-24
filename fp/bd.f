        subroutine bd(temp,gamma,step,mxstps,
     1    mlsti,npri,nens,pointr,
     2    nselec,divms,grdp,grdcmx,grdcmy,grdcmz,
     3    grdlx,grdly,grdlz,
     4    fact2,udata,urcrd,uwfpt,uwcrd,nwcrd,
     5    irand,debug,scalar,sigma,nequ,
     6    nlist,shakm,
     7    shakb,shakl,tolshk,uene,
     8    prp,prn,plp,pln,pmp,pmn,ntest,orth)
        
        implicit none
        double precision temp,step,fptdp
        integer mlsti,npri,nens,nselec,nwcrd,iens,mxstps
        integer udata,urcrd,uwfpt,uwcrd,uene
        integer irand,nlist,ntest,nequ
        logical debug
        logical orth
        double precision tmp_rand(6)

c       temp   -  assigned temperature
c       gamma  -  friction coefficient
c       step   -  time integration step (in AMKA units)
c       mxstsp -  maximum number of steps
c       tfac   -  conversion factor for time from PS to AMKA units
c       numpth    current number along reaction path
c       npri   -  print some useful(?) data each NPRI steps
c       nens  -  ensemble size
c       nselec -  number of selected constrained particles
c       urcrd  -  file containing traj starting points
c       uwfpt  -  (write) fpt file
c       uwcrd  -  (write) fp trajectory file (one trajectory written)
c       nwcrd  -  num timesteps between fp traj write
c       mxstps -  maximum number of integration steps to make,
c                 NOT COUNTING equilibration steps
c       ntest  -  check constraints every ntest steps
c       nequ   -  number of pre-equilibration steps
c       udata  -  write info on the run on unit UDATA
c       irand  - a seed for a random nmber generator
c       debug  - if .true. print a LOT of debugging info
c       orth  - constrain to hyperplane iff orth
        
        double precision divms(*)
        double precision grdcmx(3,*),grdcmy(3,*),grdcmz(3,*)
        double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
        double precision nranf, nsysf
        double precision fact2(*)
        double precision scalar(*),sigma(*),grdp(3,*)
c       double precision cm(3)

c       @

c       divms - double precision 1/mass vector(1/m(i=1,npt)
c       
c       grdcm[x-z] -  gradient of center of mass constraints
c       grdl[x-z]  -  gradient of infitesimal rotation constraints
c       coor - coordinates (coor(i,j) i=x,y,z j=1,...,npt )
c       velo - velocities (velo(i,j) i=x,y,z j=1,...,npt)
c       dpot -  forces  (dpot(i,j) i=x,y,z j=1,..,npt)
c       fact[1-2]  - constant vectors useful for integration (see freee)

c       pointr - a pointer to the selected particles
c       which are subject to const.
        integer pointr(*)

c       next common blocks are used in the energy calculation.
c       list, debugging, coordinate transfer.. etc
c       
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
c       deb----for mshk of TIP3 ---------------------
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK' 
c       
c       local
c       
        double precision ixx,iyy,izz,ixy,iyz,ixz
        double precision a1,a2,a3
        double precision tcm,xcm,ycm,zcm,xtmp,ytmp,ztmp
        double precision DET,R,Q,theta
        double precision ia1,ibb1,ic1,ia2,ibb2,ic2
        double precision exp1,expf,escond,fff
c       local
        integer i,j,npt3,smlp,idyna,istp
c       deb--------------------------------
        integer ii, ndegf, iat1, iat2
        logical shakm,shakb,shakl
        double precision tolshk
c       deb-------------------------------- 
        double precision hot_umbr,kinet,enetot,vfac,sigmav(7)
        double precision pi

        
c       dotl (dotr, dotm): left (right, mid) first passage dot product
        double precision dotl, dotr, dotm

c       plp, pln (prp, prn (pmp, pmn)) : left (right (mid)) fp plane
c       point and normal
        double precision plp(3,maxpt), pln(3,maxpt)
        double precision prp(3,maxpt), prn(3,maxpt)
        double precision pmp(3,maxpt), pmn(3,maxpt)

c       check initial structures
        double precision ss

c       fpt: in units of time step
        integer fpt

c       for file read
        double precision e0

c       track configurations
        double precision dst, dst2, coorstart(3,maxpt), dummy(3,maxpt)

c       functions
        double precision vscalar
        double precision hp, gamma, sigd

c       @
        integer k

c       --------------------------------------------------

        sigd = sqrt(2 * kboltzmann * temp * step / gamma)
        hp = step / gamma
c       write(*,*) 'hp sigd', hp, sigd

        pi = 4.d0*datan(1.d0)
        npt3  = npt*3
        smlp  = 1

c       Data for velocities constraints.
c       sigmav for the velocities is zero (in contrast to the coordinates)
        do 1 i=1,7
           sigmav(i)=0.d0
 1      continue


c       non-bonded list
        call nbondm()


        if (nlist .ne. 0) smlp = nlist
        if ( mxstps .gt. 0 ) mxstps = mxstps + nequ

c       BEGIN ENSEMBLE LOOP
        do 7900 iens = 1,nens

           fpt = 0
           fptdp = 0.d0
           write (*,*) '---------------------------'
           write (*,*) 'IENS: ', iens

c          read next structure into coor (PATH format)
           read(urcrd) e0, ((coorstart(j,i),i=1,npt),j=1,3)

           call vdcopy(coorstart, coor, npt3)

c          BEGIN INTEGRATION LOOP
           idyna = 1


c       ------------------- "while" loop begins here -------------------
c       do internal update each inbfrq for the non-bonded AND the images
c       (if applicable). A single non-bonded list is generated according
c       to the first structure.

 110       if (nlist.ne.0) call nbondm()

           
           do 100 istp = idyna,idyna+smlp-1

              if ((mxstps .gt. -1) .and. (istp .gt. mxstps)) then
                 write(*,*) 'DONE: istp, mxstps:', istp, mxstps
                 goto 7800
              end if


c             COMPUTE FORCE
              call eforce()
              if (istp/npri*npri .eq. istp) then
                 call wener(uene)
                 write(uene, *) 'nranf, nsysf', nranf, nsysf
              end if


c       write the FIRST trajectory to file (works for both oeq and fp)
              if ( iens .eq. 1 ) then
         if (((nwcrd.eq.1) .or. ((istp/nwcrd*nwcrd).eq.(istp-1)))
     2        .and. (istp.gt.nequ)) then
c                   write(*,*) 'WRITING fp...istp = ', istp
                    write(uwcrd) e_total, (coor(1,i),i=1,npt),
     1   (coor(2,i),i=1,npt), (coor(3,i),i=1,npt)
                 end if
              end if


c       Track progress, i.e., distance from starting point
              call dist2(coorstart, coor, npt3, dst2)
              call vecmin_fp(coor, pmp, npt3, dummy)
              dotm = vscalar(pmn, dummy, npt3)
              dst = sqrt(dst2 / npt3)
c             write(*,*) 'dst: ', dst
c       Check center of mass every once in a while
c             if ((istp-1)/1000*1000 .eq. (istp-1)) then
c             write(*,*) 'istp dist', istp, dst
c             call com(3, npt, coor, ptms, cm)
c             write(*,*) 'com', cm(1), cm(2), cm(3)


              if (.not. orth) then

c             Check for left first passage (unless at 1st milestone)
                 if (mlsti .gt. 1) then
                    call vecmin_fp(coor, plp, npt3, dummy)
                    dotl = vscalar(pln, dummy, npt3)
c                   write(*,*) 'dotl: ', dotl
                    if (dotl .lt. 0) then
                       fpt = - (istp - 1)
                       if ( fpt .eq. 0 ) fptdp = -1.d-3
                       goto 7800
                    end if
                 end if
                 
c                write(*,*) 'dotm: ', dotm

c            Check for right first passage
                 call vecmin_fp(coor, prp, npt3, dummy)
                 dotr = vscalar(prn, dummy, npt3)
c                write(*,*) 'dotr: ', dotr
                 if (dotr .gt. 0) then
                    fpt = istp - 1
                    if ( fpt .eq. 0 ) fptdp = 1.d-3
                    goto 7800
                 end if

              end if
              


c       BROWNIAN STEP
c       This loop is inefficient because gauss2 is
c       called twice too often. But leave as is for now. -TF

c             nsysf, nranf are the norms of the systematic
c             and random forces, respectively
              nsysf = 0
              nranf = 0
              do 647 ii = 1,inofrz
                 i = nofreez(ii)
                 call gauss2(6, tmp_rand)
                 velo(1,i) = -hp*dpot(1,i) + sigd*tmp_rand(1)
                 velo(2,i) = -hp*dpot(2,i) + sigd*tmp_rand(2)
                 velo(3,i) = -hp*dpot(3,i) + sigd*tmp_rand(3)
                 nsysf = nsysf + (hp * dpot(1,i)) ** 2
                 nsysf = nsysf + (hp * dpot(2,i)) ** 2
                 nsysf = nsysf + (hp * dpot(3,i)) ** 2
                 nranf = nranf + (sigd * tmp_rand(1)) ** 2
                 nranf = nranf + (sigd * tmp_rand(2)) ** 2
                 nranf = nranf + (sigd * tmp_rand(3)) ** 2
 647          continue
              nsysf = sqrt(nsysf)
              nranf = sqrt(nranf)


c             SHAKE (mostly copied from dyna.f)
              if (shakl.or.shakb .and. nshak.gt.0) then
                 do 649 k=1,nshak
                    iat1 = ishak1(k)
                    iat2 = ishak2(k)
                    cooref(1,k) = coor(1,iat1) - coor(1,iat2)
                    cooref(2,k) = coor(2,iat1) - coor(2,iat2)
                    cooref(3,k) = coor(3,iat1) - coor(3,iat2)
 649             continue
                 call shakept(100)
              end if

              
c             take a step
              do 648 ii = 1,inofrz
                 i = nofreez(ii)
                 coor(1,i) = coor(1,i) + velo(1,i)
                 coor(2,i) = coor(2,i) + velo(2,i)
                 coor(3,i) = coor(3,i) + velo(3,i)
 648          continue


c             correct the new position in order to satisfy the constraints
c             call Correct Rigid Body Motion (CRBM)
c             AFTER THIS CALL TO CRBM, BONDS WILL BE SLIGHTLY DISTORTED
              call crbm(coor,scalar,sigma,grdcmx,grdcmy,grdcmz,
     1             grdlx,grdly,grdlz,divms,grdp,
     2             npt,nselec,pointr,istp,ntest,orth,debug,udata)


 100       continue
           
           idyna = idyna + smlp

           goto 110
c          END INTEGRATION LOOP


 7800      if (.not. orth) then
              write(*,*) 'mlsti fpt dst', mlsti, fpt, dst
              write(*,*) 'dotl dotm dotr', dotl, dotm, dotr
              if ( fptdp .ne. 0.d0 ) then
                 write(uwfpt,*) fptdp
              else
                 write(uwfpt,*) fpt
              end if
           end if
 7900   continue
c       END ENSEMBLE LOOP


        return
        end
