        subroutine dyna_free(temp,tmpr,step,step2,beta,partit,
     1          exp2,qrot,igrid,numpth,npri,nstep,nsvel,pointr,
     2          nselec,dmass,divms,grdp,pstep,grdcmx,grdcmy,grdcmz,
     3          grdlx,grdly,grdlz,
     4          dpold,fact1,fact2,udata,ucrd,
     4          nwcrd,ntest,icol,newv,irand,debug,scalar,sigma,
     5          nlist,massw,ufrc,shakm,tolmshk,toiyes,uene)
Cdeb-- for mshk of TIP3--shakm and tolmshk added------
Cdeb tfac was changed to step2 as tfac was not being used.
c
        double precision temp,tmpr,step,step2,beta,partit,exp2,qrot
        integer igrid,numpth,npri,nstep,nsvel,nselec,ufrc
        integer udata,ucrd,nwcrd,ntest,icol,newv,irand,nlist
        logical debug,massw
        integer uene
c
c temp   -  assigned temperature
c tmpr   -  actually calculated temperature (from kinetic energy)
c               using function hot_umbr
c step   -  time integration step (in AMKA units)
C step2  -  step*step/2.d0
c tfac   -  conversion factor for time from PS to AMKA units
c igrid  -  number of segments
c numpth    current number along reaction path
c npri   -  print some useful(?) data each NPRI steps
c nstep  -  number of integration steps 
c nsvel  -  number of steps between scaling of the velocities
c nselec -  number of selected particles, on which the chain constraints
c               are imposed.
c ucrd   -  unit number of file on which the coordinates (path format)
c               are written.
c udata  -  write info on the run on unit UDATA
c nwcrd  -  write coordinates on ucrd each NWCRD steps
c ntest  -  test the constraints each NTEST steps
c icol   -  use ICOL steps for equilibratio (before collecting data)
c newv   -  assign new velocities each NEWV steps
c irand  - a seed for a random nmber generator
c debug  - if .true. print a LOT of debugging info
c
c
        double precision dmass(*),divms(*)
        double precision grdcmx(3,*),grdcmy(3,*),grdcmz(3,*)
        double precision grdlx(3,*),grdly(3,*),grdlz(3,*)
        double precision dpold(3,*)
        double precision fact1(*),fact2(*)
        double precision scalar(*),sigma(*),pstep(3,*),grdp(3,*)
c
c dmass - double precision mass vector (m(i=1,npt)
c divms - double precision 1/mass vector(1/m(i=1,npt)
c
c grdcm[x-z] -  gradient of center of mass constraints
c grdl[x-z]  -  gradient of infitesimal rotation constraints
c coor - coordinates (coor(i,j) i=x,y,z j=1,...,npt )
c velo - velocities (velo(i,j) i=x,y,z j=1,...,npt)
c dpot -  forces  (dpot(i,j) i=x,y,z j=1,..,npt)
c dpold - old forces 
c fact[1-2]  - constant vectors useful for integration (see freee)
c
        integer pointr(*)
c
c pointr - a pointer to the selected particles
c          which are subject to const.
c

c And now a LOT of neccessary common blocks
c

C
C  next common blocks are used in the energy calculation.
C  list, debugging, coordinate transfer.. etc
C
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/LINE.BLOCK'
Cdeb----for mshk of TIP3 ---------------------
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK' 
c
c local
c
c parameters to caluculate tensor on intertia at the prime and the
c remote points.
c i[x-z][x-z]  - the tensor of inertia before diagonalization
c a1,a2,a3     - coefficient for the third order resulting equation
c ia[a-c][1-2] - the principal moment of inertia [a-c] of the prime
c       system [1] and the remote system [2].
c The average is now taking the form :
c
c DELTA F = < exp {-b(V2-V1)} [(ia2 ibb2 ic2)/(ia1 ibb1 ic1)] ^1/2 > (1)
c
c note that this formula is restricted to systems which are
c rotationally and translationally invariant on the average.
c
        double precision ixx,iyy,izz,ixy,iyz,ixz
        double precision a1,a2,a3
        double precision tcm,xcm,ycm,zcm,xtmp,ytmp,ztmp
        double precision DET,R,Q,theta
        double precision ia1,ibb1,ic1,ia2,ibb2,ic2
        double precision exp1,expf,escond,fff
c local
        integer i,j,npt3,smlp,idyna,istp,lend
Cdeb--------------------------------
        integer ii, ndegf
        logical shakm,toiyes
        double precision tolmshk, tmp
Cdeb-------------------------------- 
        double precision hot_umbr,kinet,enetot,vfac,sigmav(7)
        double precision pi


        pi = 4.d0*datan(1.d0)
        npt3  = npt*3
        smlp  = 1
Cdeb---------Substract shake degrees of freedom-------------
        if (shakm) then
            ndegf = 3*inofrz - nshakm 
        else
            ndegf = 3*inofrz
        endif
Cdeb--------------------------------------------------------
c
c Data for velocities constraints.
c sigmav for the velocities is zero (in contrast to the coordinates)
c
                do 1 i=1,7
                 sigmav(i)=0.d0
1               continue

Cdeb-------------------------------------------------------
c
         if (esymyes) call squeeze()
         call nbondm()
         if (esymyes) call syminit()


Cdeb     call wener(stdo)
Cdeb            call nbondm()
                call eforce()
Cdeb------------------------------------------------------
C@@@@@@@@@@
C               call wener(stdo)
C               write(*,*)' pstep '
C@@@@@@@@@@
c
c calculate the energy for the remote point along the path
c (we jump back and forth from coor since the energy
c routines are working ONLY on coor which they get through common)
c
        do 6 j=1,nselec
                i=pointr(j)
                coor(1,i)=coor(1,i)+pstep(1,j)
                coor(2,i)=coor(2,i)+pstep(2,j)
                coor(3,i)=coor(3,i)+pstep(3,j)
C@@@@@@@@@@
C               write(*,*)j, pstep(1,j),pstep(2,j),pstep(3,j)
C@@@@@@@@@@
6       continue
      if (debug) then
        write(*,*)' x '
        write(*,1000)(coor(1,pointr(j)),j=1,nselec)
        write(*,*)' y '
        write(*,1000)(coor(2,pointr(j)),j=1,nselec)
        write(*,*)' z '
        write(*,1000)(coor(3,pointr(j)),j=1,nselec)
      end if
1000  format(1x,8(f9.6)) 
1001  format(20x,3(f15.6))
c ***************************************************
c CALL ENERGY 
c ***************************************************


        call eforce()
C@@@@@@@@@@
Cdeb    call wener(stdo)
c       stop
C@@@@@@@@@@
      if (debug) then
        write(*,*)' FIRST POINT'
        write(*,*)' dx '
        write(*,1000)(dpot(1,pointr(j)),j=1,nselec)
        write(*,*)' dy '
        write(*,1000)(dpot(2,pointr(j)),j=1,nselec)
        write(*,*)' dz '
        write(*,1000)(dpot(3,pointr(j)),j=1,nselec)
      end if
      escond=e_total
      
c
c change the remote point back to normal
c
        do 7 j=1,nselec
                i=pointr(j)
                coor(1,i)=coor(1,i)-pstep(1,j)
                coor(2,i)=coor(2,i)-pstep(2,j)
                coor(3,i)=coor(3,i)-pstep(3,j)
7       continue
c *************************************************************
c CALL ENERGY
c *************************************************************


        call eforce()

      if (debug) then
        write(*,*)' SECOND POINT'
        write(*,*)' dx '
        write(*,1000)(dpot(1,pointr(j)),j=1,nselec)
        write(*,*)' dy '
        write(*,1000)(dpot(2,pointr(j)),j=1,nselec)
        write(*,*)' dz '
        write(*,1000)(dpot(3,pointr(j)),j=1,nselec)
      end if
      expf=-beta*(escond-e_total)
      write(*,*) ' escond -e_total expf ' , escond -e_total,expf
      if (abs(expf).gt.10) then
        write(udata,*)'*** warning from free energy path (freee)'
        write(udata,*)'*** expf is too large = ',expf
        write(udata,*)'*** consider using more points for the path'
      end if

      write(udata,*)' *** END OF TRAJECTORY INITIALIZATION ***'
c
c       
        if (nlist .ne. 0 ) smlp=nlist
c
c Start integration loop
c
        write(udata,2420)
2420    format(5x,'Kin E  : Pot (1) : DelEne :
     1 Tot Ene :  Tempr   : expo  :   Q ')

        do  100  idyna = 1,nstep,smlp
c
c do internal update each inbfrq for the non-bonded AND the images
c (if applicable). A single non-bonded list is generated according
c to the first structure.
c
         if (nlist.ne.0) then
c
Cdeb-------------------------------------------------------
         if (esymyes) call squeeze()
         call nbondm()
         if (esymyes) call syminit()
Cdeb     call nbondm()
Cdeb------------------------------------------------------
c
         end if
         lend=min(nstep,smlp-1+idyna)
         do 100 istp = idyna,lend
c
c copy old forces to temporary data
c
                do 9 i=1,npt
                        dpold(1,i)=dpot(1,i)
                        dpold(2,i)=dpot(2,i)
                        dpold(3,i)=dpot(3,i)
9               continue
c
c calculate a "free" (unconstrained) step
c
Cdeb-----------------FOR mshk of TIP3-------------------------
Cdeb            do 10 i=1,npt
Cdeb            do 10 ii = 1,inofrz
Cdeb             i = nofreez(ii)  
Cdeb             coor(1,i)=coor(1,i)+velo(1,i)*step-dpold(1,i)*fact1(i)
Cdeb             coor(2,i)=coor(2,i)+velo(2,i)*step-dpold(2,i)*fact1(i)
Cdeb             coor(3,i)=coor(3,i)+velo(3,i)*step-dpold(3,i)*fact1(i)
Cdeb10          continue
c
c Make a step for unfrozen particles (velocity is used to store the step)
c dividing finally by dt will recover the magnitude required for velocity
c calculations
          do 10 ii = 1,inofrz
           i = nofreez(ii)
           velo(1,i) =  velo(1,i)*step - fact1(i)*dpold(1,i)
           velo(2,i) =  velo(2,i)*step - fact1(i)*dpold(2,i)
           velo(3,i) =  velo(3,i)*step - fact1(i)*dpold(3,i)
10        continue

c call matrix shake if (shakm)
          if (shakm) then
              call mshakpt(step,step2)
          end if
          tmp = 1.d0/step
          do 115 ii = 1,inofrz
           i = nofreez(ii)
           coor(1,i) = coor(1,i) + velo(1,i)
           coor(2,i) = coor(2,i) + velo(2,i)
           coor(3,i) = coor(3,i) + velo(3,i)
           velo(1,i) = velo(1,i)*tmp
           velo(2,i) = velo(2,i)*tmp
           velo(3,i) = velo(3,i)*tmp
115       continue
c
Cdeb----------------------------------------------------------
c
c correct the new position in order to satisfy the constraints
c call Correct Rigid Body Motion (CRBM)
c
Cdeb------------Peptide Frozen--------------------------------------
Cdeb      call crbm_free(coor,scalar,sigma,grdcmx,grdcmy,grdcmz,
Cdeb     1                 grdlx,grdly,grdlz,divms,grdp,
Cdeb     2                 npt,nselec,pointr,istp,ntest,debug,udata)
Cdeb
Cdeb------------Peptide Frozen--------------------------------------
c find the energy at the remote point
c
        do 17 j=1,nselec
                i=pointr(j)
                coor(1,i)=coor(1,i)+pstep(1,j)
                coor(2,i)=coor(2,i)+pstep(2,j)
                coor(3,i)=coor(3,i)+pstep(3,j)
17      continue
c ***************************************************************
c CALL ENERGY
c ***************************************************************
        
        call eforce()

        escond=e_total
        if(debug) then
      write(*,*) 'current energy at the remote point - escond' ,escond
        end if
Cdeb-----------------------------
       if (toiyes) then
c
c calculate the tensor of inertia elements
c
c set first the center of mass position to zero
c
        xcm=0.d0
        ycm=0.d0
        zcm=0.d0
        tcm=0.d0
        do 921 j=1,nselec
                i=pointr(j)
                tcm=tcm+dmass(i)
                xcm=xcm+dmass(i)*coor(1,i)
                ycm=ycm+dmass(i)*coor(2,i)
                zcm=zcm+dmass(i)*coor(3,i)
921     continue
        xcm=xcm/tcm
        ycm=ycm/tcm
        zcm=zcm/tcm
        if (debug) then
         write(*,*)' Center of mass position xcm ycm zcm'
         write(*,*)xcm,ycm,zcm
        end if
        ixx=0.d0
        iyy=0.d0
        izz=0.d0
        ixy=0.d0
        ixz=0.d0
        iyz=0.d0
        do 922 j=1,nselec
                i=pointr(j)
                xtmp=coor(1,i)-xcm
                ytmp=coor(2,i)-ycm
                ztmp=coor(3,i)-zcm
                ixx=ixx+dmass(i)*(ytmp*ytmp+ztmp*ztmp)
                iyy=iyy+dmass(i)*(xtmp*xtmp+ztmp*ztmp)
                izz=izz+dmass(i)*(ytmp*ytmp+xtmp*xtmp)
                ixy=ixy-dmass(i)*xtmp*ytmp
                ixz=ixz-dmass(i)*xtmp*ztmp
                iyz=iyz-dmass(i)*ytmp*ztmp
922     continue
        if (debug) then
         write(*,*)' Tensor of inertia elements xx yy zz xy xz yz '
         write(*,*)ixx,iyy,izz
         write(*,*)ixy,ixz,iyz
        end if
c
c now calculate the coefficient of the third order polynomial
c (note that the solution must be real which makes life easier)
c
        a1=-(ixx+iyy+izz)
        a2=ixx*iyy+iyy*izz+ixx*izz-iyz*iyz-ixy*ixy-ixz*ixz
        a3=-(ixx*iyy*izz-ixx*iyz*iyz-iyy*ixz*ixz-izz*ixy*ixy
     1   +2.d0*ixy*iyz*ixz)
        if (debug) then
         write(*,*)' coefficients for the cubic equation '
         write(*,*)' X^3 + a1 X^2 + a2 X + a3 = 0 '
         write(*,*)a1,a2,a3
        end if
c
c solve the cubic equation
c
        Q = ( 3.d0*a2 - a1*a1 ) / 9.d0
        R = ( 9.d0*a1*a2 - 27.d0*a3 - 2.d0*a1*a1*a1 ) / 54.d0
        DET = Q*Q*Q + R*R


        if (debug) then
         write(*,*)' Q R DET'
         write(*,*)Q,R,DET
        end if
        if (DABS(DET).lt.1.d-6) then
         if (DABS(R).gt.1.d-10) then
          R = R / DABS(R) * (DABS(R)**(1.d0/3.d0))
         end if
         if (debug) write(*,*) ' R^1/3 = ',R
         ia2 =  2.d0*R - 1.d0/3.d0*a1
         ibb2 =  -R - 1.d0/3.d0*a1
         ic2 = ibb2
         if (debug) then
          write(*,*)'**SOLUTION** (DET=0) ia ib ic'
          write(*,*)ia2,ibb2,ic2
         end if
        else if (DET.lt.0.d0) then
         if (Q.gt.0.d0) then
          write(*,*)' brrr... Q>0 at TIO '
          stop
         end if
         R = R / dsqrt( -Q*Q*Q )
         theta = dacos(R) / 3.d0
         if (debug) write(*,*)' R theta ',R,theta
         ia2 = 2.d0 * DSQRT(-Q) * dcos(theta) - 1.d0/3.d0*a1
         ibb2 = 2.d0 * DSQRT(-Q) * dcos(theta + 2.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         ic2 = 2.d0 * DSQRT(-Q) * dcos(theta + 4.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         if (debug ) then
          write(*,*)'**SOLUTION** (DET<0) ia2 ibb2 ic2'
          write(*,*)ia2,ibb2,ic2
         end if
        else if (DET.gt.0.d0) then 
         if (debug) then
          write(*,*)' error in parameters '
          write(*,*)' tensor of inertia has two imagionary eigenvalues'
         end if
         stop
        else
         if (debug) then
          write(*,*)' DET is not on the real axis ????'
          write(*,*)' DET = ',DET
         end if
        end if

       endif
Cdeb------------------------------------------------------------
c
c go back to the major point
c
        do 18 j=1,nselec
                i=pointr(j)
                coor(1,i)=coor(1,i)-pstep(1,j)
                coor(2,i)=coor(2,i)-pstep(2,j)
                coor(3,i)=coor(3,i)-pstep(3,j)
18      continue
c
c find the new force
c

        call eforce()
Cdeb----------Skip the calculation of force component
        go to 182    
c calculate the component of the force that is parallel to the RC
c
        fff  = 0.d0
        do 181 j=1,nselec
         i   =pointr(j)
         fff = fff + dpot(1,i)*grdp(1,j) + dpot(2,i)*grdp(2,j)
     1           + dpot(3,i)*grdp(3,j)
181     continue
        write(ufrc,*)istp*step*4.2532764d-2,fff
        if(debug) then
        write(*,*) ' current energy at the original point ' , e_total
        end if
182     continue
Cdeb------------------------------------------------------------------
       if (toiyes) then 
c
c calculate the tensor of inertia elements
c
c set first the center of mass position to zero
c
        xcm=0.d0
        ycm=0.d0
        zcm=0.d0
        tcm=0.d0
        do 931 j=1,nselec
                i=pointr(j)
                tcm=tcm+dmass(i)
                xcm=xcm+dmass(i)*coor(1,i)
                ycm=ycm+dmass(i)*coor(2,i)
                zcm=zcm+dmass(i)*coor(3,i)
931     continue
        xcm=xcm/tcm
        ycm=ycm/tcm
        zcm=zcm/tcm
        if (debug) then
         write(*,*)' Center of mass position xcm ycm zcm'
         write(*,*)xcm,ycm,zcm
        end if
        ixx=0.d0
        iyy=0.d0
        izz=0.d0
        ixy=0.d0
        ixz=0.d0
        iyz=0.d0
        do 932 j=1,nselec
                i=pointr(j)
                xtmp=coor(1,i)-xcm
                ytmp=coor(2,i)-ycm
                ztmp=coor(3,i)-zcm
                ixx=ixx+dmass(i)*(ytmp*ytmp+ztmp*ztmp)
                iyy=iyy+dmass(i)*(xtmp*xtmp+ztmp*ztmp)
                izz=izz+dmass(i)*(ytmp*ytmp+xtmp*xtmp)
                ixy=ixy-dmass(i)*xtmp*ytmp
                ixz=ixz-dmass(i)*xtmp*ztmp
                iyz=iyz-dmass(i)*ytmp*ztmp
932     continue
        if (debug) then
         write(*,*)' Tensor of inertia elements xx yy zz xy xz yz '
         write(*,*)ixx,iyy,izz
         write(*,*)ixy,ixz,iyz
        end if
c
c now calculate the coefficient of the third order polynomial
c (note that the solution must be real whichmakes life easier)
c
        a1=-(ixx+iyy+izz)
        a2=ixx*iyy+iyy*izz+ixx*izz-iyz*iyz-ixy*ixy-ixz*ixz
        a3=-(ixx*iyy*izz-ixx*iyz*iyz-iyy*ixz*ixz-izz*ixy*ixy
     1   +2.d0*ixy*iyz*ixz)
        if (debug) then
         write(*,*)' coefficients for the cubic equation '
         write(*,*)' X^3 + a1 X^2 + a2 X + a3 = 0 '
         write(*,*)a1,a2,a3
        end if
c
c solve the cubic equation
c
        Q = ( 3.d0*a2 - a1*a1 ) / 9.d0
        R = ( 9.d0*a1*a2 - 27.d0*a3 - 2.d0*a1*a1*a1 ) / 54.d0
        DET = Q*Q*Q + R*R


        if (debug) then
         write(*,*)' Q R DET'
         write(*,*)Q,R,DET
        end if
        if (DABS(DET).lt.1.d-6) then
         if (DABS(R).gt.1.d-10) then
          R = R / DABS(R) * (DABS(R)**(1.d0/3.d0))
         end if
         if (debug) write(*,*) ' R^1/3 = ',R
         ia1 =  2.d0*R - 1.d0/3.d0*a1
         ibb1 =  -R - 1.d0/3.d0*a1
         ic1 = ibb1
         if (debug) then
          write(*,*)'**SOLUTION** (DET=0) ia1 ibb1 ic1'
          write(*,*)ia1,ibb1,ic1
         end if
        else if (DET.lt.0.d0) then
         if (Q.gt.0.d0) then
          write(*,*)' brrr... Q>0 at TIO '
          stop
         end if
         R = R / dsqrt( -Q*Q*Q )
         theta = dacos(R) / 3.d0
         ia1 = 2.d0 * DSQRT(-Q) * dcos(theta) - 1.d0/3.d0*a1
         ibb1 = 2.d0 * DSQRT(-Q) * dcos(theta + 2.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         ic1 = 2.d0 * DSQRT(-Q) * dcos(theta + 4.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         if (debug ) then
          write(*,*)'**SOLUTION** (DET<0) ia1 ibb1 ic1'
          write(*,*)ia1,ibb1,ic1
         end if
        else if (DET.gt.0.d0) then 
         if (debug) then
          write(*,*)' error in parameters '
          write(*,*)' tensor of inertia has two imagionary eigenvalues'
         end if
         stop
        else
         if (debug) then
          write(*,*)' DET is not on the real axis ????'
          write(*,*)' DET = ',DET
         end if
        end if

       endif
Cdeb----------------------------------------------------------------------
c
c calculate "free" velocities
c
Cdeb            do 20 i=1,npt
                do 20 ii = 1,inofrz
                 i = nofreez(ii) 
Cdeb             velo(1,i)=velo(1,i)-fact2(i)*(dpot(1,i)+dpold(1,i))
Cdeb             velo(2,i)=velo(2,i)-fact2(i)*(dpot(2,i)+dpold(2,i))
Cdeb             velo(3,i)=velo(3,i)-fact2(i)*(dpot(3,i)+dpold(3,i))
Cdeb---------Since velo was used to store increment in coor and then was divided by
Cdeb---------step size the term fact2*old_force is already there.
                 velo(1,i)=velo(1,i)-fact2(i)*dpot(1,i)
                 velo(2,i)=velo(2,i)-fact2(i)*dpot(2,i)
                 velo(3,i)=velo(3,i)-fact2(i)*dpot(3,i)
20              continue
c
Cdeb for mshk of TIP3 ----velocity correction--------------
c
 
        if (shakm) call mshakvl(step)

 
Cdeb-------------------------------------------------------
c
c correct the velocities to satisfy the constraints
c
Cdeb--------------Peptide is frozen------------------------------
Cdeb     call crbm_free(velo,scalar,sigmav,grdcmx,grdcmy,grdcmz,
Cdeb     1               grdlx,grdly,grdlz,divms,grdp,
Cdeb     2               npt,nselec,pointr,istp,ntest,debug,udata)
Cdeb--------------Peptide is frozen------------------------------
c
 
c Two options:
c if the run is normal MD, do the following:
c if the temperature is more than 20 degrees different from the
c assigned temperature: scale the velocities
c
         if (nsvel.ne.0) then
          if  (istp/nsvel*nsvel.eq.istp) then
Cdeb          tmpr = hot_umbr(velo,dmass,npt,npt3,massw)
              tmpr = hot_umbr(velo,dmass,npt,ndegf,massw)
Cdeb            if ((dabs(tmpr-temp).gt.20.d0))then
                  vfac=dsqrt(temp/tmpr)
Cdeb              write(udata,*)' *** scaling velocities at step ',istp
Cdeb              write(udata,*)' current temperature is ',tmpr
Cdeb              write(udata,*)' desired temperature is ',temp
Cdeb              do 989 i=1,npt
                  do 989 ii = 1,inofrz
                        i = nofreez(ii)
                        velo(1,i)=velo(1,i)*vfac
                        velo(2,i)=velo(2,i)*vfac
                        velo(3,i)=velo(3,i)*vfac
989               continue
Cdeb            end if
          end if
         end if
c
c the time to reassign velocities?!
c
        if (istp/newv*newv.eq.istp) then

         write(udata,*)' *** assigning new velocities at step ',istp
         write(udata,*)' irand = ',irand
c
c set velocities - initial conditions
c
         call velinit(temp,1,0)

Cdeb-------------Peptide is frozen----------------------------------------
c
c correct the newly assigned velocities to satisfy the constraints
c
Cdeb           call crbm_free(velo,scalar,sigmav,grdcmx,grdcmy,grdcmz,
Cdeb     1               grdlx,grdly,grdlz,divms,grdp,
Cdeb     2           npt,nselec,pointr,istp,ntest,debug,udata)
Cdeb-------------Peptide is frozen----------------------------------------
         end if
c
c if after thermalization collect data
c
        if (istp.gt.icol) then
         expf=-beta*(escond-e_total)
Cdeb     exp1=dsqrt((ia2*ibb2*ic2)/(ia1*ibb1*ic1))
Cdeb     qrot=qrot+exp1
Cdeb     if (debug) write(*,*) ' TOI RATIO ' , exp1 
         exp1=dexp(expf)
         partit=partit+exp1
         exp2=exp2+exp1*exp1
         if(debug) then
         write(udata,*) ' exp1,partit,exp2 ' , exp1,partit,exp2
         end if
         if (abs(expf).gt.7) then
                write(udata,*)'*** warning ***',' at step ',istp
                write(udata,*)'*** expf is too large = ',expf
        end if
        end if
c 
c printing out some useful information
c
        if ((istp/npri*npri.eq.istp).and.(istp.gt.icol)) then
        kinet=0.d0
Cdeb    do 88 i=1,npt
        do 88 ii = 1,inofrz
        i = nofreez(ii)   
        kinet=kinet+dmass(i)*(velo(1,i)*velo(1,i)+
     1          velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
88      continue
        kinet=kinet/2.d0
        enetot=e_total+kinet
        write(uene,2410) kinet,e_total,escond-e_total,enetot,tmpr,
     1          expf,partit
2410    format(1x,f9.2,1x,f9.2,1x,f7.3,1x,e9.3,1x,f9.3,1x,f9.3,1x,f9.3)
        end if
c
c each nwcrd steps write coordinates in path format on ucrd
c
        if (debug) write(*,*) ' istp nwcrd icol ' , istp,nwcrd,icol
        if (istp/nwcrd*nwcrd.eq.istp) then
         write(ucrd)e_total,(coor(1,i),i=1,npt),(coor(2,i),i=1,npt),
     1          (coor(3,i),i=1,npt)
        end if
100     continue
c
c ************ END INTEGRATION LOOP ***********************
c

        partit=partit/(nstep-icol)
        exp2=exp2/(nstep-icol)
        exp2=dsqrt((exp2-partit*partit)/(nstep-icol))
c *********************************************************
c PARTITION FUCTION AND FREE ENERGY DIFFERNCE ARE DONE
c *********************************************************

        return
        end

