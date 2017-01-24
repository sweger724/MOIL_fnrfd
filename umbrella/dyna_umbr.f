        subroutine dyna_umbr(temp,tmpr,step,tfac,xfc,
     1          q0,igrid,numpth,npri,nstep,nsvel,pointr,
     2          nselec,dmass,divms,grdp,pstep,grdcmx,grdcmy,grdcmz,
     3          grdlx,grdly,grdlz,
     4          dpold,fact1,fact2,udata,ucrd,
     4          nwcrd,ntest,icol,newv,irand,debug,scalar,sigma,
     5          nlist,massw)
c
        double precision temp,tmpr,step,tfac,xfc,q0
        integer igrid,numpth,npri,nstep,nsvel,nselec
        integer udata,ucrd,nwcrd,ntest,icol,newv,irand,nlist
        logical debug,massw
c
c temp   -  assigned temperature
c tmpr   -  actually calculated temperature (from kinetic energy)
c               using function hot_umbr
c step   -  time integration step (in AMKA units)
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
c udata  -  A unit to write q
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
c coor(3,*)  - coordinates
c velo(3,*)  - velocities
c dpot(3,*)  -  forces
c dpold(3,*) - old forces 
c fact[1-2]  - constant vectors useful for integration (see freee)
c
        integer pointr(*)
c
c pointr - a pointer to the selected particles
c          (which are subject to chain const.
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
c
c local
c
c parameters to caluculate tensor on intertia at the prime and the
c remote points.
c i[x-z][x-z]  - the tensor of inertia before diagonalization
c a1,a2,a3     - coefficient for the third order resulting equation
c ia[a-c][1-2] - the principal moment of inertia [a-c] of the prime
c       system [1] and the remote system [2].
c
        double precision ixx,iyy,izz,ixy,iyz,ixz
        double precision a1,a2,a3
        double precision tcm,xcm,ycm,zcm,xtmp,ytmp,ztmp
        double precision DDD,R,Q,theta,pi
        double precision ia1,ibb1,ic1
c local
        integer i,j,jj,jjj,npt3,smlp,idyna,istp,lend
        double precision hot_umbr,kinet,enetot,vfac,sigmav(7)
        double precision dtemp,tempi
        double precision xfc2,eu

        tempi = 1.d0
        dtemp = (temp-tempi)/icol
        temp  = tempi+dtemp
        call velinit(temp,1,0)
        pi    = 4.d0*datan(1.d0)
        npt3  = npt*3
        smlp  = 1
c
c Data for velocities constraints.
c sigmav for the velocities is zero (in contrast to the coordinates)
c
                do 1 i=1,7
                 sigmav(i)=0.d0
1               continue
                call nbondm()
                call eforce()

        eu = 0.d0
        xfc2 = 2.d0*xfc
        do 311 i=1,nselec
                j=pointr(i)
                eu = eu + ((coor(1,j)-coor2(1,j))*grdp(1,i) +
     1          (coor(2,j)-coor2(2,j))*grdp(2,i) +
     2          (coor(3,j)-coor2(3,j))*grdp(3,i)) 
311     continue
        do 312 i=1,nselec
                j=pointr(i)
                dpot(1,j) = dpot(1,j) + xfc2 * eu *grdp(1,i)
                dpot(2,j) = dpot(2,j) + xfc2 * eu *grdp(2,i)
                dpot(3,j) = dpot(3,j) + xfc2 * eu *grdp(3,i)
312     continue
        eu = xfc*eu*eu

                e_total=e_total+eu

                scalar(1)=0.d0
                do 888 j=1,nselec
                i=pointr(j)
                scalar(1)=scalar(1)+grdp(1,j)*(coor(1,i)-coor2(1,i))
     1          +grdp(2,j)*(coor(2,i)-coor2(2,i))
     2          +grdp(3,j)*(coor(3,i)-coor2(3,i))
888             continue
                write(udata,*) scalar(1)+q0


c
      write(stdo,*)' *** END OF TRAJECTORY INITIALIZATION ***'
c
        
        if (nlist .ne. 0 ) smlp=nlist
c
c Start integration loop
c

        do  100  idyna = 1,nstep,smlp
c
c do internal update each inbfrq for the non-bonded AND the images
c (if applicable). A single non-bonded list is generated according
c to the first structure.
c
         if (nlist.ne.0) then
          call nbondm()
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
                do 10 i=1,npt
                 coor(1,i)=coor(1,i)+velo(1,i)*step-dpold(1,i)*fact1(i)
                 coor(2,i)=coor(2,i)+velo(2,i)*step-dpold(2,i)*fact1(i)
                 coor(3,i)=coor(3,i)+velo(3,i)*step-dpold(3,i)*fact1(i)
10              continue
c
c correct the new position in order to satisfy the constraints
c call Correct Rigid Body Motion (CRBM)
c
          call crbm_umbr(coor,scalar,sigma,grdcmx,grdcmy,grdcmz,
     1             grdlx,grdly,grdlz,divms,grdp,
     2             npt,nselec,pointr,istp,ntest,debug,stdo)

c
c find the new force
c
        call eforce()


        eu = 0.d0
        do 211 i=1,nselec
                j=pointr(i)
                eu = eu + ((coor(1,j)-coor2(1,j))*grdp(1,i) +
     1          (coor(2,j)-coor2(2,j))*grdp(2,i) +
     2          (coor(3,j)-coor2(3,j))*grdp(3,i)) 
211     continue
        do 2211 i=1,nselec
                j=pointr(i)
                dpot(1,j) = dpot(1,j) + xfc2 * eu *grdp(1,i)
                dpot(2,j) = dpot(2,j) + xfc2 * eu *grdp(2,i)
                dpot(3,j) = dpot(3,j) + xfc2 * eu *grdp(3,i)
2211     continue
        eu = xfc*eu*eu
        e_total=e_total+eu
                if (istp.gt.icol) then
                 do 889 j=1,nselec
                  i=pointr(j)
                  jj=j+nselec
                  jjj=jj+nselec
                  scalar(1)=scalar(1)+grdp(1,j)*(coor(1,i)-coor2(1,i))
     1            +grdp(2,j)*(coor(2,i)-coor2(2,i))
     2            +grdp(3,j)*(coor(3,i)-coor2(3,i))
889              continue
                 write(udata,*) scalar(1)+q0,e_total-eu,eu
                end if
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
        DDD = Q*Q*Q + R*R


        if (debug) then
         write(*,*)' Q R DDD'
         write(*,*)Q,R,DDD
        end if
        if (dabs(DDD).lt.1.d-6) then
         if (dabs(R).gt.1.d-10) then
          R = R / dabs(R) * (dabs(R)**(1.d0/3.d0))
         end if
         if (debug) write(*,*) ' R^1/3 = ',R
         ia1 =  2.d0*R - 1.d0/3.d0*a1
         ibb1 =  -R - 1.d0/3.d0*a1
         ic1 = ibb1
         if (debug) then
          write(*,*)'**SOLUTION** (DDD=0) ia1 ibb1 ic1'
          write(*,*)ia1,ibb1,ic1
         end if
        else if (DDD.lt.0.d0) then
         if (Q.gt.0.d0) then
          write(*,*)' brrr... Q>0 at TIO '
          stop
         end if
         R = R / dsqrt( -Q*Q*Q )
         theta = dacos(R) / 3.d0
         ia1 = 2.d0 * dsqrt(-Q) * dcos(theta) - 1.d0/3.d0*a1
         ibb1 = 2.d0 * dsqrt(-Q) * dcos(theta + 2.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         ic1 = 2.d0 * dsqrt(-Q) * dcos(theta + 4.d0*pi/3.d0 )
     1     - 1.d0/3.d0*a1
         if (debug ) then
          write(*,*)'**SOLUTION** (DDD<0) ia1 ibb1 ic1'
          write(*,*)ia1,ibb1,ic1
         end if
        else if (DDD.gt.0.d0) then 
         if (debug) then
          write(*,*)' error in parameters '
          write(*,*)' tensor of inertia has two imagionary eigenvalues'
         end if
         stop
        else
         if (debug) then
          write(*,*)' DDD is not on the real axis ????'
          write(*,*)' DDD = ',DDD
         end if
        end if
c
c calculate "free" velocities
c
                do 20 i=1,npt
                 velo(1,i)=velo(1,i)-fact2(i)*(dpot(1,i)+dpold(1,i))
                 velo(2,i)=velo(2,i)-fact2(i)*(dpot(2,i)+dpold(2,i))
                 velo(3,i)=velo(3,i)-fact2(i)*(dpot(3,i)+dpold(3,i))
20              continue
c
c correct the velocities to satisfy the constraints
c
         call crbm_umbr(velo,scalar,sigmav,grdcmx,grdcmy,grdcmz,
     1           grdlx,grdly,grdlz,divms,grdp,
     2           npt,nselec,pointr,istp,ntest,debug,stdo)
c
c check conservation of energy
c
        if (istp/ntest*ntest.eq.istp)  then
                kinet = 0.d0
                do 8 i = 1,npt
                        kinet = kinet + dmass(i)*(velo(1,i)*velo(1,i)+
     1           velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
8               continue
                kinet = kinet*0.5d0
                enetot = e_total + kinet
        write(stdo,*)' at step ',istp,' total energy is ',enetot        
        write(stdo,*)' tot. pot. ',e_total,' total kine. ',kinet        
        end if
c
c Two options:
c if the run is normal MD, do the following:
c if the temperature is more than 20 degrees different from the
c assigned temperature: scale the velocities
c
         if (nsvel.ne.0 .or. istp.lt.icol) then
          if (istp.lt.icol) temp = temp+dtemp
          if  (istp/nsvel*nsvel.eq.istp .or. istp.lt.icol) then
                tmpr = hot_umbr(velo,dmass,npt,npt3,massw)
                if ((dabs(tmpr-temp).gt.20.d0))then
                  vfac=dsqrt(temp/tmpr)
                  write(stdo,*)' *** scaling velocities at step ',istp
                  write(stdo,*)' current temperature is ',tmpr
                  write(stdo,*)' desired temperature is ',temp
                  write(stdo,*)' scaling factor = ',vfac
                  do 989 i=1,npt
                        velo(1,i)=velo(1,i)*vfac
                        velo(2,i)=velo(2,i)*vfac
                        velo(3,i)=velo(3,i)*vfac
989               continue
         end if
        end if
       end if
c
c the time to reassign velocities?!
c
        if (istp/newv*newv.eq.istp) then

         write(stdo,*)' *** assigning new velocities at step ',istp
         write(stdo,*)' irand = ',irand
c
c set velocities - initial conditions
c
         call velinit(temp,1,0)

c
c correct the newly assigned velocities to satisfy the constraints
c
          call crbm_umbr(velo,scalar,sigmav,grdcmx,grdcmy,grdcmz,
     1           grdlx,grdly,grdlz,divms,grdp,
     2           npt,nselec,pointr,istp,ntest,debug,stdo)

         if (debug) then
          write(stdo,*)' reassigning velocities at ',istp,' irand ',irand
       write(stdo,103)((velo(j,i),j=1,3),i=1,npt)
103    format(1x,8(f9.4,1x),/)
         end if
        end if
c 
c printing out some useful information
c
        if(istp/npri*npri.eq.istp) then
         call wener(stdo)
         kinet=0.d0
         do 88 i=1,npt
          kinet=kinet+dmass(i)*(velo(1,i)*velo(1,i)+
     1          velo(2,i)*velo(2,i)+velo(3,i)*velo(3,i))
88       continue
         kinet=0.5d0*kinet
         enetot=e_total+kinet
         write(stdo,241) kinet,e_total,enetot
241   format(1x,' Kinetic Energy ',e11.4,/,1x,' Potential Energy (1)'
     2  ,f9.3,/,1x,' Total Energy ' ,e11.4)
         tmpr = hot_umbr(velo,dmass,npt,npt3,massw)
         write(stdo,242) tmpr
242      format(1x, 'current temperature is ',f9.3)
        end if
c
c each nwcrd steps write coordinates in path format on ucrd
c
        if (debug) write(*,*) ' istp nwcrd icol ' , istp,nwcrd,icol
        if (istp/nwcrd*nwcrd.eq.istp) then
         write(ucrd)e_total,((coor(j,i),i=1,npt),j=1,3)
        end if
100     continue
c
c ************ END INTEGRATION LOOP ***********************
c

        return
        end

