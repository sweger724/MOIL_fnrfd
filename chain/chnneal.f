        subroutine chnneal(igrid,npri,nstep,temper,
     1          pointr,ipick,nselec,d0,e0,e1,grdcmx,grdcmy,grdcmz,
     2          grdlx,grdly,grdlz,r,dv,udata,ucrd,
     3          nwcrd,ntest,debug,scalar,sigma,sigmav,
     4          gamma,rho,lambda,fixend,nlist,step)
        implicit none
c
        double precision gamma,rho,lambda
        double precision temper
        integer igrid,npri,nstep,nselec,nlist
        integer ucrd,udata,nwcrd,ntest
        logical debug,fixend
c
c igrid  -  number of chain monomers
c npri   -  print some useful(?) data each NPRI steps
c nstep  -  number of minimization steps 
c nselec -  number of selected particles, on which the chain constraints
c               are imposed.
c ucrd   -  unit number of file on which the coordinates (path format)
c               are written.
c udata  -  write info on the run on unit UDATA
c nwcrd  -  write coordinates on ucrd each NWCRD steps
c ntest  -  test the constraints each NTEST steps
c debug  - if .true. print a LOT of debugging info
c
        double precision grdcmx(*),grdcmy(*),grdcmz(*)
        double precision grdlx(*),grdly(*),grdlz(*)
        double precision r(3,*),dv(3,*)
        double precision d0(*),e0(*),e1(*)
        double precision scalar(*),sigma(*),sigmav(*)
c
c dmass - double precision mass vector (m(i=1,npt)(j),j=1,3*igrid)
c divms - double precision 1/mass vector(1/m(i=1,npt)(j),j=1,3*igrid)
c
c d0         -  vecotor of distances between i,i+1 & i,i+2 pairs
c e0         -  vecotor of the individual energies of the monomers
c e1         -  work vector for the polymer energy routine (sd)
c grdcm[x-z] -  gradient of center of mass constraints
c grdl[x-z]  -  gradient of infitesimal rotation constraints
c r - coordinates (rx(i=1,npt),ry(i=1,npt),rz(i=1,npt)(j=1,igrid))
c dv- forces      (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c
        integer pointr(*),ipick(*)
c
c pointr - a pointer to the selected particles
c          (which are subject to chain const.
c


C
C  next common blocks are used in the energy calculation.
C  list, debugging, coordinate transfer.. etc
C
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
c
c local
        integer i,j,k,npt3,ndegf,middle
        integer ncall
        double precision sener
        double precision velocity(3,lgrid*maxpt)
        double precision dv_old(3,lgrid*maxpt)
        double precision dtemper,temp_init
        double precision step,factor

        ncall  = 0
        npt3  = npt*3
        ndegf = (igrid-2)*npt3
        middle= igrid/2*npt
        temp_init = temper
        dtemper = temper/(nstep+1)

        write(udata,*)' Linear annealing will be performed'


c initialise velocities
        call draw_normal(velocity,ndegf,temper)
c make first energy/force call
          call sds(sener,dv,r,d0,e0,e1,nselec,pointr,ipick,npt,
     1          ndegf+2*npt3,gamma,rho,lambda,fixend,debug)
c project out from gradient rigid body forces
          j = 1
          do 1 i=1,igrid-2
           j = j + npt
           call crbm(dv(1,j),scalar,sigmav,grdcmx,grdcmy,grdcmz,
     1          grdlx,grdly,grdlz,ptms,npt,igrid,nselec,pointr,
     2          debug,stdo)
1       continue

c start moving
        do 10 k=1,nstep

         do 2 i=npt+1,(igrid-1)*npt
c       save old forces
           dv_old(1,i) = dv(1,i)
           dv_old(2,i) = dv(2,i)
           dv_old(3,i) = dv(3,i)
c       generate new coordiantes
           r(1,i) = r(1,i) + velocity(1,i-npt)*step
     1           - 0.5*step*step*dv(1,i)
           r(2,i) = r(2,i) + velocity(2,i-npt)*step
     1           - 0.5*step*step*dv(2,i)
           r(3,i) = r(3,i) + velocity(3,i-npt)*step
     1           - 0.5*step*step*dv(3,i)
2         continue
c       compute new forces
          call sds(sener,dv,r,d0,e0,e1,nselec,pointr,ipick,npt,
     1          ndegf+2*npt3,gamma,rho,lambda,fixend,debug)
c project out from gradient rigid body forces
          j = 1
          do 3 i=1,igrid-2
           j = j + npt
           call crbm(dv(1,j),scalar,sigmav,grdcmx,grdcmy,grdcmz,
     1          grdlx,grdly,grdlz,ptms,npt,igrid,nselec,pointr,
     2          debug,stdo)
3       continue
c current temperature scaling facture
         factor = dsqrt(1.d0 - dtemper/temper)
         temper = temper - dtemper
c compute new velocities
          do 4 i=1,(igrid-2)*npt
           velocity(1,i) = velocity(1,i) - 
     1          0.5d0*step*(dv_old(1,i+npt)+dv(1,i+npt))
           velocity(1,i) = velocity(1,i)*factor
           velocity(2,i) = velocity(2,i) - 
     1          0.5d0*step*(dv_old(2,i+npt)+dv(2,i+npt))
           velocity(2,i) = velocity(2,i)*factor
           velocity(3,i) = velocity(3,i) - 
     1          0.5d0*step*(dv_old(3,i+npt)+dv(3,i+npt))
           velocity(3,i) = velocity(3,i)*factor
4         continue
          write(udata,*)' step #   spw ',k,sener
10      continue

c
c
c save in unit udata the energies of individual monomers
c
           write(udata,*)' Writing coordinate at step ',step
           k  = - npt 
           do 5 i=1,igrid
                k = k + npt
                write(ucrd)e0(i),(r(1,j),j=k+1,k+npt),
     1                  (r(2,j),j=k+1,k+npt),(r(3,j),j=k+1,k+npt)
5          continue
        return
        end
