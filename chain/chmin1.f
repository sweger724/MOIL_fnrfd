        subroutine chmin1(igrid,npri,nstep,
     1          pointr,ipick,nselec,d0,e0,e1,grdcmx,grdcmy,grdcmz,
     2          grdlx,grdly,grdlz,r,dv,tolg,udata,ucrd,
     3          nwcrd,ntest,debug,scalar,sigma,sigmav,
     4          gamma,rho,lambda,fixend,nlist,estred)
c
        double precision gamma,rho,lambda,estred
        double precision tolg
        double precision tmp1, tmp2
        integer igrid,npri,nstep,nselec,nlist
        integer ucrd,udata,nwcrd,ntest
        logical debug,fixend
        integer finished
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
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'

C  SGB MODIFICATIONS
        include 'COMMON/SGB.BLOCK'
C  END SGB MODIFICATIONS
c
c local
        integer i,j,k,npt3,ndegf,smlp,imin,middle
        integer ncall
        double precision sener

        ncall  = 0
        npt3  = npt*3
        ndegf = (igrid-2)*npt3
        middle= igrid/2*npt
        smlp  = 1
c
c Make a first energy call. Note that at first a call to nbondm is made
c 
          do 1 i=1,npt
                j = middle + i
                coor(1,i) = r(1,j)
                coor(2,i) = r(2,j)
                coor(3,i) = r(3,j)
1         continue
          call nbondm()
          if (nlist .eq. 0 ) then
           write(udata,*)' No update of non-bonded list '
           call pwl_chn(ndegf,nstep,tolg,estred,ncall,
     1          r(1,1),r(1,npt+1),dv(1,1),dv(1,npt),
     2          sener,pointr,ipick,igrid,nselec,
     3          d0,e0,e1,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,
     4          scalar,sigma,sigmav,gamma,rho,lambda,
     5          fixend,finished,npri)
          write(udata,100)ncall,sener
          write(udata,*)' Stru #    dist(i,i+1) dist(i,i+2)   e '
          do 15 i=1,igrid
           if (i.le.igrid-2) then
            tmp1 = d0(i)/dsqrt(dfloat(npt))
            tmp2 = dsqrt(d0(igrid-1+i)/dfloat(npt))
           else if (i.le.igrid-1) then
            tmp1 = d0(i)/dsqrt(dfloat(npt))
            tmp2 = 0.d0
           else
            tmp1 = 0.d0
            tmp2 = 0.d0
           end if
           write(udata,101)i,tmp1,tmp2,e0(i)
15        continue
          tmp1 = 0.d0
          do 16 j=1,igrid-2
           k = j*npt
           do 16 i=1,npt
            tmp1 = tmp1 + dv(1,i+k)*dv(1,i+k) 
     1          + dv(2,i+k)*dv(2,i+k) + dv(3,i+k)*dv(3,i+k)
16        continue
          tmp1 = dsqrt(tmp1/dfloat(ndegf))
          write(udata,102)tmp1
c
c
c save in unit udata the energies of individual monomers
c
          write(udata,*)' Writing coordinate at step ',ncall
          k  = - npt3
          do 17 i=1,igrid
                k = k + npt3
                write(ucrd)e0(i),(r(1,j),j=1,npt),
     1                  (r(2,j),j=1,npt),(r(3,j),j=1,npt)
17        continue
          else
           smlp=nlist
c
c Start minimization loop
c

           do  10  imin = 1,nstep,smlp
C SGB MODIFICATION
              if (imin.eq.1) then
                 do i = 1,igrid
                    alphac(i) = 1
                 end do
              else 
                 do i = 1,igrid
                    alphac(i) = 1
                 end do
              end if

C END SGB MODIFICATION
c
c regenerate nonbonded list according to middle structure eaxh smlp
c steps.
c
          do 2 i=1,npt
                j = middle + i
                coor(1,i) = r(1,j)
                coor(2,i) = r(2,j)
                coor(3,i) = r(3,j)
2         continue
          write(udata,*)' The non-bonded list will be updated each '
     1    ,smlp,' steps'
          call nbondm()
          finished=0
          call pwl_chn(ndegf,smlp,tolg,estred,ncall,
     1          r(1,1),r(1,npt+1),dv(1,1),dv(1,npt+1),
     2          sener,pointr,ipick,igrid,nselec,
     3          d0,e0,e1,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,
     4          scalar,sigma,sigmav,gamma,rho,lambda,
     5          fixend,finished,npri)
c
c each npri steps printout some info
c
           if (finished.eq.1) then
              goto 200
           end if
          write(udata,100)ncall,sener
100       format(1x,'chain energy after ',i6,' energy calls = ',f16.5)
          write(udata,*)' Stru #    dist(i,i+1) dist(i,i+2)   e '
          do 3 i=1,igrid
           if (i.le.igrid-2) then
            tmp1 = d0(i)/dsqrt(dfloat(npt))
            tmp2 = dsqrt(d0(igrid-1+i)/dfloat(npt))
           else if (i.le.igrid-1) then
            tmp1 = d0(i)/dsqrt(dfloat(npt))
            tmp2 = 0.d0
           else
            tmp1 = 0.d0
            tmp2 = 0.d0
           end if
           write(udata,101)i,tmp1,tmp2,e0(i)
101        format(1x,i6,1x,3(f14.5,1x))
3         continue
          tmp1 = 0.d0
          do 4 j=1,igrid-2
           k = j*npt
           do 4 i=1,npt
            tmp1 = tmp1 + dv(1,i+k)*dv(1,i+k) 
     1          + dv(2,i+k)*dv(2,i+k) + dv(3,i+k)*dv(3,i+k)
4         continue
          tmp1 = dsqrt(tmp1/dfloat(ndegf))
          write(udata,102)tmp1
102       format(//,1x,' Current gradient ',f10.5,//)
10      continue

 200    continue
c
c
c save in unit udata the energies of individual monomers
c
           write(udata,*)' Writing coordinate at step ',ncall
           k  = - npt 
           do 5 i=1,igrid
                k = k + npt
                write(ucrd)e0(i),(r(1,j),j=k+1,k+npt),
     1                  (r(2,j),j=k+1,k+npt),(r(3,j),j=k+1,k+npt)
5          continue
        end if
        return
        end
