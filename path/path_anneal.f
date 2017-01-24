        subroutine path_anneal(ncalls,temper,constant_tmp,
     1          pointr,ipick,nselec,d0,grdcmx,grdcmy,grdcmz,
     2          grdlx,grdly,grdlz,udata,ucrd,
     3          debug,scalar,sigmav,
     4          nlist,s,dtopt,clo,dclo,npri)
c
        implicit none

        double precision temper
        integer ncalls,nselec,nlist,jloc
        integer ucrd,udata,npri
        logical debug

        character*80 filename

c npri   -  print some useful(?) data each NPRI steps
c ncalls  -  number of minimization steps
c nselec -  number of selected particles, on which the chain constraints
c               are imposed.
c ucrd   -  unit number of file on which the coordinates (path format)
c               are written.
c udata  -  write info on the run on unit UDATA
c nwcrd  -  write coordinates on ucrd each NWCRD steps
c debug  - if .true. print a LOT of debugging info
c
        double precision grdcmx(*),grdcmy(*),grdcmz(*)
        double precision grdlx(*),grdly(*),grdlz(*)
        double precision d0(*)
        double precision scalar(*),sigmav(*)
c
c dmass - double precision mass vector (m(i=1,npt)(j),j=1,3*(pseg+2))
c divms - double precision 1/mass vector(1/m(i=1,npt)(j),j=1,3*(pseg+2))
c
c d0         -  vecotor of distances between i,i+1 & i,i+2 pairs
c grdcm[x-z] -  gradient of center of mass constraints
c grdl[x-z]  -  gradient of infitesimal rotation constraints
c r - coordinates (rx(i=1,npt),ry(i=1,npt),rz(i=1,npt)(j=1,pseg+2))
c dv- forces      (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,pseg+2))
c
        integer pointr(*),ipick(*)
c
c pointr - a pointer to the selected particles
c          (which are subject to chain const.
c

      double precision time

C  next common blocks are used in the energy calculation.
C  list, debugging, coordinate transfer.. etc

        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/SGB.BLOCK'
        include 'COMMON/PATH.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
c
c local
        integer i,j,k,ndegf
        integer step_index,ii,l,n,npick
        integer ncall, utmp
        double precision s,dtopt,ss(4)
        double precision dr(3,lgrid*maxpt)
        double precision dtemper,temp_init
        logical constant_tmp
        double precision step,factor,drr,err
        double precision rms,clo,dclo,S0,S1,Snew

        write (6,*) 'in sdel_test pseg', pseg

        ncall  = 0
c        npt = 2
c        igrid=3
        ndegf = 3*npt*(pseg+2)

        call UpdateLists(r,npt,pseg)

        call Get_Dls(pseg,npt,r,d0,ipick)
        call GetDave(dave,d0)

c start moving
        do 10 step_index=1,ncalls
            call Communicate_Positions(npt,r,.true.)
            call draw_normal(dr,ndegf,1.d-10)
            drr = 0.d0
            do 6 j=1,npt
              do 6 l=1,3
               dr(l,j)=0.d0
C              coor(l,j)=r(l,j+npt)*massfac(j)
               dr(l,j+(pseg+1)*npt)=0.d0
6           continue

C           if (esymyes) call squeeze()
C            call nbondm()
C           if (esymyes) call syminit()

            do 7 j=2,pseg+1
               k=(j-1)*npt
               do 7 i=1,npt
                 drr = drr + dr(1,i+k)**2 + dr(2,i+k)**2 + dr(3,i+k)**2
7           continue
            drr = dsqrt(drr)
            write(6,*)'drr: ',drr
            call sds(S,d0,ipick,npt,clo)
            S0=S
C           call eforce()
C           S0 = e_total
C           call Get_EnergyDerivatives(pseg,npt,r,dv,e0)
C           S0=e0(2)
            do 5 j=2,pseg+1
               k=(j-1)*npt
               do 4 i=1,npt
                 r(1,i+k) = r(1,i+k) + dr(1,i+k)
                 r(2,i+k) = r(2,i+k) + dr(2,i+k)
                 r(3,i+k) = r(3,i+k) + dr(3,i+k)
C                coor(1,i)=r(1,i+npt)*massfac(i)
C                coor(2,i)=r(2,i+npt)*massfac(i)
C                coor(3,i)=r(3,i+npt)*massfac(i)
4              continue
5           continue
            call Get_Dls(pseg,npt,r,d0,ipick)
            call GetDave(dave,d0)
            call Communicate_Positions(npt,r,.true.)
            Snew=S0
            
            do 11 j=2,pseg+1
              k=(j-1)*npt
              do 12 i=1,npt
                Snew=Snew+dsall(1,i+k)*dr(1,i+k)+dsall(2,i+k)*dr(2,i+k)
     1              + dsall(3,i+k)*dr(3,i+k)
C                  Snew=Snew+(dpot(1,i)*dr(1,i+k)+dpot(2,i)*dr(2,i+k)
C     1                  + dpot(3,i)*dr(3,i+k))*massfac(i)
12               continue            
11           continue
            call sds(S,d0,ipick,npt,clo)    
            S1=S
C           call Get_EnergyDerivatives(pseg,npt,r,dv,e0)
C           S1=e0(2)
C           call eforce()
C           S1=e_total
            drr = drr**2


C           do i=1,(pseg+2)*npt
C             do l = 1 , 3
C                 r(l,i) = r(l,i) - dr(l,i)
C                 dv(l,i) = dv(l,i) - dr(l,i)
C             end do 
C           end do

            if (paral) then
                  ss(1)=S0
                  ss(2)=S1
                  ss(3)=Snew
                  ss(4)=drr
                  if(.not.first) then
                     call Send_Double(ss,4,0,procID)
                  endif
                  if (first) then
                     do n=1,numprocs-1
                        call Recv_double(ss,4,n,n)
                        S0=S0+ss(1)
                        S1=S1+ss(2)
                        Snew=Snew+ss(3)
                        drr = drr + ss(4)
                     enddo
                  endif
               endif
               drr = dsqrt(drr)
               if (first) then
                   err=dabs(Snew-S1)/drr**2
                   write(6,*)'S0,Snew,S1',S0,Snew,S1,err
               endif
 10     continue


C summarizing linear cooling, printing out a path
             
        return
        end
