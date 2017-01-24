        subroutine sds(sener,dv,r,d0,e0,e1,nselec,pointr
     1          ,ipick,npt,ndegf,gamma,rho,lambda,fixend
     2          ,debug)
c
c  a subroutine to calculate the energy of a polymer of system copies
c  The polymer energy (sener) is giving by
c
c  S =        SUM V(i) +              gamma SUM (d(i,i+1) - <d>)^2 + 
c      monomer internal energies    nearest monomers harmonic attraction
c
c       +  rho/lambda SUM exp( -lambda  d(i,i+2)^2 / <d>^2  )
c        next nearest monomers repulsion
c
c for more details see Czerminski and Elber, Int.J.Quant.Chem
c
c
c *** Note - in the current implementation the selection does not work!

        double precision sener,gamma,rho,lambda
        integer nselec,npt,ndegf

        double precision dv(3,*),r(3,*),d0(*),e0(*),e1(*)
        integer pointr(*),ipick(*)
        logical fixend,debug

c
c common block for COORDinates and potential ENERGY derivatives
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/UNITS.BLOCK'

C      SGB MODIFICATIONS
        include 'COMMON/SGB.BLOCK'
C      END SGB MODIFICATIONS
c
c Local
c
        integer igrid,npt2,i,j,k,l
        double precision tmp,tmp1,d2w,gam2,rhol,rho2,dave,dave2
        double precision tmpx,tmpy,tmpz
        logical debug1
c       double precision dv1(10000)
c
c initialiaztion
c
        debug1 = .false.
        npt2  = 2*npt
        igrid = ndegf/(3*npt)
        gam2  = 2.d0*gamma
        rhol  = rho/lambda
        rho2  = 2.d0*rho
c@
C       debug1 = .true.

        sener = 0.d0
        do 1 i = 1,igrid*npt
                dv(1,i) = 0.d0
                dv(2,i) = 0.d0
                dv(3,i) = 0.d0
1       continue

        do 2 i = 1,2*igrid-3
                d0(i) = 0.d0
2       continue


        do  5 j = 2,igrid-1
C SGB MODIFICATIONS
           chainno = j
C END SGB MODIFICATIONS
                k = (j-1)*npt
                do 3 i = 1,npt
                        coor(1,i) =  r(1,k+i)
                        coor(2,i) =  r(2,k+i)
                        coor(3,i) =  r(3,k+i)
3               continue
c
c energy call for the internal energy of monomer j
c
                if (esymyes) call squeeze()
                call nbondm()
                if (esymyes) call syminit()
                call eforce()
c
c e0 stores the internal energy of the monomers
c e_total is obtained from ENER common block
c
                e0(j) = e_total
c         write(6,*) 'pote',j,e0(j)
                if (debug1) then
                 call init_wre(stdo)
                 call wener(stdo)
                end if
                sener = sener + e_total
c
c copy internal energy derivatives to the total gradient
c
C$DOIT IVDEP
                do 4 i = 1,npt
                        dv(1,k+i) = dpot(1,i)
                        dv(2,k+i) = dpot(2,i)
                        dv(3,k+i) = dpot(3,i)
4               continue

c *** End of internal energy calculations

5       continue
c
c calculate now distances between all i,i+1 and i,i+2 pairs
c and the corresponding energy terms
c
C$DOIT SCALAR
        do 56 j = 1,igrid-1
                k = (j-1)*npt
                l = k + npt
C$DOIT IVDEP
                do 55 i = 1,npt
                 if (ipick(i).gt.0) then
                 tmpx = r(1,i+l) - r(1,i+k)
                 tmpy = r(2,i+l) - r(2,i+k)
                 tmpz = r(3,i+l) - r(3,i+k)
                 d0(j) = d0(j) + tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
                 end if
55              continue
56              continue
        dave = 0.d0
        do 6 j = 1,igrid-1
                d0(j) = dsqrt(d0(j))
                dave  = dave + d0(j)
6       continue
        dave = dave / (igrid-1)
        dave2 = 1.d0/(dave*dave)
        
C$DOIT SCALAR
        do 8 j = 1,igrid-2
                k = (j-1)*npt
                l = k + 2*npt
C$DOIT IVDEP
                do 7 i = 1,npt
                 if (ipick(i).gt.0) then
                 tmpx = r(1,i+l) - r(1,i+k)
                 tmpy = r(2,i+l) - r(2,i+k)
                 tmpz = r(3,i+l) - r(3,i+k)
                 d0(j+igrid-1) = d0(j+igrid-1) +
     1                   tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
                 end if
7               continue
8       continue
        if (debug1) then
        write(*,*)' Inside sds coordinates: '
        do 81 j=1,igrid
          k = npt*(j-1)
         do 81 i=1,npt
          write(*,*)r(1,k+i),r(2,k+i),r(3,k+i)
81      continue

        write(*,*)(d0(j+igrid-1),j=1,3)
        end if
c
c       calculate contribution from nearest neighbours
c
        tmp = 0.d0
        do 9 j = 1,igrid-1
                tmp1 = (d0(j) - dave)
                tmp  = tmp + tmp1*tmp1
9       continue
        sener = sener + gamma*tmp

c
c       calculate contribution from next nearest neighbours
c       also save the exponents for derivatives calculations
c
        tmp  = 0.d0
        tmp1 = 0.d0
        d2w  =0.d0
        do 10 j = igrid,2*igrid-3
                e1(j) = dexp( - lambda*d0(j)/(dave*dave) )
                tmp   = tmp  + e1(j)
                d2w   = d2w + e1(j)*d0(j)
                e1(j) = e1(j)*dave2
10      continue
        sener = sener + rho/lambda*tmp
        d2w = rho2*dave2/((igrid-1)*dave)*d2w
c
c       calcuate average distance multiply by gam2 divided by specific
c       distances to save time later
c
        do 11 j=1,igrid-1
                e1(j) = gam2*(1.d0 - dave/d0(j))
11      continue
c
c calculate derivatives from nearest neighbour distance terms.
c
C$DOIT SCALAR
        do 12 j = 1,igrid-1
         k  = (j-1)*npt
         l  = k +npt
C$DOIT IVDEP
         do 12 i = 1,npt
          if (ipick(i).gt.0) then
          tmpx = e1(j)*(r(1,k+i)-r(1,l+i))
          tmpy = e1(j)*(r(2,k+i)-r(2,l+i))
          tmpz = e1(j)*(r(3,k+i)-r(3,l+i))
          dv(1,i+k)      = dv(1,i+k)      + tmpx
          dv(2,i+k)      = dv(2,i+k)      + tmpy
          dv(3,i+k)      = dv(3,i+k)      + tmpz
          dv(1,i+l)      = dv(1,i+l)      - tmpx
          dv(2,i+l)      = dv(2,i+l)      - tmpy
          dv(3,i+l)      = dv(3,i+l)      - tmpz
          end if
12      continue
c
c calculate derivatives from second nearest neighbours (repulsion)
c
        do 13 j = 1,igrid-2
         k = (j-1)*npt
         l = k+npt2
         tmp1 = rho2*e1(j+igrid-1)
C$DOIT IVDEP
         do 13 i = 1,npt
         if (ipick(i).gt.0) then
          tmpx = (r(1,i+k)-r(1,i+l))*tmp1
          tmpy = (r(2,i+k)-r(2,i+l))*tmp1
          tmpz = (r(3,i+k)-r(3,i+l))*tmp1
          dv(1,i+k) = dv(1,i+k) - tmpx
          dv(2,i+k) = dv(2,i+k) - tmpy
          dv(3,i+k) = dv(3,i+k) - tmpz
          dv(1,i+l) = dv(1,i+l) + tmpx
          dv(2,i+l) = dv(2,i+l) + tmpy
          dv(3,i+l) = dv(3,i+l) + tmpz
          end if
13      continue
        do 14 j = 1,igrid-1
          k = (j-1)*npt
          l = k + npt
          tmp1 = d2w/d0(j)
          do 14 i = 1,npt
          if (ipick(i).gt.0) then
          tmpx  = (r(1,i+k)-r(1,i+l))*tmp1
          tmpy  = (r(2,i+k)-r(2,i+l))*tmp1
          tmpz  = (r(3,i+k)-r(3,i+l))*tmp1
          dv(1,i+k)      = dv(1,i+k)    + tmpx
          dv(2,i+k)      = dv(2,i+k)    + tmpy
          dv(3,i+k)      = dv(3,i+k)    + tmpz
          dv(1,i+l) = dv(1,i+l) - tmpx
          dv(2,i+l) = dv(2,i+l) - tmpy
          dv(3,i+l) = dv(3,i+l) - tmpz
          end if
14      continue 
        if (debug1) then
        write(*,*)' structure prtc dx dy dx'
        do 15 j=1,igrid-1
         k= (j-1)*npt
         do 15 i=1,npt
         write(*,100)j,i,dv(1,i+k),dv(2,i+k),dv(3,i+k)
100             format(1x,2(i7,1x),3(f10.5,1x))
15      continue
        end if
c@@@@@@@@@@@@@@@@@
c       write(*,*)' Analytical derivatives '
c       do 20 j=1,igrid
c        k = (j-1)*npt
c        do 20 i=1,npt
c          write(*,100)dv(1,k+i), dv(2,k+i), dv(3,k+i)
c20     continue
c100       format(1x,3(f12.7,1x))
c@@@@@@@@@@@@@@@@@
        return
        end

