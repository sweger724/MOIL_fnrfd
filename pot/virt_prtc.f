c this file contains several subroutines relevant for calculations using
c virtual particles   -  05.97 jm
c
c since finally vp_list and vp_ex should be created on the con level all
c the structures are defined in CONNECT.BLOCK
c
        subroutine vp_init()

c
c initializing necessary structures (virt. prts exclusion list) and seting
c up the initial coordinates
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'

        integer i,k,kprev
        character*3 ptname

c
c prepare vp_list, vp_ex and vp_point
c
c for a given virt. prtc k included in vp_list, vp_ex contains numbers of
c real particles related to it (assumed to be just before VP in CRD file)
c
c as these related real particles are having numbers k-1, k-2 ...
c it is actually obsolate now - however, assuming that in the future one may
c move this initialization to con level, it is better to use explicitly vp_ex
c that may be in the future prepared by the connectivity module
c
        do 20 i=1,nvptmax
           vp_list(i)=0
 20     continue   
        nvpt=0

        k=0
        do 100 i=1,npt
           ptname = ptnm(i)(1:3)
           if (ptname(1:2).eq.'VP') then
              nvpt = nvpt + 1
              vp_list(nvpt) = i
              kprev = k
c temporarily VP is assumed to be VP@
              if ((ptname.eq.'VP2').or.(ptname.eq.'VP ')) then
                 k = k+2
                 vp_ex(k-1) = i-2
                 vp_ex(k) = i-1
              end if  
              if (ptname.eq.'VP3') then
                 k = k+3
                 vp_ex(k-2) = i-3
                 vp_ex(k-1) = i-2
                 vp_ex(k) = i-1
              end if
              if (k.eq.kprev) then
                 write(6,*)' Unidentified VP detected : ',ptname
                 stop
              end if  
              vp_point(nvpt) = k  
           end if   
           if (nvpt.gt.nvptmax) then
              write(6,*)' nvptmax too small - enlarge it'
              stop
           end if   
 100    continue   

        call vp_locate()

        return
        end
c-----------------------------------------------------------------
        subroutine vp_locate()

c
c get the virtual particles coordinates from the related real particles coord.
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'


        integer i,j,istart,iend,n,ii
        double precision denom,term1,term2,term3

        istart = 1
        do 100 i=1,nvpt
           ii = vp_list(i)
           iend = vp_point(i)
           denom = 0.d0
           do 50 j=istart,iend
              n = vp_ex(j)
              if(com_flag) then
                 denom = denom + ptms(n)
              end if  
              if(gcnt_flag) then
                 denom = denom + 1.d0
              end if  
 50        continue
           denom = 1.d0/denom
           term1 = 0.d0
           term2 = 0.d0
           term3 = 0.d0
           do 75 j=istart,iend
              n = vp_ex(j)
              if(com_flag) then
                 term1 = term1 + ptms(n)*coor(1,n)
                 term2 = term2 + ptms(n)*coor(2,n)
                 term3 = term3 + ptms(n)*coor(3,n)
              end if  
              if(gcnt_flag) then
                 term1 = term1 + coor(1,n)
                 term2 = term2 + coor(2,n)
                 term3 = term3 + coor(3,n)
              end if  
 75        continue
           coor(1,ii) = denom*term1
           coor(2,ii) = denom*term2
           coor(3,ii) = denom*term3
           istart = iend + 1   
 100    continue   


        return
        end
c-----------------------------------------------------------------
        subroutine vp_fdistrib()

c
c distribute the forces acting on virtual particles among the related real 
c particles according to appropriate for a given case (com, gcnt ...) scheme
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/RESTART.BLOCK'


        integer i,j,istart,iend,n,ii
        double precision denom,term1,term2,term3




        if(gcnt_flag) then
          write(6,*)' I am sorry but you must first add force 
     &                redistribution scheme in vp_fdistrib '
          stop
        end if 


        istart = 1
        do 100 i=1,nvpt
           ii = vp_list(i)
           iend = vp_point(i)
           denom = 0.d0
           do 50 j=istart,iend
              n = vp_ex(j)
              if(com_flag) then
                 denom = denom + ptms(n)
              end if  
 50        continue
           denom = 1.d0/denom
           term1 = 0.d0
           term2 = 0.d0
           term3 = 0.d0
           do 75 j=istart,iend
              n = vp_ex(j)
              if(com_flag) then
                 term1 = denom*ptms(n)*dpot(1,ii)
                 term2 = denom*ptms(n)*dpot(2,ii)
                 term3 = denom*ptms(n)*dpot(3,ii)
              end if  
              dpot(1,n) = dpot(1,n) + term1
              dpot(2,n) = dpot(2,n) + term2
              dpot(3,n) = dpot(3,n) + term3
 75        continue
           dpot(1,ii) = 0.d0
           dpot(2,ii) = 0.d0
           dpot(3,ii) = 0.d0
           istart = iend + 1   
 100    continue   


        return
        end

c-----------------------------------------------------------------
        subroutine vp_cdie()

c
c remove electrostatic contributions to energy and forces coming from
c the presence of the VP among the neighbours of related to it real particles
c i.e. remove what is undesirably added in cdie.f 
c
c if A is related to VP and B is chemically bonded to A, one may wish to
c remove also elctrostatic iteraction VP-B : use NBLIST for this purpose ...
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
c       include 'COMMON/VECTOR.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/SPECL.BLOCK'


        integer i,j,k,istart,iend,n
        double precision epstmp
        double precision rx,ry,rz,r2,s2,e2,q,df2
        double precision s,tmp
        double precision qi
        double precision xi,yi,zi,dxi,dyi,dzi

        if (eelyes) then
                epstmp = kofdie/eps
        else
                epstmp =  0.0d0
        end if

        istart = 1
        do 100 k=1,nvpt
           i = vp_list(k)
           iend = vp_point(k)
           tmp = 1.d0/ptwei(i)
           qi  = ptchg(i)*epstmp
           xi  = coor(1,i)
           yi  = coor(2,i)
           zi  = coor(3,i)
           dxi = 0.d0
           dyi = 0.d0
           dzi = 0.d0

           do 50 n=istart,iend
              j = vp_ex(n)

              rx = xi - coor(1,j)
              ry = yi - coor(2,j)
              rz = zi - coor(3,j)
              r2=rx*rx+ry*ry+rz*rz
              s2=1.0d0/r2
              q = qi*ptchg(j)
              s = dsqrt(s2)

              if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  then
                 q = q*tmp
              end if

              e2 = - q*s
              df2 = - e2*s2
c
c NOTE - df2 = -e2*s2 = +q*s*s2 which is opposite of df2 in cdie.f
c

              rx = df2*rx
              ry = df2*ry
              rz = df2*rz
              dxi = dxi + rx
              dyi = dyi + ry
              dzi = dzi + rz
              dpot(1,j) = dpot(1,j) - rx
              dpot(2,j) = dpot(2,j) - ry
              dpot(3,j) = dpot(3,j) - rz
              e_el = e_el + e2

 50        continue

           dpot(1,i) = dpot(1,i) + dxi
           dpot(2,i) = dpot(2,i) + dyi
           dpot(3,i) = dpot(3,i) + dzi

           istart = iend + 1   

 100    continue   


        return
        end

c-----------------------------------------------------------------
        subroutine vp_cdie_ewald()
        implicit none

c
c remove electrostatic contributions to energy and forces coming from
c the presence of the VP among the neighbours of related to it real particles
c i.e. remove what is undesirably added to e_dir in cdie_ewald.f 
c
c remove also the corresponding reciprocal contributions due to VP and its
c related real particles added to e_recip
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/EWALD.BLOCK'


        integer i,j,k,istart,iend,n,ind
        double precision epstmp
        double precision rx,ry,rz,r2,s2,e2,q,df,df1,df2
        double precision s,tmp
        double precision qi,spi,rij,xerfc,erf,ovt,del,term
        double precision xi,yi,zi,dxi,dyi,dzi,derfc,dx,erfcc

        if (eelyes) then
                epstmp = kofdie/eps
        else
                epstmp =  0.0d0
        end if

c
        spi=1.0d0/pi
        spi=dsqrt(spi)
        term=2.0d0*spi*ewaldcof
        ovt=1.0d0/3.0d0
        del=1.0d0/erftbdns

        istart = 1
        do 100 k=1,nvpt
           i = vp_list(k)
           iend = vp_point(k)
           tmp = 1.d0/ptwei(i)
           qi  = ptchg(i)*epstmp
           xi  = coor(1,i)
           yi  = coor(2,i)
           zi  = coor(3,i)
           dxi = 0.d0
           dyi = 0.d0
           dzi = 0.d0

           do 50 n=istart,iend
              j = vp_ex(n)

              rx = xi - coor(1,j)
              ry = yi - coor(2,j)
              rz = zi - coor(3,j)
              r2=rx*rx+ry*ry+rz*rz
              s2=1.0d0/r2
              q = qi*ptchg(j)
              rij=dsqrt(r2)
              s  = 1.0d0/rij

              if ((lesid(i).ne.0) .and.
     *                      (lesid(i) .eq.lesid(j)))  then
                 q = q*tmp
              end if
c WARNING - VP + LES + PME has not been tested
c           
c remember also, if VP belongs to multiplied fragment each VP copy should be
c placed after related real particles of a given copy - at present the order
c of particles in CRD file matters
c

              xerfc = ewaldcof*rij
c cubic spline on erfc,derf
              ind = int(erftbdns*xerfc) + 1
              dx = xerfc - (ind-1)*del
              derfc = -erf_arr(2,ind) - dx*( erf_arr(3,ind) 
     $                + 0.5d0*dx*erf_arr(4,ind)) 
              erfcc = dx*( erf_arr(2,ind) + 0.5d0* dx*( erf_arr(3,ind)
     $                    + dx*erf_arr(4,ind)*ovt) ) + erf_arr(1,ind)
              
              erf = 1.d0 - erfcc

              e2 = - q*s*erf
              df1 = q*s*(s2*erf - s*ewaldcof*derfc)
              e_corr = e_corr + e2

c non-interpolated erfc
c                   call erfcfun(xerfc,erfc)
c                   erf = 1.d0 - erfc
c                   e2 = - q*s*erf
c                   df2 = q*s*(s2*erf - s*term*exp(-xerfc**2))

              e2 = - q*s*erfcc
              df2 = q*s*(s2*erfcc + s*ewaldcof*derfc)
              e_el = e_el + e2

c non-interpolated erfc
c                       call erfcfun(xerfc,erfc)
c                       e2 = q*s*erfc
c                       df2 = -q*s*(s2*erfc + s*term*exp(-xerfc**2))

              df = df1 + df2

              rx = df*rx
              ry = df*ry
              rz = df*rz
              dxi = dxi + rx
              dyi = dyi + ry
              dzi = dzi + rz
              dpot(1,j) = dpot(1,j) - rx
              dpot(2,j) = dpot(2,j) - ry
              dpot(3,j) = dpot(3,j) - rz

 50        continue

           dpot(1,i) = dpot(1,i) + dxi
           dpot(2,i) = dpot(2,i) + dyi
           dpot(3,i) = dpot(3,i) + dzi

           istart = iend + 1   

 100    continue   


        return
        end
