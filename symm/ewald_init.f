c
c This file contains several subroutines relevant for the evaluation
c of the long range interactions using PME algorithm. Some of them
c are directly taken from the authors of PME algorithm. JM X.1996
c
c--------------------------------------------------------------------
c
        subroutine ewald_init()
c
c a subroutine to set up parameters for ewald summation
c its output contains: ewaldcoef, cutd and cutr for direct and
c receiprocal sums, number of necessary receiprocal space vectors
c
c it also prepares LES particles list
c
c cutoff and a,b,c come from line_loop.f through NBLIST and SYMM resp.
c dtol comes from line_loop.f through EWALD.BLOCK
c
      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

      double precision epsilon,dtolopt
      double precision cutoff,cut2,dtolmin,dtolmax,scfac
      double precision epstmp,qi,spi,term,e_self_les
      double precision totchg,epsilon1  
      integer i,j,mxerftab,k,m,namel,level
        integer npt_elec

        character*10 name

        data name/'ewald_init'/
        data namel/10/

        if (metalyes) then
                npt_elec = 2 * npt
        else
                npt_elec = npt
        end if

      pi=4.d0*datan(1.d0)
c
      epsilon=1.d-18
      epsilon1=1.d-7
      dtolopt=5.2d-05   
      dtolmin=1.d-10
      dtolmax=0.1d0

c      cutoff=rmax
      cutoff=sqrt(cutele2)
      totchg=0.d0

        do 33 i=1,npt
           totchg = totchg + ptchg(i)
 33     continue        

        if (abs(totchg).gt.epsilon1 
     6          .and. (.not.metalyes)) then
           write(*,*)'The charge of the system is',totchg
           level = 1
           call alert(name,namel,
     *          'Add counterions and try again',29,level)
        endif


      write(stdo,*)' -------- PME - setting up --------'
c
c adjusting dtol for optimal choice of ewaldcof
      if (dtol.lt.dtolmin) then
         dtol = dtolmin
      else if (dtol.gt.dtolmax) then 
         dtol = dtolopt
      end if

c For this value of cutoff and direct space tolerance dtol, first find
c the ewald coefficient ewaldcof

      call find_ewaldcof(cutoff,dtol,ewaldcof)

c Next find direct space cutoff for "exact" direct sum with given value
c of ewaldcof - this is only for referencing yourself. 

      call find_cutoff(cut2,epsilon,ewaldcof)
c
      write(stdo,11) ewaldcof,cutoff
 11   format('  ewald coeff. = ',f8.6,'   actuall cutoff = ',f6.2)      
      write(stdo,13) dtol
      write(stdo,14) cut2
 13     format('  contrib. to direct sum on cutoff sphere= ',e12.4e2)   
 14     format('  desired direct space cutoff for "exact" = ',f8.2)     

c prepare table of values of erfc
c note - the maximum distance for nonbond interaction is assumed to be
c        less than 1.5*cutoff; scfac = 1.5d0 
c        erftbdns is number of grid points per unit int.

c high accuracy
      erftbdns = 2.0d0*perft
      erftbdns = 0.5d0*perft
      scfac = 1.5d0     
      mxerftab = int(ewaldcof*erftbdns*cutoff*scfac)
      if ( mxerftab .gt. MAXTAU )then
         level = 1
         call alert(name,namel,
     *          'erfc needs more room - increase MAXTAU',38,level)
      endif
c@
c
      write(*,*)' before 115 mxerftab = ',mxerftab
      do 115 i=1,mxerftab
        tau(i)=0.d0
        erf_arr(1,i)=0.d0
        erf_arr(2,i)=0.d0
        erf_arr(3,i)=0.d0
        erf_arr(4,i)=0.d0
 115  continue   
      call fill_erf_table(erftbdns,mxerftab,erf_arr,tau)

c
c Calculate self-energy
c
      if (my_pe .eq. 0) then
      spi=1.0d0/pi
      spi=dsqrt(spi)
      term=2.0d0*spi*ewaldcof
      epstmp = kofdie/eps
      e_self = 0.0d0
      do 117 i=1,npt    
c       qi  = ptchg(i)*epstmp
        qi = ptchg(i) 
        e_self = e_self + qi*qi
 117  continue  
c      e_self = - ewaldcof*spi*e_self
      e_self = - ewaldcof*spi*epstmp*e_self
        if (metalyes) e_self = e_self + e_self
      e_self_les = 0.d0 
      tmp_self = e_self
      else
        tmp_self = 0.d0
        e_self = 0.d0
      end if
        
c
c Next setup the reciprocal space calculations.
c

      reclng(1) = a
      reclng(2) = b
        if (metalyes) reclng(2) = reclng(2) + reclng(2)
      reclng(3) = c
      do 133 i = 1,3
       do 132 j = 1,3
        recip(i,j) = 0.d0
132    continue
       recip(i,i) = 1.d0/reclng(i)
133   continue

      if (nfft1.lt.1 .or. nfft2.lt.1 .or. nfft3.lt.1) then
         call adjust_grid(reclng(1),reclng(2),reclng(3),
     $          nfft1,nfft2,nfft3,sgridx,sgridy,sgridz)
      end if     

       call pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,npt_elec,intrpord,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      if ( siz_Q .gt. MAXT )then
        level = 1
        call alert(name,namel,
     *          'PME needs more room - increase MAXT',35,level)
      endif
c
      if ( intrpord .gt. MAXORD )then 
        level = 1
        call alert(name,namel,
     * 'max grid interpolation order MAXORD too small',45,level)
      endif
c
      call pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,intrpord,sizfftab,sizffwrk)
c 
      write(stdo,15) nfft1,nfft2,nfft3,intrpord-1
 15   format('  Grid dimensions in PME: x = ',i3,' y = ',i3,
     &       ' z = ',i3,' ; intrpord = ',i2)
c 
      e_ew_receip = 0.d0
c
      write(stdo,*)' -------- PME - end of setting up --------'
c
c Prepare now LES particles list
c
      k=0       
      numles=0  
      if (lestyp.gt.0) then
         lesflag=.true.
         do 20 i=1,npt
            if (lesid(i).gt.0) then
               k=k+1
               if (k.gt.maxlespt) then
                level = 1
                call alert(name,namel,'MAXLESPT too small',18,level)
               end if     
               lespt(k)=i
c              now reorder lespt
               do 17 j=1,k-1
                  m=k-j
                  if (i.lt.lespt(m)) then
                     lespt(m+1)=lespt(m)
                     lespt(m)=i
                  end if   
 17            continue   
            end if   
 20      continue   
         numles=k

         do 25 i=1,numles
            k=lespt(i)
            lesqm(i,1)=ptchg(k)
            lesqm(i,3)=ptms(k)
            e_self_les = e_self_les + ptchg(k)*ptchg(k)
 25      continue       
         e_self_les = - ewaldcof*spi*epstmp*e_self_les
         e_self = e_self - e_self_les   
         tmp_self = e_self
      end if    

c Fix the size of the LES hole: if two LES prts occupy the same position
c they are replaced by one particle in receip
c
      leshole=0.1d0     
c
      return
      end
c
c------------------------------------------------------------------
      subroutine adjust_grid(a,b,c,nfft1,nfft2,nfft3,
     $                       sgridx,sgridy,sgridz)
      implicit none
      double precision a,b,c
      double precision sgridx,sgridy,sgridz
      integer i,j,k,n,i2,i3,nfft1,nfft2,nfft3
c
c note - it is desirable to have nnft at least of the nbox
c dimension and moreover, powers of 2,3,5 (the latter not
c supported here) are desirable for efficiency of fft;
c this crude routine may not give you the optimal choice, but
c it does not require explicit decomposition of a,b,c in terms
c of powers of 2 and 3
c
c by sgx, -y, -z one may change this, but try first sgx=0.5,2,4 aso
c

      i=int(log(a))
      j=int(log(b))
      k=int(log(c))

      i2=int(2**(i+1))
      i3=int(3**(i))
      n=i2
      if (abs(i3-a).lt.abs(i2-a)) then
         n=i3
      end if
      if (i3.lt.a)  n=int(2**(i-1))*i2
      nfft1=int(n*sgridx)

      i2=int(2**(j+1))
      i3=int(3**(j))
      n=i2
      if (abs(i3-b).lt.abs(i2-b)) then
         n=i3
      end if
      if (i3.lt.b)  n=int(2**(j-1))*i2
      nfft2=int(n*sgridy)

      i2=int(2**(k+1))
      i3=int(3**(k))
      n=i2
      if (abs(i3-c).lt.abs(i2-c)) then
         n=i3
      end if
      if (i3.lt.c)  n=int(2**(k-1))*i2
      nfft3=int(n*sgridz)
c
      return
      end
c
c
c------------------------------------------------------------------
      subroutine find_ewaldcof(cutoff,dtol,ewaldcof)
      implicit none
      double precision cutoff,dtol,ewaldcof

      integer i,n
      double precision pi,term,x,xlo,xhi,y,erfc

c first get direct sum tolerance. How big must ewaldcof be to get
c terms outside the cutoff below tol
      pi = 3.14159265358979323846d0

      x = 0.5d0
      i = 0
10    x = 2.d0 * x
      i = i + 1
      y = x*cutoff
      call erfcfun(y,erfc)
      term = erfc/cutoff
      if ( term .ge. dtol)goto 10
c binary search tolerance is 2 to the -60th
      n = i + 60
      xlo = 0.d0
      xhi = x
      do 20 i = 1,n
        x = (xlo+xhi)/2
        y = x*cutoff
        call erfcfun(y,erfc)
        term = erfc/cutoff
        if ( term .ge. dtol )then
           xlo = x
        else 
           xhi = x
        endif
20    continue
      ewaldcof = x
      write(6,*)"ewald_cof:",ewaldcof
      return
      end
c----------------------------------------------
      subroutine find_cutoff(cutoff,dtol,ewaldcof)
      implicit none
      double precision cutoff,dtol,ewaldcof

      integer i,n
      double precision pi,term,x,xlo,xhi,y,erfc

c first get direct sum tolerance. How big must ewaldcof be to get
c terms outside the cutoff below tol
      pi = 3.14159265358979323846d0

      x = 0.5d0
      i = 0
10    x = 2.d0 * x
      i = i + 1
      y = x*ewaldcof
      call erfcfun(y,erfc)
      term = erfc/x
      if ( term .ge. dtol)goto 10
c binary search tolerance is 2 to the -60th
      n = i + 60
      xlo = 0.d0
      xhi = x
      do 20 i = 1,n
        x = (xlo+xhi)/2
        y = x*ewaldcof
        call erfcfun(y,erfc)
        term = erfc/x
        if ( term .ge. dtol )then
           xlo = x
        else 
           xhi = x
        endif
20    continue
      cutoff = x

      return
      end
c-------------------------------------------------------------
      subroutine pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      implicit none
      integer nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack

c INPUT  
c      nfft1,nfft2,nfft3,numatoms,order
c      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
c      numatoms is number of atoms
c      order is the order of B-spline interpolation

c OUTPUT
c      sizfftab,sizffwrk,siztheta,siz_Q
c      sizfftab is permanent 3d fft table storage
c      sizffwrk is temporary 3d fft work storage
c      siztheta is size of arrays theta1-3 dtheta1-3
c      sizheap is total size of permanent storage
c      sizstack is total size of temporary storage


c This routine computes the above output parameters needed for 
c heap or stack allocation.

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      siztheta = numatoms*order
      siz_Q = 2*nfftdim1*nfftdim2*nfftdim3
      sizheap = nfft1+nfft2+nfft3+sizfftab
      sizstack = siz_Q+6*siztheta+sizffwrk+3*numatoms
      write(6,*)' total HEAP storage needed for PME = ',sizheap
      write(6,*)' total STACK storage needed for PME= ',sizstack

      return
      end
c----------------------------------------------------
c NEXT PART IS FOR INTERPOLATION OF ERFC
c----------------------------------------------------
      subroutine fill_erf_table(erftbdns,mxerftab,
     $       erf_arr,tau)

      implicit none
      integer mxerftab
      double precision erftbdns,erf_arr(4,*),tau(*)

      double precision del,x,pi,fac,erf
      integer i,ibcbeg,ibcend

      del = 1.d0 / erftbdns
      pi = 3.14159265358979323846d0
      fac = 2.d0 / sqrt(pi)

       write(6,*)'del  erftbdns ',del,erftbdns
       write(*,*)' fac '

      erf_arr(2,1) = -fac
      x = (mxerftab-1)*del
      erf_arr(2,mxerftab) = -fac * exp(-x*x)
      do 100 i = 1,mxerftab
        x = del*(i-1)
        call erfcfun(x,erf)
        tau(i) = x
        erf_arr(1,i) = erf
c       if (i.lt.10) write(6,*)'i imax x erf ',i,mxerftab,x,erf
100   continue
      ibcbeg = 1
      ibcend = 1
      call cubspl ( tau, erf_arr, mxerftab, ibcbeg, ibcend )
      return
      end
c-------------------------------------------------------------
c the code below is from netlib
c-------------------------------------------------------
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
c  from  * a practical guide to splines *  by c. de boor    
c     ************************  input  ***************************
c     n = number of data points. assumed to be .ge. 2.
c     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
c        data points. tau is assumed to be strictly increasing.
c     ibcbeg, ibcend = boundary condition indicators, and
c     c(2,1), c(2,n) = boundary condition information. specifically,
c        ibcbeg = 0  means no boundary condition at tau(1) is given.
c           in this case, the not-a-knot condition is used, i.e. the
c           jump in the third derivative across tau(2) is forced to
c           zero, thus the first and the second cubic polynomial pieces
c           are made to coincide.)
c        ibcbeg = 1  means that the slope at tau(1) is made to equal
c           c(2,1), supplied by input.
c        ibcbeg = 2  means that the second derivative at tau(1) is
c           made to equal c(2,1), supplied by input.
c        ibcend = 0, 1, or 2 has analogous meaning concerning the
c           boundary condition at tau(n), with the additional infor-
c           mation taken from c(2,n).
c     ***********************  output  **************************
c     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
c        of the cubic interpolating spline with interior knots (or
c        joints) tau(2), ..., tau(n-1). precisely, in the interval
c        (tau(i), tau(i+1)), the spline f is given by
c           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
c        where h = x - tau(i). the function program *ppvalu* may be
c        used to evaluate f or its derivatives from tau,c, l = n-1,
c        and k=4.
      implicit none
      integer ibcbeg,ibcend,n,   i,j,l,m
      double precision c(4,n),tau(n),   divdf1,divdf3,dtau,g
c****** a tridiagonal linear system for the unknown slopes s(i) of
c  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(i) ending up in c(2,i), all i.
c     c(3,.) and c(4,.) are used initially for temporary storage.
      l = n-1
compute first differences of tau sequence and store in c(3,.). also,
compute first divided difference of data and store in c(4,.).
      do 10 m=2,n
         c(3,m) = tau(m) - tau(m-1)
   10    c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
construct first equation from the boundary condition, of the form
c             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     go to 12
c     no condition at left end and n = 2.
      c(4,1) = 1.d0
      c(3,1) = 1.d0
      c(2,1) = 2.d0*c(4,2)
                                        go to 25
c     not-a-knot condition at left end and n .gt. 2.
   12 c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.d0*c(3,1))*c(4,2)*c(3,3)+
     $          c(3,2)**2*c(4,3))/c(3,1)
                                        go to 19
c     slope prescribed at left end.
   15 c(4,1) = 1.d0
      c(3,1) = 0.d0
                                        go to 18
c     second derivative prescribed at left end.
   16 c(4,1) = 2.d0
      c(3,1) = 1.d0
      c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
   18 if(n .eq. 2)                      go to 25
c  if there are interior knots, generate the corresp. equations and car-
c  ry out the forward pass of gauss elimination, after which the m-th
c  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   19 do 20 m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.d0 * 
     $           (c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
   20    c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
construct last equation from the second boundary condition, of the form
c           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since c array happens to be set up just right for it
c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
c     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
c     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.d0*g)*c(4,n)*c(3,n-1)
     *            + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
c     knot at left end point).
   22 c(2,n) = 2.d0*c(4,n)
      c(4,n) = 1.d0
                                        go to 28
c     second derivative prescribed at right endpoint.
   24 c(2,n) = 3.d0*c(4,n) + c(3,n)/2.d0*c(2,n)
      c(4,n) = 2.d0
                                        go to 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                go to 22
c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
   28 g = -1.d0/c(4,n-1)
complete forward pass of gauss elimination.
   29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
carry out back substitution
   30 j = l 
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
c****** generate cubic coefficients in each interval, i.e., the deriv.s
c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
         c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
   50    c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
                                        return
      end

c----------------------------------------------------
c END OF PART FOR INTERPOLATION OF ERFC
c----------------------------------------------------
