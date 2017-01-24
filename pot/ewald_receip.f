c
c This file contains several subroutines relevant for the evaluation
c of the long range interactions using PME algorithm. Some of them
c are directly taken from the authors of PME algorithm. JM X.1996
c
c--------------------------------------------------------------------
c
        subroutine ewald_receip()
c
c this subroutine computes the receiprocal sum using PME 
c algorithm by Tom Darden et. al. JCP 103, No 19, pg. 8577
c 
c
c      implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/PRESS.BLOCK'

      double precision volume,e_self_les,spi,epstmp
      integer maxblock  

      volume = reclng(1)*reclng(2)*reclng(3)
      maxblock=maxlespt
        pi = 4.d0*datan(1.d0)
      spi=1.0d0/pi
      spi=dsqrt(spi)
      epstmp = kofdie/eps

c
c call the modified Darden's subroutine for receiprocal sum 
c

      call do_pmesh_kspace(metalyes,b,
     $   npt,coor,ptchg,recip,volume,ewaldcof,
     $   intrpord,nfft1,nfft2,nfft3,e_ew_receip,dpot,lesdpot,
     $   sizfftab,sizffwrk,siztheta,siz_Q,
     $   bsp_mod1,bsp_mod2,bsp_mod3,fftable,QQ,ffwork,
     $   theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,
     $   fr1,fr2,fr3,epstmp,e_self_les,fpme,tmples,
     $   maxblock,lespoint,leshole,killer,
     $   leskill,lesqm,lesflag,numles,lespt,
     $   pressON,V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ,
     $   xrel,yrel,zrel,V_PIrelXX,V_PIrelYY,V_PIrelZZ)

c see do_pmesh_kspace for explanation of arguments


c modify self-energy due to les holes

      if (lesflag) then 
         e_self_les = - ewaldcof*spi*epstmp*e_self_les
         e_self = tmp_self + e_self_les 
c        write(6,*) 'self-hole ',e_self_les
      end if

      return
      end
c-------------------------------------------------------
      subroutine fix_dim(nfft,nfftdim)
      integer nfft,nfftdim,n
      nfftdim = nfft
      n = nfft/2
      if ( nfft .eq. 2*n )nfftdim = nfft+1
      return
      end
c----------------------------------------------------
      subroutine pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)

c      implicit none

c  see DO_PMESH_KSPACE for explanation of arguments

      integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)
      double precision fftable(sizfftab),ffwork(sizffwrk)
   
      double precision dummy
      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      call fft_setup(dummy,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      return
      end
c----------------------------------------------------
      subroutine do_pmesh_kspace(metalyes,y_step,
     $   numatoms,coor,ptchg,recip,volume,ewaldcof,
     $   order,nfft1,nfft2,nfft3,eer,dpot,lesdpot,
     $   sizfftab,sizffwrk,siztheta,siz_Q,
     $   bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q,ffwork,
     $   theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,
     $   fr1,fr2,fr3,epstmp,e_self_les,fpme,tmples,
     $   maxblock,lespoint,leshole,killer,
     $   leskill,lesqm,lesflag,numles,lespt,
     $   pressON,V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ,
     $   xrel,yrel,zrel,V_PIrelXX,V_PIrelYY,V_PIrelZZ)
c      implicit none

c INPUT 
c       numatoms:  number of atoms
c       coor: x,y,z  atomic coords
c       ptchg:  atomic charges
c       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
c       volume: the volume of the unit cell
c       ewaldcof:   ewald convergence parameter
c       order: the order of Bspline interpolation. E.g. cubic is order 4
c          fifth degree is order 6 etc. The order must be an even number 
c          and at least 4.
c       nfft1,nfft2,nfft3: the dimensions of the charge grid array
c
      integer numatoms,order,nfft1,nfft2,nfft3
      double precision coor(3,numatoms),volume,eer,dpot(3,numatoms),
     $       ptchg(numatoms),recip(3,3),ewaldcof,epstmp

c
c OUTPUT
c       eer:  ewald reciprocal or k-space  energy
c       dpot: - derivatives incremented by k-space sum

c SIZES of some arrays
      integer   sizfftab,sizffwrk,siztheta,siz_Q,maxblock

c HEAP STORAGE:  These arrays need to be preserved throughout simulation
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     $                 bsp_mod3(nfft3),fftable(sizfftab)
c STACK STORAGE: These arrays can be tossed after leaving this routine
      double precision Q(siz_Q),ffwork(sizffwrk),theta1(siztheta),
     $          theta2(siztheta),theta3(siztheta),dtheta1(siztheta),
     $          dtheta2(siztheta),dtheta3(siztheta),
     $          fr1(numatoms),fr2(numatoms),fr3(numatoms)

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
c metal flag
        logical metalyes
        double precision y_step
c LES variables
      double precision leshole,lesqm(maxblock,4),e_self_les,fpme
      double precision lesdpot(3,maxblock)
        integer tmples(maxblock)        
      integer lespoint(2*maxblock),killer(maxblock)
      logical lesflag,leskill(maxblock,2)       
      integer numles,lespt(maxblock),iblock,i   
c pressure variables
      double precision V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ
      double precision V_PIrelXX,V_PIrelYY,V_PIrelZZ
      double precision xrel(*),yrel(*),zrel(*)
      logical pressON

c  get some integer array dimensions
      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
c  prepare fractional coordinates of particles, prepare LES holes      
      call get_scaled_fractionals(metalyes,y_step,
     $         numatoms,coor,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3,
     $         e_self_les,tmples,lesdpot,maxblock,lespoint,iblock,
     $         leshole,ptchg,
     $         leskill,lesqm,lesflag,numles,lespt,killer)
c  get the coefficients of b-splimes
      call get_bspline_coeffs(metalyes,
     $         numatoms,fr1,fr2,fr3,order,maxblock,iblock,lespoint,
     $         theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
c  prepare the grid table
      call fill_charge_grid(metalyes,
     $         numatoms,ptchg,theta1,theta2,theta3,fr1,fr2,fr3,epstmp,
     $         maxblock,iblock,lespoint,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)

      call fft_back(
     $         Q,fftable,ffwork,nfft1,nfft2,nfft3,
     $         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
c  compute the receiprocal energy
      call scalar_sum(
     $         Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,
     $         pressON,V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ)
      call fft_forward(
     $         Q,fftable,ffwork,nfft1,nfft2,nfft3,
     $         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
c  compute the receiprocal contribution to derivatives
      call grad_sum(metalyes,
     $         numatoms,ptchg,recip,theta1,theta2,theta3,dpot,lesdpot,
     $         dtheta1,dtheta2,dtheta3,fr1,fr2,fr3,epstmp,lesflag,fpme,
     $         maxblock,iblock,lespoint,lespt,lesqm,killer,numles,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,
     $         xrel,yrel,zrel,V_PIrelXX,V_PIrelYY,V_PIrelZZ,pressON)


c  for LES particles get the original charges
      if ((lesflag).and.((iblock.gt.1).or.
     & ((iblock.eq.1).and.(lespoint(maxblock+1).ne.numatoms)))) then
         do 125 i=1,numles
           ptchg(lespt(i))=lesqm(i,1)
 125     continue       
      end if
c

      return
      end
c----------------------------------------------------------------------
      subroutine get_scaled_fractionals(metalyes,y_step,
     $           numatoms,coor,recip,nfft1,nfft2,nfft3,
     $           fr1,fr2,fr3,e_self_les,tmples,lescoor,
     $           maxblock,lespoint,iblock,leshole,ptchg,
     $           leskill,lesqm,lesflag,numles,lespt,killer)
c      implicit none

c INPUT:
c      numatoms: number of atoms
c      x,y,z: arrays of cartesian coords
c      recip: the 3x3 array of reciprocal vectors stored as columns
c OUTPUT:
c     fr1,fr2,fr3 the scaled and shifted fractional coords
c     pointers to killed LES particles not to be included in receip

      integer numatoms,nfft1,nfft2,nfft3,maxblock
      integer tmples(maxblock)
c     use lesdpot array for lescoor
      double precision lescoor(3,maxblock),xcm,zcm,ycm,r2,x,y,z
      double precision recip(3,3),coor(3,numatoms),ptchg(numatoms)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)

      integer lespoint(2*maxblock)
      double precision  lesqm(maxblock,4),qj,e_self_les
      integer numles,lespt(maxblock),iblock,killer(maxblock)    
      logical lesflag,leskill(maxblock,2)       

      integer n,k,i,j,lastprt,jprev,numkill,m,p
      double precision w,leshole2,rx,ry,rz,mi,mj,mij,leshole
c metal flag (NOTE THAT AT PRESENT LES DOES NOT WORK WITH METALS)
        logical metalyes
        integer n_metal
        double precision y_step


      do 100 n = 1,numatoms
        x=coor(1,n) 
        y=coor(2,n) 
        z=coor(3,n) 
        w = x*recip(1,1)+y*recip(2,1)+z*recip(3,1)
        fr1(n) = nfft1*(w - anint(w) + 0.5d0)
        w = x*recip(1,2)+y*recip(2,2)+z*recip(3,2)
        fr2(n) = nfft2*(w - anint(w) + 0.5d0)
        w = x*recip(1,3)+y*recip(2,3)+z*recip(3,3)
        fr3(n) = nfft3*(w - anint(w) + 0.5d0)
 100  continue

        if (metalyes) then
                do 101 n=1,numatoms
                 x = coor(1,n)
                 z = coor(3,n)
                 if (coor(2,n).lt.0.d0) then
                        y = -y_step - coor(2,n)
                 else
                 y =  y_step - coor(2,n)
                 end if
                 w = x*recip(1,1)+y*recip(2,1)+z*recip(3,1)
                 n_metal = n+numatoms
                 fr1(n_metal) = nfft1*(w - anint(w) + 0.5d0)
                 w = x*recip(1,2)+y*recip(2,2)+z*recip(3,2)
                 fr2(n_metal) = nfft2*(w - anint(w) + 0.5d0)
                 w = x*recip(1,3)+y*recip(2,3)+z*recip(3,3)
                 fr3(n_metal) = nfft3*(w - anint(w) + 0.5d0)
101             continue
        end if
c now check for LES particles on the top of each other
c prepare necessary pointers, recalculate relevant fr and charges

        do 120 n=1,maxblock
           tmples(n)=0
           lespoint(n)=0
           lespoint(n+maxblock)=0
 120    continue 

        do 125 i=1,numles
c        killing
         leskill(i,1)=.false.
c        being killed
         leskill(i,2)=.false.
         killer(i)=i
         lesqm(i,2)=lesqm(i,1)
         lesqm(i,4)=lesqm(i,3)
 125    continue        
  
        lespoint(1)=1
        lespoint(maxblock+1)=numatoms
        iblock=1
        e_self_les = 0.d0

        if (lesflag) then

        lastprt=numatoms   
        leshole2=leshole*leshole
        numkill=0

        do 200 n=1,numles

           if (leskill(n,2)) go to 200
           i=lespt(n)
           lescoor(1,n)=coor(1,i)
           lescoor(2,n)=coor(2,i)
           lescoor(3,n)=coor(3,i)
           x  = lescoor(1,n)
           y  = lescoor(2,n)
           z  = lescoor(3,n)
           mi = lesqm(n,4)

           do 190 k=1,numles
              if (leskill(k,2)) go to 190
              j=lespt(k)
              if (.not.(i.lt.j)) go to 190
              rx = x - coor(1,j)
              ry = y - coor(2,j)
              rz = z - coor(3,j)
              mj = lesqm(k,4)
              r2=rx*rx+ry*ry+rz*rz
              if (r2.lt.leshole2) then
                 leskill(n,1)=.true.
                 leskill(k,2)=.true.
                 killer(k)=n
c                collect charges and masses
                 qj=ptchg(j)
                 lesqm(n,2)=lesqm(n,2)+qj
                 mij=mi+mj
                 lesqm(n,4)=mij
                 mij=1.d0/mij
                 xcm=mij*(mi*x+mj*coor(1,j))
                 ycm=mij*(mi*y+mj*coor(2,j))
                 zcm=mij*(mi*z+mj*coor(3,j))
                 lescoor(1,n)=xcm
                 lescoor(2,n)=ycm
                 lescoor(3,n)=zcm
c                put killed prts into tmples and reorder it
                 numkill=numkill+1
                 tmples(numkill)=j
                 do 170 p=1,numkill-1
                   m=numkill-p
                   if (j.lt.tmples(m)) then
                     tmples(m+1)=tmples(m)
                     tmples(m)=j
                   end if   
 170             continue        
              end if     
c             end of if leshole

 190       continue

           if ((.not.(leskill(n,2))).and.(leskill(n,1))) then
              x=lescoor(1,n) 
              y=lescoor(2,n) 
              z=lescoor(3,n) 
              w = x*recip(1,1)+y*recip(2,1)+z*recip(3,1)
              fr1(i) = nfft1*(w - anint(w) + 0.5d0)
              w = x*recip(1,2)+y*recip(2,2)+z*recip(3,2)
              fr2(i) = nfft2*(w - anint(w) + 0.5d0)
              w = x*recip(1,3)+y*recip(2,3)+z*recip(3,3)
              fr3(i) = nfft3*(w - anint(w) + 0.5d0)
              ptchg(i)=lesqm(n,2)
           end if  

           e_self_les = e_self_les + ptchg(i)*ptchg(i)
 
 200    continue  

c       now prepare les pointer 
        if (numkill.gt.0) then
           jprev=0
           do 300 i=1,numkill
              j=tmples(i)
c             write(6,*)'killed .. ',j
              if (j.gt.jprev+1) then
                    lespoint(maxblock+iblock)=j-1
                    iblock=iblock+1
                    lespoint(iblock)=j+1
              else
                    lespoint(iblock)=j+1                    
              end if   
              if (j.eq.lastprt) then
                    iblock=iblock-1
                    if (j.gt.jprev+1) lastprt=j-1
                    go to 300
              end if   
              lespoint(maxblock+iblock)=lastprt
              jprev=j
 300       continue  
        end if   

c       do 400 i=1,iblock
c          write(6,410)lespoint(i),lespoint(maxblock+i)
c 400   continue 
c 410   format(1x,'lespoint start - end : ',2i6)

        end if
c       with respect to lesflag

      return
      end
c---------------------------------------------------------------
      subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
c      implicit none
      integer nfft1,nfft2,nfft3,order
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)

      integer MAXORDER
      parameter (MAXORDER=25)
      integer MAXNFFT
      parameter (MAXNFFT=4000)
      double precision array(MAXORDER),darray(MAXORDER),w
      double precision bsp_arr(MAXNFFT)
      integer i,maxn

c this routine loads the moduli of the inverse DFT of the B splines
c bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
c Order is the order of the B spline approx.

      if ( order .gt. MAXORDER )then
       write(6,*)'order too large - check on MAXORDER'
       stop
      endif
      maxn = max(nfft1,nfft2,nfft3)
      if ( maxn .gt. MAXNFFT )then 
       write(6,*)'nfft1-3 too large - check on MAXNFFT'
       stop
      endif
      w = 0.d0
      call fill_bspline(w,order,array,darray)
      do 100 i = 1,maxn
        bsp_arr(i) = 0.d0
100   continue
      do 150 i = 2,order+1
       bsp_arr(i) = array(i-1)
150   continue
      call DFTMOD(bsp_mod1,bsp_arr,nfft1)
      call DFTMOD(bsp_mod2,bsp_arr,nfft2)
      call DFTMOD(bsp_mod3,bsp_arr,nfft3)
      return
      end
c------------------------------------------------------------------------
      subroutine DFTMOD(bsp_mod,bsp_arr,nfft)
c      implicit none
      integer nfft
      double precision bsp_mod(nfft),bsp_arr(nfft)
c Computes the modulus of the discrete fourier transform of bsp_arr,
c  storing it into bsp_mod

      integer j,k
      double precision sum1,sum2,twopi,arg,tiny
      twopi = 2.d0*3.14159265358979323846d0
      tiny = 1.d-7
      do 300 k = 1,nfft
       sum1 = 0.d0
       sum2 = 0.d0
       do 250 j = 1,nfft
         arg = twopi*(k-1)*(j-1)/nfft
         sum1 = sum1 + bsp_arr(j)*cos(arg)
         sum2 = sum2 + bsp_arr(j)*sin(arg)
250    continue
       bsp_mod(k) = sum1**2 + sum2**2
300   continue
      do 400 k = 1,nfft
       if ( bsp_mod(k) .lt. tiny )
     $     bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
400   continue
      return
      end
c------------------------------------------------------------------------
      subroutine get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
c      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $       nfftable,nffwork,sizfftab,sizffwrk
      integer n,nfftmax

      nfftmax = max(nfft1,nfft2,nfft3)
      nfftdim1 = nfft1
      n = nfft1/2
      if ( nfft1 .eq. 2*n )nfftdim1 = nfft1+1
      nfftdim2 = nfft2
      n = nfft2/2
      if ( nfft2 .eq. 2*n )nfftdim2 = nfft2+1
      nfftdim3 = nfft3
      n = nfft3/2
      if ( nfft3 .eq. 2*n )nfftdim3 = nfft3+1

      nfftable = 4*nfftmax + 15
      nffwork = nfftmax
      sizfftab = 3*nfftable
      sizffwrk  = 2*nfftmax

c#ifdef SGIFFT
c      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
c      nffwork = 0
c      sizfftab = nfftable
c      sizffwrk  = nffwork
c#endif
c#ifdef PUBFFT
c      nfftable = 4*nfftmax + 15
c      nffwork = nfftmax
c      sizfftab = 3*nfftable
c      sizffwrk  = 2*nfftmax
c#endif
      return
      end
c---------------------------------------------------------------
      subroutine fft_setup(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
c      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork

c      write(6,*)'using public domain fft code'
      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)


c#ifdef SGIFFT
c      call ZFFT3DI(nfft1,nfft2,nfft3,fftable)
c#endif
c#ifdef PUBFFT
c      write(6,*)'using public domain fft code'
c      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)
c#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_forward(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
c      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

      integer isign
      integer nfftable,nffwork

      isign = 1

      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)


c#ifdef SGIFFT
c      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
c     $   nfftdim1,nfftdim2,fftable)
c#endif
c#ifdef PUBFFT
c      call pubz3d(isign,nfft1,nfft2,nfft3,array,
c     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
c#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_back(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
c      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork

      integer isign

      isign = -1

      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

c#ifdef SGIFFT
c      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
c     $   nfftdim1,nfftdim2,fftable)
c#endif
c#ifdef PUBFFT
c      call pubz3d(isign,nfft1,nfft2,nfft3,array,
c     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
c#endif
      return
      end
c-----------------------------------------------------------
      subroutine fill_charge_grid(metalyes,
     $         numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,epstmp,
     $         maxblock,iblock,lespoint,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
c---------------------------------------------------------------------
c INPUT:
c      numatoms:  number of atoms
c      charge: the array of atomic charges
c      theta1,theta2,theta3: the spline coeff arrays
c      fr1,fr2,fr3 the scaled and shifted fractional coords
c      nfft1,nfft2,nfft3: the charge grid dimensions
c      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
c      order: the order of spline interpolation
c OUTPUT:
c      Q the charge grid
c---------------------------------------------------------------------
c      implicit none

      integer numatoms,order,nfft1,nfft2,nfft3
      integer nfftdim1,nfftdim2,nfftdim3
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),charge(numatoms)
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)

      integer n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision prod,epstmp,fc,fc1,prod1
      integer maxblock,iblock,lespoint(2*maxblock),kk   

c metal flag
        logical metalyes

      ntot = 2*nfftdim1*nfftdim2*nfftdim3
      call clearQ(Q,ntot)

      fc = dsqrt(epstmp)        

      do 777 kk=1,iblock        
         if (.not.metalyes) then
        do 300 n=lespoint(kk),lespoint(kk+maxblock)     
c       do 300 n = 1,numatoms
         k0 = int(fr3(n)) - order
         fc1 = charge(n)*fc
         do 200 ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          prod = theta3(ith3,n)*fc1
          do 150 ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod1 = theta2(ith2,n)*prod
           i0 = int(fr1(n)) - order
           do 100 ith1 = 1,order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod1
 100       continue
 150      continue
 200     continue
 300    continue
         else
          do 310 n=1,numatoms
         k0 = int(fr3(n)) - order
         fc1 = charge(n)*fc
         do 210 ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          prod = theta3(ith3,n)*fc1
          do 151 ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod1 = theta2(ith2,n)*prod
           i0 = int(fr1(n)) - order
           do 110 ith1 = 1,order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod1
 110       continue
 151      continue
 210     continue
 310      continue
c
c here come the image charges on the metal.
c
        do 320 n=numatoms+1,2*numatoms
         k0 = int(fr3(n)) - order
         fc1 = -charge(n-numatoms)*fc
         do 220 ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          prod = theta3(ith3,n)*fc1
          do 152 ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod1 = theta2(ith2,n)*prod
           i0 = int(fr1(n)) - order
           do 120 ith1 = 1,order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod1
 120       continue
 152      continue
 220     continue
 320    continue
         end if
 777  continue
c       
      return
      end
c-----------------------------------------------------------
      subroutine clearQ(Q,ntot)
      integer ntot
      double precision Q(ntot)
      integer i
      do 10 i = 1,ntot
        Q(i) = 0.d0
10    continue
      return
      end
c-----------------------------------------------------------
      subroutine scalar_sum(
     $         Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,
     $         pressON,V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ)
c      implicit none

      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume
      double precision recip(3,3),eer,pivol,nom
      double precision r11,r12,r13,r21,r22,r23,r31,r32,r33

      double precision pi,fac,denom,eterm,energy
      integer k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3
      double precision mhat1,mhat2,mhat3,msq,struc2
 
c pressure
      double precision V_PIrec,V_PIrecXX,V_PIrecYY,V_PIrecZZ
      double precision XX,YY,ZZ
      logical pressON

      indtop = nfft1*nfft2*nfft3
      pi = 3.14159265358979323846d0
      fac = pi**2/ewaldcof**2
      pivol = pi*volume

      r11=recip(1,1)    
      r12=recip(1,2)    
      r13=recip(1,3)    
      r21=recip(2,1)    
      r22=recip(2,2)    
      r23=recip(2,3)    
      r31=recip(3,1)    
      r32=recip(3,2)    
      r33=recip(3,3)    

      nff = nfft1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      energy = 0.d0

      do 100 ind = 1,indtop-1

c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1

       k3 = ind/nff + 1
       jnd = ind - (k3-1)*nff
       k2 = jnd/nfft1 + 1
       k1 = jnd - (k2-1)*nfft1 +1
       m1 = k1 - 1
       if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
       m2 = k2 - 1
       if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
       m3 = k3 - 1
       if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
       mhat1 = r11*m1+r12*m2+r13*m3
       mhat2 = r21*m1+r22*m2+r23*m3
       mhat3 = r31*m1+r32*m2+r33*m3
       msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
       denom = pivol*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
       nom = 1.d0/denom
       eterm = exp(-fac*msq)
       eterm = nom * eterm
       struc2 = Q(1,k1,k2,k3)*Q(1,k1,k2,k3)
       struc2 = struc2 + Q(2,k1,k2,k3)*Q(2,k1,k2,k3)
       energy = energy + eterm * struc2
       Q(1,k1,k2,k3) = eterm * Q(1,k1,k2,k3)
       Q(2,k1,k2,k3) = eterm * Q(2,k1,k2,k3)
c pressure
       if (pressON) then
        V_PIrec = 0.5d0*eterm*struc2
        XX = 1.d0-2.d0*(1.d0+fac*msq)/msq*mhat1*mhat1
        YY = 1.d0-2.d0*(1.d0+fac*msq)/msq*mhat2*mhat2
        ZZ = 1.d0-2.d0*(1.d0+fac*msq)/msq*mhat3*mhat3
        V_PIrecXX = V_PIrecXX + V_PIrec*XX
        V_PIrecYY = V_PIrecYY + V_PIrec*YY
        V_PIrecZZ = V_PIrecZZ + V_PIrec*ZZ
       endif
100   continue
      eer = 0.5d0 * energy

      return
      end
c
c-----------------------------------------------------------
      subroutine grad_sum(metalyes,
     $         numatoms,charge,recip,theta1,theta2,theta3,dpot,lesdpot,
     $         dtheta1,dtheta2,dtheta3,fr1,fr2,fr3,epstmp,lesflag,fpme,
     $         maxblock,iblock,lespoint,lespt,lesqm,killer,numles,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,
     $         xrel,yrel,zrel,V_PIrelXX,V_PIrelYY,V_PIrelZZ,pressON)
c      implicit none

      integer numatoms,order,nfft1,nfft2,nfft3
      integer nfftdim1,nfftdim2,nfftdim3
      double precision recip(3,3),fc,epstmp,term1,t1a,t2a,t3a,fpme
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),charge(numatoms)
      double precision dtheta1(order,numatoms),dtheta2(order,numatoms),
     $     dtheta3(order,numatoms)
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)

      integer n,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision f1,f2,f3,term,fxn,fyn,fzn,dpot(3,numatoms),w
      double precision r11,r12,r13,r21,r22,r23,r31,r32,r33,t1,t2,t3
      integer maxblock,iblock,lespoint(2*maxblock),kk,
     $     killer(maxblock),numles,lespt(maxblock)

      double precision lesdpot(3,maxblock),lesqm(maxblock,4),zero,qj
      logical lesflag   
c metal flag
        logical metalyes
c pressure
      double precision xrel(*),yrel(*),zrel(*)
      double precision V_PIrelXX,V_PIrelYY,V_PIrelZZ
      logical pressON

      fc = dsqrt(epstmp)
      zero = 1.d-6      

c use lesdpot to save values of dpot for LES prts
      if ((lesflag).and.((iblock.gt.1).or.
     & ((iblock.eq.1).and.(lespoint(maxblock+1).ne.numatoms)))) then
       do 10 i=1,numles
         j=lespt(i)
         lesdpot(1,i)=dpot(1,j)
         lesdpot(2,i)=dpot(2,j)
         lesdpot(3,i)=dpot(3,j)
         qj=abs(charge(j))
         if (qj.lt.zero) then
            charge(j)=1.d0
            lesqm(i,2)=charge(j)
         end if   
 10    continue  
      end if
        
      r11=recip(1,1)    
      r12=recip(1,2)    
      r13=recip(1,3)    
      r21=recip(2,1)    
      r22=recip(2,2)    
      r23=recip(2,3)    
      r31=recip(3,1)    
      r32=recip(3,2)    
      r33=recip(3,3)    

      fpme=0.d0 

C$DOACROSS LOCAL(f1,f2,f3,k0,k,j0,j,i0,i,term,n,ith1,ith2,ith3),
C$&  SHARE(numatoms,fr1,fr2,fr3,charge,Q,fx,fy,fz,recip,order,
C$&   nfft1,nfft2,nfft3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
c Summation is done only on real particles
c Even for metal.
c (no need to compute forces on image particles)
c
      do 777 kk=1,iblock        
       do 400 n=lespoint(kk),lespoint(kk+maxblock)      
c      do 400 n = 1,numatoms
        f1 = 0.d0
        f2 = 0.d0
        f3 = 0.d0
        k0 = int(fr3(n)) - order
        term1 = fc*charge(n)
        do 200 ith3 = 1,order
         k0 = k0 + 1
         k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
         j0 = int(fr2(n)) - order
         t3a = nfft3 * dtheta3(ith3,n)
         t2a = nfft2 * theta3(ith3,n)
         t1a = nfft1 * theta3(ith3,n)
         do 150 ith2 = 1,order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          i0 = int(fr1(n)) - order
          t3 = t3a * theta2(ith2,n)
          t2 = t2a * dtheta2(ith2,n)
          t1 = t1a * theta2(ith2,n)
          do 100 ith1 = 1,order
           i0 = i0 + 1
           i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
           term = term1*Q(1,i,j,k)
c force is -grad
           f1 = f1 - term * dtheta1(ith1,n) * t1
           f2 = f2 - term * theta1(ith1,n) * t2
           f3 = f3 - term * theta1(ith1,n) * t3
100       continue
150      continue
200     continue
        fxn = r11*f1+r12*f2+r13*f3
        fyn = r21*f1+r22*f2+r23*f3
        fzn = r31*f1+r32*f2+r33*f3
c
c@
c
        dpot(1,n) = dpot(1,n) - fxn
        dpot(2,n) = dpot(2,n) - fyn
        dpot(3,n) = dpot(3,n) - fzn
c pressure
        if (pressON) then
         V_PIrelXX = V_PIrelXX + fxn*xrel(n)
         V_PIrelYY = V_PIrelYY + fyn*yrel(n)
         V_PIrelZZ = V_PIrelZZ + fzn*zrel(n)
        endif
        fpme = fpme + fxn*fxn + fyn*fyn + fzn*fzn
 400   continue
 777  continue  

c now redistribute receip. contr. to derivatives among killed LES copies

      if ((lesflag).and.((iblock.gt.1).or.
     & ((iblock.eq.1).and.(lespoint(maxblock+1).ne.numatoms)))) then
       do 800 i=1,numles
         j=lespt(i)
         lesdpot(1,i)=dpot(1,j)-lesdpot(1,i)
         dpot(1,j)=dpot(1,j)-lesdpot(1,i)
         lesdpot(2,i)=dpot(2,j)-lesdpot(2,i)
         dpot(2,j)=dpot(2,j)-lesdpot(2,i)
         lesdpot(3,i)=dpot(3,j)-lesdpot(3,i)
         dpot(3,j)=dpot(3,j)-lesdpot(3,i)
 800   continue  
       do 850 i=1,numles
         j=lespt(i)
         k=killer(i)
         w=1.d0/lesqm(k,2)
         w=lesqm(i,1)*w
         dpot(1,j)= dpot(1,j) + w*lesdpot(1,k)
         dpot(2,j)= dpot(2,j) + w*lesdpot(2,k)
         dpot(3,j)= dpot(3,j) + w*lesdpot(3,k)
         fpme=fpme-lesdpot(1,k)*lesdpot(1,k)-
     &             lesdpot(2,k)*lesdpot(2,k)-lesdpot(2,k)*lesdpot(3,k)
         fpme=fpme+dpot(1,j)*dpot(1,j)+dpot(2,j)*dpot(2,j)+
     &             dpot(3,j)*dpot(3,j)
 850   continue  

      end if    

      fpme = dsqrt( fpme/(3.d0*numatoms))       
c
      return
      end
c-----------------------------------------------------------
c---------------------------------------------------------------------
      subroutine get_bspline_coeffs(metalyes,
     $           numatoms,fr1,fr2,fr3,order,maxblock,iblock,lespoint,
     $           theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
c---------------------------------------------------------------------
c INPUT:
c      numatoms: number of atoms
c      fr1,fr2,fr3 the scaled and shifted fractional coords
c      order: the order of spline interpolation
c OUTPUT
c      theta1,theta2,theta3: the spline coeff arrays
c      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
c---------------------------------------------------------------------
c      implicit none
      integer numatoms,order,maxblock,iblock,lespoint(2*maxblock)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),dtheta1(order,numatoms),
     $     dtheta2(order,numatoms),dtheta3(order,numatoms)

      double precision w
      integer n, kk
c metal flag
        logical metalyes

      do 777 kk=1,iblock        
         if (.not.metalyes) then
       do 100 n=lespoint(kk),lespoint(kk+maxblock)       
c      do 100 n = 1,numatoms
        w = fr1(n)-int(fr1(n))
        call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
        w = fr2(n)-int(fr2(n))
        call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
        w = fr3(n)-int(fr3(n))
        call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
 100   continue
         else
          do 101 n=1,2*numatoms
         w = fr1(n)-int(fr1(n))
         call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
         w = fr2(n)-int(fr2(n))
         call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
         w = fr3(n)-int(fr3(n))
         call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
101       continue
         end if
 777  continue
      return
      end
c---------------------------------------------------
      subroutine fill_bspline(w,order,array,darray)
c---------- use standard B-spline recursions: see doc file
c      implicit none
      integer order
      double precision w,array(order),darray(order)

      integer k
c do linear case
      call init(array,w,order)
c compute standard b-spline recursion
      do 10 k = 3,order-1
       call one_pass(array,w,k)
10    continue
c perform standard b-spline differentiation
      call diff(array,darray,order)
c one more recursion
      call one_pass(array,w,order)
      return
      end
c---------------------------------------------------
      subroutine init(c,x,order)
c      implicit none
      integer order
      double precision c(order),x
      c(order) = 0.d0
      c(2) = x
      c(1) = 1.d0 - x
      return
      end
c-------------------------------------
      subroutine one_pass(c,x,k)
c      implicit none
      double precision c(*),x
      integer k

      double precision div
      integer j

      div = 1.d0 / (k-1)
      c(k) = div*x*c(k-1)
      do 100 j = 1,k-2
       c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
100   continue
      c(1) = div*(1-x)*c(1)
      return
      end
c-------------------------------------
      subroutine diff(c,d,order)
c      implicit none
      double precision c(*),d(*)
      integer order

      integer j
      d(1) = -c(1)
      do 10 j = 2,order
       d(j) = c(j-1) - c(j)
10    continue
      return
      end
c-------------------------------------
c---------------------------------------------------------
c   beginning of public fft
c---------------------------------------------------------
*****************************************************************************
*
*       3D (slow) Fourier Transform
*   this 1d->3d code is brute force approach
*   the 1d code is a double precision version of fftpack from netlib
*   due to Paul N Swartztrauber at NCAR Boulder Coloraso
*
*****************************************************************************

      subroutine pubz3di(n1,n2,n3,table,ntable)
c      implicit none
      integer n1,n2,n3,ntable
      double precision table(ntable,3)
c ntable should be 4*max(n1,n2,n3) +15


      call cffti(n1,table(1,1))
      call cffti(n2,table(1,2))
      call cffti(n3,table(1,3))

      return
      end
*****************************************************************************
      subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable,
     $    work,nwork)
c      implicit none

      integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
      double complex w(ld1,ld2,n3)
      double complex work( nwork)
      double precision table(ntable,3)

      integer i,j,k
c ntable should be 4*max(n1,n2,n3) +15
c nwork should be max(n1,n2,n3)
c
c   transform along X  first ...
c
      do 100 k = 1, n3
       do 90 j = 1, n2
        do 70 i = 1,n1
          work(i) = w(i,j,k)
70      continue
        if ( isign .eq. -1) call cfftf(n1,work,table(1,1))
        if ( isign .eq. 1) call cfftb(n1,work,table(1,1))
        do 80 i = 1,n1
          w(i,j,k) = work(i)
80      continue
90     continue
100   continue
c
c   transform along Y then ...
c
      do 200 k = 1,n3
       do 190 i = 1,n1
        do 170 j = 1,n2
          work(j) = w(i,j,k)
170     continue
        if ( isign .eq. -1) call cfftf(n2,work,table(1,2))
        if ( isign .eq. 1) call cfftb(n2,work,table(1,2))
        do 180 j = 1,n2
          w(i,j,k) = work(j)
180     continue
190    continue
200   continue
c
c   transform along Z finally ...
c
      do 300 i = 1, n1
       do 290 j = 1, n2
        do 270 k = 1,n3
          work(k) = w(i,j,k)
270     continue
        if ( isign .eq. -1) call cfftf(n3,work,table(1,3))
        if ( isign .eq. 1) call cfftb(n3,work,table(1,3))
        do 280 k = 1,n3
          w(i,j,k) = work(k)
280     continue
290    continue
300   continue

      return
      end
c----------------------------------------------------
      SUBROUTINE CFFTB(N,C,WSAVE)
c      implicit double precision (a-h,o-z)
      double precision C,WSAVE  
      integer N,IW1,IW2 
      DIMENSION       C(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1(N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE CFFTB1(N,C,CH,WA,IFAC)
c      implicit double precision (a-h,o-z)
c      implicit none
      double precision C,CH,WA  
      integer N,ifac,nf,na,l1,iw,k1,ip,l2,ido,idot,idl1,ix2,ix3 
      integer i,ix4,nac,n2
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE CFFTF(N,C,WSAVE)
c      implicit double precision (a-h,o-z)
c      implicit none
      double precision C,WSAVE  
      integer N,IW1,IW2 
      DIMENSION       C(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE CFFTF1(N,C,CH,WA,IFAC)
c      implicit double precision (a-h,o-z)
c      implicit none
      double precision C,CH,WA  
      integer N,ifac,nf,na,l1,iw,k1,ip,l2,ido,idot,idl1,ix2,ix3 
      integer i,ix4,nac,n2
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE CFFTI(N,WSAVE)
c      implicit double precision (a-h,o-z)
c      implicit none
      double precision WSAVE    
      integer N,IW1,IW2 
      DIMENSION       WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE CFFTI1(N,WA,IFAC)
c      implicit double precision (a-h,o-z)
c      implicit none
      double precision WA,TPI,ARGH,ARG,ARGLD,FI
c      implicit integer (i-n)
      integer N,ifac,nf,l1,k1,ip,l2,ido,idot
      integer i,ntryh,nl,j,ntry,nq,nr,ib
      integer i1,ld,ipm,ii
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959d0
      ARGH = TPI/DFLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.d0
            WA(I) = 0.d0
            LD = LD+L1
            FI = 0.d0
            ARGLD = DFLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.d0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
      SUBROUTINE PASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
c      implicit none
      integer idot,ido,ip,l1,idl1,nac,nt,ipp2,ipph,idp,jc,j,i,k,idl,inc
      integer lc,idlj,ik,l,idij,idj
      double precision CC,C1,C2,CH,CH2,WA,WAR,WAI
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     $                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     $                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
      SUBROUTINE PASSB2(IDO,L1,CC,CH,WA1)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     $                WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSB3(IDO,L1,CC,CH,WA1,WA2)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,TAUR,TAUI,CI2,CI3,CR2,CR3
      double precision DR2,DR3,DI2,DI3
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     $                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5d0,.866025403784439d0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSB4(IDO,L1,CC,CH,WA1,WA2,WA3)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,WA3,TI1,TI3,TI4,TR1,TR3,TR4
      double precision CR2,CR3,CR4,CI2,CI3,CI4
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     $                WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,WA3,TI3,TI4,TR3,TR4
      double precision CR2,CR3,CR4,CI2,CI3,CI4,WA4,TI5,TR11,TI11,TI12
      double precision CI5,CR5,TR5,DI5,DR5,TR12,DR4,DI4,DI3,DR3,DR2,DI2
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     $                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947d0,
     $ .951056516295154d0,
     $-.809016994374947d0,.587785252292473d0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
c      implicit none
      integer idot,ido,ip,l1,idl1,nac,nt,ipp2,ipph,idp,jc,j,i,k,idl,inc
      integer lc,idlj,ik,l,idij,idj
      double precision CC,C1,C2,CH,CH2,WA,WAR,WAI
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     $                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     $                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
      SUBROUTINE PASSF2(IDO,L1,CC,CH,WA1)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     $                WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSF3(IDO,L1,CC,CH,WA1,WA2)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,TAUR,TAUI,CI2,CI3,CR2,CR3
      double precision DR2,DR3,DI2,DI3
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     $                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5d0,-.866025403784439d0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSF4(IDO,L1,CC,CH,WA1,WA2,WA3)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,WA3,TI1,TI3,TI4,TR1,TR3,TR4
      double precision CR2,CR3,CR4,CI2,CI3,CI4
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     $                WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
      SUBROUTINE PASSF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
      integer k,l1,i,ido
      double precision CC,CH,WA1,TR2,TI2,WA2,WA3,TI3,TI4,TR3,TR4
      double precision CR2,CR3,CR4,CI2,CI3,CI4,WA4,TI5,TR11,TI11,TI12
      double precision CI5,CR5,TR5,DI5,DR5,TR12,DR4,DI4,DI3,DR3,DR2,DI2
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     $                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947d0,
     $ -.951056516295154d0,
     $-.809016994374947d0,-.587785252292473d0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
c---------------------------------------------------------
c   end of public fft
c---------------------------------------------------------
