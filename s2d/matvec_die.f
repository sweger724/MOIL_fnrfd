      subroutine matvec_die(vector,vectout,sepfast,npts)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the non-bonded contribution.
c (very similar to the bond contribution subroutine)
c
c vector     -  input vector
c vectout    -  output_vector (product of vector by offdiagonal)
c offdiag    - offdiagonal matrix
c ifirst     - first pair
c ilast      - last pair
c
c
      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'
      double precision vector(*),vectout(*),mass
c
      integer i,j,k
      integer ii,jj,kk,npts
      logical sepfast	

C       Skip the mass weighting on the FAST Modes

      if (sepfast) then
c
c if sepfast=.true. offdiag contribution involving fast modes (waters)
c should not be calculated
c

      do 10 k=ifirst,ilast

         i=pairs1(k)
         if (i.gt.npts) go to 10
         j=pairs2(k)
         if (j.gt.npts) go to 10
         kk = 6*(k-ifirst)+1
         ii = 3*(i-1)+1
         jj = 3*(j-1)+1
c
c Off diagonal elements: for each pair i and j, the storage is
c kk=xx kk+1=xy=yx kk+2=xz=zx kk+3=yy kk+4=yz=zy kk+5=zz
c
         vectout(ii) = vectout(ii) + offdiag(kk)*vector(jj) +
     1                  offdiag(kk+1)*vector(jj+1) +
     2                  offdiag(kk+2)*vector(jj+2)
         vectout(ii+1) = vectout(ii+1) + offdiag(kk+1)*vector(jj) +
     1                  offdiag(kk+3)*vector(jj+1) +
     2                  offdiag(kk+4)*vector(jj+2)
         vectout(ii+2) = vectout(ii+2) + offdiag(kk+2)*vector(jj) +
     1                  offdiag(kk+4)*vector(jj+1) +
     2                  offdiag(kk+5)*vector(jj+2)
         vectout(jj) = vectout(jj) + offdiag(kk)*vector(ii) +
     1                  offdiag(kk+1)*vector(ii+1) +
     2                  offdiag(kk+2)*vector(ii+2)
         vectout(jj+1) = vectout(jj+1) + offdiag(kk+1)*vector(ii) +
     1                  offdiag(kk+3)*vector(ii+1) +
     2                  offdiag(kk+4)*vector(ii+2)
         vectout(jj+2) = vectout(jj+2) + offdiag(kk+2)*vector(ii) +
     1                  offdiag(kk+4)*vector(ii+1) +
     2                  offdiag(kk+5)*vector(ii+2)
c
 10   continue
c

      else
c i.e. when sepfast=.false.

      do 20 k=ifirst,ilast

         i=pairs1(k)
         j=pairs2(k)
         kk = 6*(k-ifirst)+1
         ii = 3*(i-1)+1
         jj = 3*(j-1)+1 
c
c Off diagonal elements: for each pair i and j, the storage is
c kk=xx kk+1=xy=yx kk+2=xz=zx kk+3=yy kk+4=yz=zy kk+5=zz
c

         if (MASSWEIGHT.eq.1) then
c            mass = 1.0/sqrt(ptms(i) * ptms(j))
            mass = massfac(i)*massfac(j)
         else 
            mass = 1.0
         end if

         vectout(ii) = vectout(ii) + mass * (
     $        offdiag(kk)*vector(jj) +
     $        offdiag(kk+1)*vector(jj+1) +  
     $        offdiag(kk+2)*vector(jj+2) )

         vectout(ii+1) = vectout(ii+1) + mass * (
     $        offdiag(kk+1)*vector(jj) +
     $        offdiag(kk+3)*vector(jj+1) +  
     $        offdiag(kk+4)*vector(jj+2) )

         vectout(ii+2) = vectout(ii+2) + mass * (
     $        offdiag(kk+2)*vector(jj) +
     $        offdiag(kk+4)*vector(jj+1) +  
     $        offdiag(kk+5)*vector(jj+2) )

         vectout(jj) = vectout(jj) + mass * (
     $        offdiag(kk)*vector(ii) +
     $        offdiag(kk+1)*vector(ii+1) +  
     $        offdiag(kk+2)*vector(ii+2) )

         vectout(jj+1) = vectout(jj+1) + mass * (
     $        offdiag(kk+1)*vector(ii) +
     $        offdiag(kk+3)*vector(ii+1) + 
     $        offdiag(kk+4)*vector(ii+2) )

         vectout(jj+2) = vectout(jj+2) + mass * (
     $        offdiag(kk+2)*vector(ii) +
     $        offdiag(kk+4)*vector(ii+1) +  
     $        offdiag(kk+5)*vector(ii+2) )

c
 20   continue
c
      end if
c
      return
      end


c---------------------------------------------------------------------
      subroutine matvec_wat(vector,vectout)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the non-bonded contribution.
c (very similar to the bond contribution subroutine)
c
c vector     -  input vector
c vectout    -  output_vector (product of vector by offdiagonal)
c offdiag    - offdiagonal matrix
c ifirst     - first pair
c ilast      - last pair
c
c     
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'

      double precision vector(*),vectout(*),mass
c
      integer i,j,k
      integer ii,jj,kk
      

      do 20 k=ifirst,ilast


         i=pairs1(k)
         j=pairs2(k)

        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/sqrt(ptms(i) * ptms(j))
           mass = massfac(i)*massfac(j)
        else 
           mass = 1.0
        end if

         kk = 6*(k-ifirst)+1
         ii = 3*(i-1)+1
         jj = 3*(j-1)+1 
c
c Off diagonal elements: for each pair i and j, the storage is
c kk=xx kk+1=xy=yx kk+2=xz=zx kk+3=yy kk+4=yz=zy kk+5=zz
c
         vectout(ii) = vectout(ii) + mass * (
     $        offdiag(kk)*vector(jj) +
     $        offdiag(kk+1)*vector(jj+1) +  
     $        offdiag(kk+2)*vector(jj+2) )

         vectout(ii+1) = vectout(ii+1) + mass * (
     $        offdiag(kk+1)*vector(jj) +
     $        offdiag(kk+3)*vector(jj+1) +  
     $        offdiag(kk+4)*vector(jj+2) )

         vectout(ii+2) = vectout(ii+2) + mass * (
     $        offdiag(kk+2)*vector(jj) +
     $        offdiag(kk+4)*vector(jj+1) +  
     $        offdiag(kk+5)*vector(jj+2) )

         vectout(jj) = vectout(jj) + mass * (
     $        offdiag(kk)*vector(ii) +
     $        offdiag(kk+1)*vector(ii+1) +  
     $        offdiag(kk+2)*vector(ii+2) )

         vectout(jj+1) = vectout(jj+1) + mass * (
     $        offdiag(kk+1)*vector(ii) +
     $        offdiag(kk+3)*vector(ii+1) + 
     $        offdiag(kk+4)*vector(ii+2) )

         vectout(jj+2) = vectout(jj+2) + mass * (
     $        offdiag(kk+2)*vector(ii) +
     $        offdiag(kk+4)*vector(ii+1) +  
     $        offdiag(kk+5)*vector(ii+2) )
c
 20   continue
c
c
      return
      end





