      subroutine matvec_bond(vector,vectout)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that. 
c
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the Bond contribution.
c
c vector     -  input vector
c vectout -  output_vector (product of vector by offdiagonal)
c offdiag    - offdiagonal matrix
c ifirst     - first bond
c ilast      - last bond
c
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'

      double precision vector(*),vectout(*),mass
c
      integer k
      integer ii,jj,kk
c
c
c Bond contribution

      do 10 k=ifirst,ilast


c
         kk = 6*(k-ifirst)+1
         ii = 3*(ib1(k)-1)+1
         jj = 3*(ib2(k)-1)+1

         if (MASSWEIGHT.eq.1) then
c             mass = 1.0/( sqrt(ptms(ib1(k)) * ptms(ib2(k)) ) )
            mass = massfac(ib1(k)) * massfac(ib2(k))
         else
            mass = 1.0
         end if
         vectout(ii)   = vectout(ii) + 
     $        mass * ( offdiag(kk)*vector(jj) +
     $        offdiag(kk+1)*vector(jj+1) + 
     $        offdiag(kk+2)*vector(jj+2) )

         vectout(ii+1) = vectout(ii+1)+ 
     $        mass * ( offdiag(kk+1)*vector(jj) +
     $        offdiag(kk+3)*vector(jj+1) + 
     $        offdiag(kk+4)*vector(jj+2) )

         vectout(ii+2) = vectout(ii+2) +
     $        mass * ( offdiag(kk+2)*vector(jj) +
     $        offdiag(kk+4)*vector(jj+1) +  
     $        offdiag(kk+5)*vector(jj+2) )

	 vectout(jj)   = vectout(jj) + 
     $        mass * ( offdiag(kk)*vector(ii) +
     $        offdiag(kk+1)*vector(ii+1) +  
     $        offdiag(kk+2)*vector(ii+2) )

         vectout(jj+1) = vectout(jj+1) + 
     $        mass * ( offdiag(kk+1)*vector(ii) +
     $        offdiag(kk+3)*vector(ii+1) +  
     $        offdiag(kk+4)*vector(ii+2) )

         vectout(jj+2) = vectout(jj+2) + mass * (
     $        offdiag(kk+2)*vector(ii) +
     $        offdiag(kk+4)*vector(ii+1) +  
     $        offdiag(kk+5)*vector(ii+2) )
c

         if (debugpdq) then
c            write (6,*) 'VECTOUT ',vectout(jj),vectout(jj+1),
c     $           vectout(jj+2)
c            write (6,*) 'VECTOUT ',vectout(ii),vectout(ii+1),
c     $           vectout(ii+2)
         end if
 10   continue
c
c
      return
      end










