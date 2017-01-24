      subroutine matvec_diag(vector,vectout,npts)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the Bond contribution.
c
c vector     -  input vector
c vectout    -  output_vector (product of vector by the diagonal)
c diag       - diagonal matrix
c
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'

      double precision vector(*),vectout(*),mass
      integer i
      integer ii,jj,npts


      do 2 i=1,npts
c
c........Contribution of diagonal terms
c........Note that the diag matrix is organized as:
c........i=xx i+1=xy=yx i+2=xz=zx i+3=yy i+4=yz=zy i+5=zz
c
         ii = 3*(i-1)+1
         jj = 6*(i-1)+1
         if (MASSWEIGHT.eq.1) then
            mass = 1.0/ptms(i)
c            mass = massfac(i)*massfac(i)
         else
            mass = 1.0
         end if

         vectout(ii)   = vectout(ii)   + mass * (
     $        diag(jj)*vector(ii) 
     $        + diag(jj+1)*vector(ii+1)
     $        + diag(jj+2)*vector(ii+2) )

         vectout(ii+1) = vectout(ii+1) + mass * (
     $        diag(jj+1)*vector(ii) 
     $        + diag(jj+3)*vector(ii+1)
     $        + diag(jj+4)*vector(ii+2) )

         vectout(ii+2) = vectout(ii+2) + mass * (
     $        diag(jj+2)*vector(ii) 
     $        + diag(jj+4)*vector(ii+1)
     $        + diag(jj+5)*vector(ii+2) )
c
2     continue

      return
      end
