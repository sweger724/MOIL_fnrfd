      subroutine matvec_theta(vector,vectout)
c 
c Calculate the product d^2V/dxidxj *vector(i)
c where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the Bond-angle contribution.
c
c vector     -  input vector
c vectout    -  output_vector (product of vector by offdiagonal)
c offdiag    - offdiagonal matrix
c ifirst     - first bond-angle
c ilast      - last bond-angle
c
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/PDQ.BLOCK'
	include 'COMMON/MASSFAC.BLOCK'
      double precision mass
      double precision vector(*),vectout(*)
c       

      integer ith,iith
      integer ii,jj,kk
c       
c Angle contribution
c
      do 9 ith=ifirst,ilast
         iith= 27*(ith-ifirst)+1
         ii  = 3*(iangl1(ith)-1)+1
         jj  = 3*(iangl2(ith)-1)+1
         kk  = 3*(iangl3(ith)-1)+1
c pair i,j

	 if (MASSWEIGHT.eq.1) then
c	    mass = 1.0/sqrt(ptms(iangl1(ith))*ptms(iangl2(ith)) )
            mass = massfac(iangl1(ith)) * massfac (iangl2(ith))
	 else
	    mass = 1.0
	 end if


         vectout(ii) = vectout(ii) + mass * (
     $        offdiag(iith)*vector(jj) +
     $        offdiag(iith+1)*vector(jj+1) + 
     $        offdiag(iith+2)*vector(jj+2) )


         vectout(ii+1) = vectout(ii+1) + mass * (
     $        offdiag(iith+3)*vector(jj) +
     $        offdiag(iith+4)*vector(jj+1) +
     $        offdiag(iith+5)*vector(jj+2) )

         vectout(ii+2) = vectout(ii+2) + mass * (
     $        offdiag(iith+6)*vector(jj) +
     $        offdiag(iith+7)*vector(jj+1) +
     $        offdiag(iith+8)*vector(jj+2) )

         vectout(jj) = vectout(jj) + mass * (
     $        offdiag(iith)*vector(ii) +
     $        offdiag(iith+3)*vector(ii+1) +
     $        offdiag(iith+6)*vector(ii+2) )
         vectout(jj+1) = vectout(jj+1) + mass * (
     $        offdiag(iith+1)*vector(ii) +
     $        offdiag(iith+4)*vector(ii+1) +
     $        offdiag(iith+7)*vector(ii+2) )

         vectout(jj+2) = vectout(jj+2) + mass * (
     $        offdiag(iith+2)*vector(ii) +
     $        offdiag(iith+5)*vector(ii+1) +
     $        offdiag(iith+8)*vector(ii+2) )

c pair i,k

	 if (MASSWEIGHT.eq.1) then
c	    mass = 1.0/sqrt(ptms(iangl1(ith))*ptms(iangl3(ith)) )
            mass = massfac(iangl1(ith)) * massfac(iangl3(ith))
	 else
	    mass = 1.0
	 end if


         vectout(ii) = vectout(ii) + mass * (
     $        offdiag(iith+9)*vector(kk) +
     $        offdiag(iith+10)*vector(kk+1) + 
     $        offdiag(iith+11)*vector(kk+2) )
         
         vectout(ii+1) = vectout(ii+1) + mass * (
     $        offdiag(iith+12)*vector(kk) +
     $        offdiag(iith+13)*vector(kk+1) + 
     $        offdiag(iith+14)*vector(kk+2) )

         vectout(ii+2) = vectout(ii+2) + mass * (
     $        offdiag(iith+15)*vector(kk) +
     $        offdiag(iith+16)*vector(kk+1) + 
     $        offdiag(iith+17)*vector(kk+2) )

         vectout(kk) = vectout(kk) + mass * ( 
     $        offdiag(iith+9)*vector(ii) +
     $     offdiag(iith+12)*vector(ii+1) + 
     $        offdiag(iith+15)*vector(ii+2) )

         vectout(kk+1) = vectout(kk+1) + mass * (
     $        offdiag(iith+10)*vector(ii) +
     $        offdiag(iith+13)*vector(ii+1) + 
     $        offdiag(iith+16)*vector(ii+2))

         vectout(kk+2) = vectout(kk+2) + mass * (
     $        offdiag(iith+11)*vector(ii) +
     $        offdiag(iith+14)*vector(ii+1) + 
     $        offdiag(iith+17)*vector(ii+2) )

c pair j,k

	 if (MASSWEIGHT.eq.1) then
c	    mass = 1.0/sqrt(ptms(iangl2(ith))*ptms(iangl3(ith)) )
            mass = massfac(iangl2(ith)) * massfac(iangl3(ith))
	 else
	    mass = 1.0
	 end if


         vectout(jj) = vectout(jj) + mass * (
     $        offdiag(iith+18)*vector(kk) +
     $        offdiag(iith+19)*vector(kk+1) + 
     $        offdiag(iith+20)*vector(kk+2) )

         vectout(jj+1) = vectout(jj+1) + mass * (
     $        offdiag(iith+21)*vector(kk) +
     $        offdiag(iith+22)*vector(kk+1) + 
     $        offdiag(iith+23)*vector(kk+2) )

         vectout(jj+2) = vectout(jj+2) + mass * (
     $        offdiag(iith+24)*vector(kk) +
     $        offdiag(iith+25)*vector(kk+1) + 
     $        offdiag(iith+26)*vector(kk+2) )
         
         vectout(kk) = vectout(kk) + mass * (
     $        offdiag(iith+18)*vector(jj) +
     $        offdiag(iith+21)*vector(jj+1) + 
     $        offdiag(iith+24)*vector(jj+2) )

         vectout(kk+1) = vectout(kk+1) + mass * (
     $        offdiag(iith+19)*vector(jj) +
     $        offdiag(iith+22)*vector(jj+1) + 
     $        offdiag(iith+25)*vector(jj+2) )

         vectout(kk+2) = vectout(kk+2) + mass * (
     $        offdiag(iith+20)*vector(jj) +
     $        offdiag(iith+23)*vector(jj+1) + 
     $        offdiag(iith+26)*vector(jj+2) )

9     continue
c
      return
      end











