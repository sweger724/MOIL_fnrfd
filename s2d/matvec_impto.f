      subroutine matvec_impto(vector,vectout)
c       
c       Calculate the product d^2V/dxidxj *vector(i)
c       where d^2V is the second derivative matrix and vector(i)
c is an arbitrary vector. The length of vector MUST be 3*npt
c no test are made on that.
c This subroutine deals only with the part of the 
c off diagonal matrix that deals with the torsion-angle contribution.
c
c vector     -  input vector
c vectout    -  output_vector (product of vector by offdiagonal)
c offdiag    - offdiagonal matrix
c ifirst     - first torsional angle
c ilast      - last torsional angle
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
      integer iphi,iiphi
      integer kk,ll,mm,nn


c contribution of improper torsions
      
      do 11 iphi=ifirst,ilast
         
        kk = 3*(iimp1(iphi)-1)+1
        ll = 3*(iimp2(iphi)-1)+1
        mm = 3*(iimp3(iphi)-1)+1
        nn = 3*(iimp4(iphi)-1)+1
        iiphi = 54*(iphi-ifirst)+1

c k,l pair


        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp1(iphi)) * ptms(iimp2(iphi)) ) )
           mass = massfac(iimp1(iphi)) * massfac(iimp2(iphi))
        else
           mass = 1.0
        end if


        
        vectout(kk)   = vectout(kk)   + mass * (
     $       offdiag(iiphi)*vector(ll) + 
     $       offdiag(iiphi+3)*vector(ll+1) + 
     $       offdiag(iiphi+6)*vector(ll+2) )

        vectout(kk+1) = vectout(kk+1) + mass * (
     $       offdiag(iiphi+1)*vector(ll) +
     $       offdiag(iiphi+4)*vector(ll+1) + 
     $       offdiag(iiphi+7)*vector(ll+2) )

        vectout(kk+2) = vectout(kk+2) + mass * (
     $       offdiag(iiphi+2)*vector(ll) +
     $       offdiag(iiphi+5)*vector(ll+1) + 
     $       offdiag(iiphi+8)*vector(ll+2) )

        vectout(ll)   = vectout(ll)   + mass * (
     $       offdiag(iiphi)*vector(kk) +
     $       offdiag(iiphi+1)*vector(kk+1) +
     $       offdiag(iiphi+2)*vector(kk+2) )

        vectout(ll+1) = vectout(ll+1) + mass * (
     $       offdiag(iiphi+3)*vector(kk) +
     $       offdiag(iiphi+4)*vector(kk+1) + 
     $       offdiag(iiphi+5)*vector(kk+2) )

        vectout(ll+2) = vectout(ll+2) + mass * (
     $       offdiag(iiphi+6)*vector(kk) +
     $       offdiag(iiphi+7)*vector(kk+1) +
     $       offdiag(iiphi+8)*vector(kk+2) )


c k,m

        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp1(iphi)) * ptms(iimp3(iphi)) ) )
           mass = massfac(iimp1(iphi)) * massfac(iimp3(iphi))
        else
           mass = 1.0
        end if

        vectout(kk)   = vectout(kk)   + mass * (
     $       offdiag(iiphi+9)*vector(mm) +
     $       offdiag(iiphi+12)*vector(mm+1) +
     $       offdiag(iiphi+15)*vector(mm+2) )

        vectout(kk+1) = vectout(kk+1) + mass * (
     $       offdiag(iiphi+10)*vector(mm) +
     $       offdiag(iiphi+13)*vector(mm+1) +
     $       offdiag(iiphi+16)*vector(mm+2) )

        vectout(kk+2) = vectout(kk+2) + mass * (
     $       offdiag(iiphi+11)*vector(mm) +
     $       offdiag(iiphi+14)*vector(mm+1) + 
     $       offdiag(iiphi+17)*vector(mm+2) )

        vectout(mm)   = vectout(mm)   + mass * (
     $       offdiag(iiphi+9)*vector(kk) +
     $       offdiag(iiphi+10)*vector(kk+1) + 
     $       offdiag(iiphi+11)*vector(kk+2) )

        vectout(mm+1) = vectout(mm+1) + mass * (
     $       offdiag(iiphi+12)*vector(kk) +
     $       offdiag(iiphi+13)*vector(kk+1) + 
     $       offdiag(iiphi+14)*vector(kk+2) )

        vectout(mm+2) = vectout(mm+2) + mass * (
     $       offdiag(iiphi+15)*vector(kk) +
     $       offdiag(iiphi+16)*vector(kk+1) +
     $       offdiag(iiphi+17)*vector(kk+2) )

c k,n


        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp1(iphi)) * ptms(iimp4(iphi)) )) 
           mass = massfac(iimp1(iphi)) * massfac(iimp4(iphi))
        else
           mass = 1.0
        end if

        vectout(kk)   = vectout(kk)   + mass * (
     $       offdiag(iiphi+18)*vector(nn) +
     $       offdiag(iiphi+21)*vector(nn+1) + 
     $       offdiag(iiphi+24)*vector(nn+2) )

        vectout(kk+1) = vectout(kk+1) + mass * (
     $       offdiag(iiphi+19)*vector(nn) +
     $       offdiag(iiphi+22)*vector(nn+1) +
     $       offdiag(iiphi+25)*vector(nn+2) )

        vectout(kk+2) = vectout(kk+2) + mass * (
     $       offdiag(iiphi+20)*vector(nn) +
     $       offdiag(iiphi+23)*vector(nn+1) +
     $       offdiag(iiphi+26)*vector(nn+2) )

        vectout(nn)   = vectout(nn)   + mass * (
     $       offdiag(iiphi+18)*vector(kk) +
     $       offdiag(iiphi+19)*vector(kk+1) +
     $       offdiag(iiphi+20)*vector(kk+2) )

        vectout(nn+1) = vectout(nn+1) + mass * (
     $       offdiag(iiphi+21)*vector(kk) +
     $       offdiag(iiphi+22)*vector(kk+1) + 
     $       offdiag(iiphi+23)*vector(kk+2) )

        vectout(nn+2) = vectout(nn+2) + mass * (
     $       offdiag(iiphi+24)*vector(kk) +
     $       offdiag(iiphi+25)*vector(kk+1) + 
     $       offdiag(iiphi+26)*vector(kk+2) )

c l,m

        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp2(iphi)) * ptms(iimp3(iphi)) ) )
           mass = massfac(iimp2(iphi)) * massfac(iimp3(iphi))
        else
           mass = 1.0
        end if


        vectout(ll)   = vectout(ll)   + mass * (
     $       offdiag(iiphi+27)*vector(mm) +
     $       offdiag(iiphi+30)*vector(mm+1) + 
     $       offdiag(iiphi+33)*vector(mm+2) )

        vectout(ll+1) = vectout(ll+1) + mass * (
     $       offdiag(iiphi+28)*vector(mm) +
     $       offdiag(iiphi+31)*vector(mm+1) + 
     $       offdiag(iiphi+34)*vector(mm+2) )

        vectout(ll+2) = vectout(ll+2) + mass * (
     $       offdiag(iiphi+29)*vector(mm) +
     $       offdiag(iiphi+32)*vector(mm+1) +
     $       offdiag(iiphi+35)*vector(mm+2) )

        vectout(mm)   = vectout(mm)   + mass * (
     $       offdiag(iiphi+27)*vector(ll) +
     $       offdiag(iiphi+28)*vector(ll+1) +
     $       offdiag(iiphi+29)*vector(ll+2) )

        vectout(mm+1) = vectout(mm+1) + mass * (
     $       offdiag(iiphi+30)*vector(ll) +
     $       offdiag(iiphi+31)*vector(ll+1) + 
     $       offdiag(iiphi+32)*vector(ll+2) )

        vectout(mm+2) = vectout(mm+2) + mass * (
     $       offdiag(iiphi+33)*vector(ll) +
     $       offdiag(iiphi+34)*vector(ll+1) +
     $       offdiag(iiphi+35)*vector(ll+2) )

c l,n

        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp2(iphi)) * ptms(iimp4(iphi)) ) )
           mass = massfac(iimp2(iphi)) * massfac(iimp4(iphi))
        else
           mass = 1.0
        end if

        vectout(ll)   = vectout(ll)   + mass * (
     $       offdiag(iiphi+36)*vector(nn) +
     $       offdiag(iiphi+39)*vector(nn+1) + 
     $       offdiag(iiphi+42)*vector(nn+2) )

        vectout(ll+1) = vectout(ll+1) + mass * (
     $       offdiag(iiphi+37)*vector(nn) +
     $       offdiag(iiphi+40)*vector(nn+1) + 
     $       offdiag(iiphi+43)*vector(nn+2) )

        vectout(ll+2) = vectout(ll+2) + mass * (
     $       offdiag(iiphi+38)*vector(nn) +
     $       offdiag(iiphi+41)*vector(nn+1) +
     $       offdiag(iiphi+44)*vector(nn+2) )

        vectout(nn)   = vectout(nn)   + mass * (
     $       offdiag(iiphi+36)*vector(ll) +
     $       offdiag(iiphi+37)*vector(ll+1) + 
     $       offdiag(iiphi+38)*vector(ll+2) )

        vectout(nn+1) = vectout(nn+1) + mass * (
     $       offdiag(iiphi+39)*vector(ll) +
     $       offdiag(iiphi+40)*vector(ll+1) + 
     $       offdiag(iiphi+41)*vector(ll+2) )

        vectout(nn+2) = vectout(nn+2) + mass * (
     $       offdiag(iiphi+42)*vector(ll) +
     $       offdiag(iiphi+43)*vector(ll+1) +
     $       offdiag(iiphi+44)*vector(ll+2) )

c m,n

        if (MASSWEIGHT.eq.1) then
c           mass = 1.0/( sqrt(ptms(iimp3(iphi)) * ptms(iimp4(iphi)) ) )
           mass = massfac(iimp3(iphi))* massfac(iimp4(iphi))
        else
           mass = 1.0
        end if

        vectout(mm)   = vectout(mm)   + mass * (
     $       offdiag(iiphi+45)*vector(nn) +
     $       offdiag(iiphi+48)*vector(nn+1) + 
     $       offdiag(iiphi+51)*vector(nn+2) )

        vectout(mm+1) = vectout(mm+1) + mass * (
     $       offdiag(iiphi+46)*vector(nn) +
     $       offdiag(iiphi+49)*vector(nn+1) + 
     $       offdiag(iiphi+52)*vector(nn+2) )

        vectout(mm+2) = vectout(mm+2) + mass * (
     $       offdiag(iiphi+47)*vector(nn) +
     $       offdiag(iiphi+50)*vector(nn+1) + 
     $       offdiag(iiphi+53)*vector(nn+2) )

        vectout(nn)   = vectout(nn)   + mass * (
     $       offdiag(iiphi+45)*vector(mm) +
     $       offdiag(iiphi+46)*vector(mm+1) + 
     $       offdiag(iiphi+47)*vector(mm+2) )

        vectout(nn+1) = vectout(nn+1) + mass * (
     $       offdiag(iiphi+48)*vector(mm) +
     $       offdiag(iiphi+49)*vector(mm+1) + 
     $       offdiag(iiphi+50)*vector(mm+2) )

        vectout(nn+2) = vectout(nn+2) + mass * (
     $       offdiag(iiphi+51)*vector(mm) +
     $       offdiag(iiphi+52)*vector(mm+1) +
     $       offdiag(iiphi+53)*vector(mm+2) )

11    continue

      return
      end
