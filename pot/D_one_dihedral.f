      double precision function D_one_dihedral(i,j,k,l,coor,D_dihedral)
c
c compute a dihedral angle and its deirivative of for sequential atoms
c following the aper by Blondel and Karplus JCC 17,1132-1141(1996)
c
      implicit none
      integer i,j,k,l
      double precision coor(3,*), D_dihedral(3,4)

c local

      integer m,n
      double precision F(3),G(3),H(3),A(3),B(3)
      double precision Anorm, Bnorm, Gnorm
      double precision cosphi,sinphi
      double precision t1,t2,t3,t4
      double precision vscalar, vnorm2

      do m=1,3
       F(m) = coor(m,i) - coor(m,j)
       G(m) = coor(m,j) - coor(m,k)
       H(m) = coor(m,l) - coor(m,k)
      end do

      A(1) = F(2)*G(3)-F(3)*G(2)
      A(2) = -F(1)*G(3)+F(3)*G(1)
      A(3) = F(1)*G(2)-F(2)*G(1)

      
      B(1) = H(2)*G(3)-H(3)*G(2)
      B(2) = -H(1)*G(3)+H(3)*G(1)
      B(3) = H(1)*G(2)-H(2)*G(1)

      Anorm = vnorm2(A,3)
      Bnorm = vnorm2(B,3)
      Gnorm = vnorm2(G,3)
      if ((Anorm.eq.0) .or. (Bnorm .eq.0) .or. (Gnorm .eq.0)) then
       write(*,*)' zero vector norms A B G ',Anorm,Bnorm,Gnorm
       write(*,*) ' in routine D_one_dihedral '
       stop
      end if

      cosphi = 0.d0
      do m=1,3
       cosphi = cosphi + A(m)*B(m)
      end do
      cosphi = cosphi/(Anorm*Bnorm)

      sinphi = (B(2)*A(3)-B(3)*A(2))*G(1)
      sinphi = sinphi - (B(1)*A(3)-B(3)*A(1))*G(2)
      sinphi = sinphi + (B(1)*A(2)-B(2)*A(1))*G(3)
      sinphi = sinphi/(Anorm*Bnorm*Gnorm)
      D_one_dihedral = acos(cosphi)
      if (sinphi.lt.0) D_one_dihedral = - D_one_dihedral
      do m =1,3
       t1 = -Gnorm/(Anorm*Anorm)*A(m)
       t2 = vscalar(F,G,3)*A(m)/(Anorm*Anorm*Gnorm)
       t3 = -vscalar(H,G,3)*B(m)/(Bnorm*Bnorm*Gnorm)
       t4 = Gnorm*B(m)/(Bnorm*Bnorm)
       D_dihedral(m,1) = t1
       D_dihedral(m,2) = -t1+t2+t3
       D_dihedral(m,3) = -t3-t2-t4
       D_dihedral(m,4) = t4
      end do
      return
      end
