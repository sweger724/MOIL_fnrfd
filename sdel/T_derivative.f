      subroutine T_deriv(dAA,r,d0,p,dv,dS,pseg,npt)

      implicit none


      integer npt, pseg
      double precision dAA(3,*),r(3,*),d0(*),p(*),dv(3,*),dS(3,*)

c     
c     common block for COORDinates and potential ENERGY derivatives
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/SDEL.BLOCK'

C     Use this for calling matvec
      include 'COMMON/SCNDRV_SDEL.BLOCK'
      include 'COMMON/GETVEC.BLOCK'


c     pdq variables

C     p() carries our momenta of each sys.
C     term0,term1,term2 are temp. variables for calculating dS.dS
      double precision minuscontrib,pluscontrib,diagcontrib

c     matvec routine variable
      double precision mtvc2(3,MAXPT)

      integer j,jp1,jm1

      double precision ukjm1,ukj,ukjp1
      double precision qkjm1,qkj

      double precision umj,qmjm1,qmj

      integer i,k,l,m
      double precision minusparts(4),plusparts(4),diagparts(7)


C      write(6,*)"Entered T_derivative."

         do i = 2,pseg+1
C     get the start location of each grid point
C     for the derivatives dA/dq we need to go 2 above and two below

                  jp1 = (i  ) * npt
                  j   = (i-1) * npt
                  jm1 = (i-2) * npt
                  
C     Now calculate dS.dS

C     These are the offdiagonal and diagonal terms we calculate in our sum

                  do l = 1,4
                     minusparts(l) = 0.0d0
                     plusparts(l) = 0.0d0
                  end do

                  do l = 1,7
                     diagparts(l) = 0.0d0
                  end do

                  
                  do k = 1,npt
C     here are parts of our dS.dS sum
                     do l = 1,3

                        ukjm1 = dv(l,jm1+k) / p(i-1)
                        ukj   = dv(l,j+k  ) / p(i)
                        ukjp1 = dv(l,jp1+k) / p(i+1)

                        qkjm1 = (r(l,jm1+k) - r(l,j+k))   / d0(i-1)
                        qkj   = (r(l,j+k)   - r(l,jp1+k)) / d0(i)

C     here are parts of the contribution from the off diagonal (-1)
      minusparts(1) =minusparts(1) - 2*dS(l,jm1+k) * qkjm1
      minusparts(2) =minusparts(2) + 2*dS(l,jm1+k) * ukjm1
      minusparts(3) =minusparts(3) + 2*dS(l,jm1+k)*p(i)/d0(i-1)*qkjm1
      minusparts(4) =minusparts(4) + 2*dS(l,jm1+k)*p(i-1)/d0(i-1)*qkjm1
                           

C     Now do the offdiagonal +1 side
      plusparts(1) = plusparts(1) + 2*dS(l,jp1+k) * qkj
      plusparts(2) = plusparts(2) - 2*dS(l,jp1+k) * ukjp1
      plusparts(3) = plusparts(3) + 2*dS(l,jp1+k) * p(i)/d0(i) * qkj
      plusparts(4) = plusparts(4) + 2*dS(l,jp1+k) * p(i+1)/d0(i)*qkj

C     diagparts sum
      diagparts(1) = diagparts(1) + 2*dS(l,j+k) * (qkj - qkjm1)
      diagparts(2) = diagparts(2) + 2*dS(l,j+k) * ukj 
      diagparts(3) = diagparts(3) + 2*dS(l,j+k) * p(i)/d0(i) * qkj
      diagparts(4) = diagparts(4) + 2*dS(l,j+k) * p(i+1)/d0(i)* qkj
      diagparts(5) = diagparts(5) + 2*dS(l,j+k) * p(i)/d0(i-1)*qkjm1
      diagparts(6) = diagparts(6) + 2*dS(l,j+k) *p(i-1)/d0(i-1)*qkjm1
      diagparts(7)=  diagparts(7) + 2*dS(l,j+k) * 
     $                              ukj/p(i)*(d0(i) + d0(i-1))

C     end loop    l = 1, 3
                     end do 
C     end loop    k = 1, npt
                  end do


C     Do the multiplication
C     This is what we would do if we called getvec() (hopefully ;))

C     Now we want to call getvec so copy in our non-massweighted coordinates 
               do m = 1,npt
                 do l = 1,3
C                   coor(l,m) = r(l,j+m)*massfac(m)
                    mtvc2(l,m) = -2*(d0(i) + d0(i-1))/p(i)*dS(l,j+m)
                 end do
               end do

C               call getvec(mtvc2,.false.)
                call Getvec_Approx(mtvc2(1,1),dv(1,1+j),r(1,j+1))


C     Do the multiplication
                  
C     calculate the derivatives (this is the m loop in the notes)

        do m = 1 ,npt
          do l = 1,3
                     umj = dv(l,j+m)/p(i)
                   qmjm1 =  (r(l,jm1+m) - r(l,j+m)  ) / d0(i-1)
                   qmj   =  (r(l,j+m)   - r(l,jp1+m)) / d0(i)         
                        
c     Now lets actually add up our derivatives

C     -1 off diagonal term
      minuscontrib = minusparts(1) * umj + minusparts(2) * qmjm1 +
     $               minusparts(3) * qmjm1 + minusparts(4) * qmjm1
      
      minuscontrib = minuscontrib - 2*dS(l,jm1+m) * p(i)/d0(i-1)
     $             - 2*dS(l,jm1+m)*p(i-1)/d0(i-1)

C     +1 off diagonal term
       pluscontrib = plusparts(1) * umj + plusparts(2) * qmj +
     $               plusparts(3) * qmj + plusparts(4) * qmj
       
       pluscontrib = pluscontrib - 2*dS(l,jp1+m) * p(i)/d0(i)
     $             - 2*dS(l,jp1+m) * p(i+1)/d0(i)
                     
                     
       diagcontrib = - diagparts(1) * umj - diagparts(2) *(qmj-qmjm1)
     $               - diagparts(3) * qmj - diagparts(4) * qmj 
     $               - diagparts(5) * qmjm1 - diagparts(6) * qmjm1
     $               - diagparts(7) * umj

C     Now do Kronecker delta fn terms
      diagcontrib = diagcontrib + 2*dS(l,j+m) * p(i)/d0(i) 
     $            + 2*dS(l,j+m) * p(i+1)/d0(i)
     $            + 2*dS(l,j+m) * p(i)/d0(i-1)
     $            + 2*dS(l,j+m) * p(i-1)/d0(i-1)


                      diagcontrib = diagcontrib + mtvc2(l,m)

C     Add result 
                  dAA(l,j+m) = dAA(l,j+m) + minuscontrib 
     $                       + pluscontrib + diagcontrib

         end do
       end do
                  
C     end loop i = 2,pseg+1
      end do
               
C    make sure that the remainder is 0         
               do j = 1,npt
                 do l = 1,3
                  dAA(l,j) = 0.0d0
                  dAA(l,(pseg+1)*npt+j) = 0.0d0
                 end do 
               end do

789        return
            end
