      subroutine vecmin_upw(x, y, z,n, xmy)

c     vector subtraction by upwind 
c     IN
c     x, y: equal-length vectors
c     n: length of x & y
c
c     OUT
c     xmy: x - y

      integer i, n
      double precision x(3,*), y(3,*), z(3,*),xmy(3,*)
      double precision fact , e1, e2,tmp1,tmp2
      real sigma 

c    x next mst 
c    z current mst 
c    y prev mst  
c    xmy tangent vector 


      do 100 i = 1, n
         do j=1,3


        if(x(j,i).ge.z(j,i).and.z(j,i).ge.y(j,i)) then
            xmy(j,i)=x(j,i)-z(j,i)

        else if(x(j,i).le.z(j,i).and.z(j,i).le.y(j,i)) then
            xmy(j,i)=z(j,i)-y(j,i)
  
        else
          tmp1=dabs(y(j,i)-z(j,i))
          tmp2=dabs(z(j,i)-x(j,i))
          e1=dmax1(tmp1,tmp2)
          e2=dmin1(tmp1,tmp2)
          if(x(j,i).ge.y(j,i)) then
            
              xmy(j,i)=e1*(x(j,i)-z(j,i))+
     .             e2*(z(j,i)-y(j,i))
  4         continue
          else
          
              xmy(j,i)=e2*(x(j,i)-z(j,i))
     .             +e1*(z(j,i)-y(j,i))
  
          endif
        endif


        enddo 
c         fact=(x(j,i)-z(j,i))*z(j,i)
c          if ( fact.ge.0) then 
c             xmy(j,i)= x(j,i) - z(j,i)
c          else
c             xmy(j,i)= y(j,i) - z(j,i) 
c          endif 
         
 100  continue

      end
