      subroutine network(g,dg,list1,list2,rreq,g_ref,yesalpha)
      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ELASTIC.BLOCK'

c calculate bond energies and forces.


c note that the vector.block contains temporary vectors used for
c vectorization. 
c
      integer list1(0:*),list2(*)
      double precision g, dg(3,*),rreq(*),g_ref
      double precision rx,ry,rz,r2,s,r,db,df,e
      logical yesalpha
      integer i,j,jj

c e_bond=total bond energy 
c rx,ry,rz = distance squared between particles in that axis
c r2 = distance squared between particles
c r = distance between particles
c s = 1/r
c db,df temporary variables
c e = bond energy (non-acumulated)
c also used : ichunk for vectorizacion length
c xtmp,ytmp,ztemp for vectorization purposes


      g = 0.d0
      if (yesalpha) g = g + enm_alpha
      
      do i = 1, npt
         dg(1,i)=0.d0
         dg(2,i)=0.d0
         dg(3,i)=0.d0
      end do

c initialize loop over bonds in ichunk chunks
C      write(6,*)'Next structure:'
      do 100 i=1,npt-1
C            write(6,*)i,list1(i),list2(list1(i)),rreq(list1(i))
            do 200 jj=list1(i-1)+1,list1(i)
            
                  j = list2(jj)
                  rx=coor(1,i)-coor(1,j)
                  ry=coor(2,i)-coor(2,j)
                  rz=coor(3,i)-coor(3,j)
                  r2=rx*rx + ry*ry + rz*rz
                  r=dsqrt(r2)
                  db=r-rreq(jj)
C                  write(6,*)'ddb: ',db
                  df=0.5d0*enm_gamma*db
                  e=df*db
c@@@
C               write(6,*)' bond between i j ebij ',i,j,e
c@@@
                  g=g + e
                  if (g .gt. (g_ref*2.5d0)) return 
                  
                  if (r .gt. 1.d-6) then
                     s=1.d0/r
                     df=2.d0*df*s
                     dg(1,i) = dg(1,i) + rx*df
                     dg(2,i) = dg(2,i) + ry*df
                     dg(3,i) = dg(3,i) + rz*df
                     dg(1,j) = dg(1,j) - rx*df
                     dg(2,j) = dg(2,j) - ry*df
                     dg(3,j) = dg(3,j) - rz*df
                  endif   

200         continue

100   continue

      return
      end
