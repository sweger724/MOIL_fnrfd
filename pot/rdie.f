      subroutine rdie()

C     rdie eps(r)=1/(eps0*r)
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      double precision pick 
      double precision epstmp
      double precision rx,ry,rz,r2,s2,a,b,e1,e2,q,df,df1,df2,ai,bi
      double precision s,s6,tmp
      double precision qi
      double precision xi,yi,zi,dxi,dyi,dzi
c      
      integer i,j,k,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
      integer ptbeg,ptend
c
      if (eelyes) then
         epstmp = kofdie/eps
      else
         epstmp =  0.0d0
      end if
      if (nocut) then
         cutvdw2 = 10000.d0
         cutele2 = 10000.d0
      end if
      
      if (evdyes) then
         pick = 1.0d0
      else
         pick = 0.0d0
      end if
      

c -------------------------------------------------------------
c yael

      if (prll_on_off) then	
         ptbeg = poipt(monp(my_pe))+1
         if (my_pe.eq.(num_pes-1)) then
            ptend = npt-1
         else
            ptend = poipt(monp(my_pe+1))
         endif   
         if (my_pe.eq.0) then
            jbeg1 = 1
            jbeg2 = 1
            jbeg3 = 1
         else
            jbeg1=point1(ptbeg-1)+1
            jbeg2=point2(ptbeg-1)+1
            jbeg3=point3(ptbeg-1)+1
         endif
      else
         ptbeg = 1
         ptend = npt-1
         jbeg1=1
         jbeg2=1
         jbeg3=1
      end if	
c     
c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
c


      do 400 i=ptbeg,ptend
         
         
         jend1 = point1(i)
         tmp = 1.d0/ptwei(i)
         ai  = epsgm12(i)*pick
         bi  = epsgm6(i)*pick
         qi  = ptchg(i)*epstmp
         xi  = coor(1,i)
         yi  = coor(2,i)
         zi  = coor(3,i)
         dxi = 0.d0
         dyi = 0.d0
         dzi = 0.d0
         
         if (jbeg1.le.jend1) then
            
            do 200 k=jbeg1,jend1
               j=list1(k)
               rx = xi - coor(1,j)
               ry = yi - coor(2,j)
               rz = zi - coor(3,j)
               r2=rx*rx+ry*ry+rz*rz
               s2=1.0d0/r2
c     
c               if (r2.lt.1.0)write(*,*)'problems at',i,j
c     
               q = qi*ptchg(j)
               s = dsqrt(s2)
c               
               if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
                  df1 = 0.d0
                  if ((lesid(i).ne.0) .and.
     *                 (lesid(i) .eq.lesid(j))) then
                     q = q*tmp
                  end if
                  e2  = q*s2
                  df2 = -2.d0*e2*s2
                  go to 100
               end if
               
               a= ai*epsgm12(j)
               b= bi*epsgm6(j)
               if ((lesid(i).ne.0) .and.
     *              (lesid(i) .eq.lesid(j)))  then
                  a = a*tmp
                  b = b*tmp
                  q = q*tmp
               end if
               
               e2  = q*s2
               df2 = -2.d0*e2*s2
               s6=s2*s2*s2
               a = a*s6*s6
               b = b*s6
               e1 = a - b
               df1 = -6.0d0*s2*(a+e1)
               e_vdw = e_vdw + e1 
               
c the electrostatic part should be done without buffering
 100           continue
               df = df1 + df2
               
               rx = df*rx
               ry = df*ry
               rz = df*rz
               dxi = dxi + rx
               dyi = dyi + ry
               dzi = dzi + rz
               dpot(1,j) = dpot(1,j) - rx
               dpot(2,j) = dpot(2,j) - ry
               dpot(3,j) = dpot(3,j) - rz
               e_el = e_el + e2
               
 200        continue
         end if
         jbeg1 = jend1 + 1
         
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c


         jend2 = point2(i)
         
         
         if (jbeg2.le.jend2) then
            do 205 k=jbeg2,jend2
               j=list2(k)
               rx = xi - coor(1,j)
               ry = yi - coor(2,j)
               rz = zi - coor(3,j)
               r2=rx*rx+ry*ry+rz*rz
               if (r2.gt.cutvdw2) go to 205
c     
c               if (r2.lt.1.0)write(*,*)'problems at',i,j
c     
               a = ai*epsgm12(j)
               b = bi*epsgm6(j)
               if ((lesid(i).ne.0) .and.
     *              (lesid(i) .eq.lesid(j)))  then
                  a = a*tmp
                  b = b*tmp
               end if
               s2=1.0d0/r2
               s6=s2*s2*s2
               a = a*s6*s6
               b = b*s6
               e1 = a - b
               df = -6.0d0*s2*(a+e1)
               
               rx = df*rx
               ry = df*ry
               rz = df*rz
               dxi = dxi + rx
               dyi = dyi + ry
               dzi = dzi + rz
               dpot(1,j) = dpot(1,j) - rx
               dpot(2,j) = dpot(2,j) - ry
               dpot(3,j) = dpot(3,j) - rz
               
               e_vdw = e_vdw + e1 
 205        continue
         end if
         jbeg2 = jend2 + 1
         
         
         
c start THIRD loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c


         jend3 = point3(i)
         

         if (jbeg3.le.jend3) then
            
            do 210 k=jbeg3,jend3
               j=list3(k)
               rx = xi - coor(1,j)
               ry = yi - coor(2,j)
               rz = zi - coor(3,j)
               r2=rx*rx+ry*ry+rz*rz
               if (r2.gt.cutele2) go to 210
c     
c               if (r2.lt.1.0)write(*,*)'problems at',i,j
c     
               q=qi*ptchg(j)
               
c chen
               if ((lesid(i).ne.0) .and.
     *              (lesid(i) .eq.lesid(j)))  q = q * tmp
               s2=1.0d0/r2
               s = dsqrt(s2)
               e2  = q*s2
               df2 = -2.d0*e2*s2
               rx = df2*rx
               ry = df2*ry
               rz = df2*rz
               dxi = dxi + rx
               dyi = dyi + ry
               dzi = dzi + rz
               dpot(1,j) = dpot(1,j) - rx
               dpot(2,j) = dpot(2,j) - ry
               dpot(3,j) = dpot(3,j) - rz
               e_el = e_el + e2
               
 210        continue
         end if
         jbeg3 = jend3 + 1
         
         
         
         dpot(1,i) = dpot(1,i) + dxi
         dpot(2,i) = dpot(2,i) + dyi
         dpot(3,i) = dpot(3,i) + dzi
 400  continue
      
      return 
      end
      














