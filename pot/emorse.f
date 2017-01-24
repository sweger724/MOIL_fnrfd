        subroutine emorse(n)
        include 'COMMON/LENGTH.BLOCK'   
        include 'COMMON/ENERGY.BLOCK'  
        include 'COMMON/COORD.BLOCK'   
        include 'COMMON/CONNECT.BLOCK'   
        include 'COMMON/SPECL.BLOCK'   
        
        double precision R,e1,e2,tmpe,tmp,tempV
        double precision tempx,tempy,tempz                         
        integer n,i ,j,level,namel
        character*6 name

        name = 'emorse'
        namel = 6
         if(D(n).le.0. .or. alpha(n).le.0.) then
          level=1
          call alert(name,namel,'Morse parameter Dmor or
     1 alph is 0',16,level)
         end if
         e_morseb(n) = 0.d0
        tote_morsb = 0.d0
         i=imb1(n)
         j=imb2(n)
         tempx=coor(1,i)-coor(1,j)
         tempy=coor(2,i)-coor(2,j)
         tempz=coor(3,i)-coor(3,j)
   
         R=DSQRT(tempx*tempx+tempy*tempy+tempz*tempz)   
         dist(n) = R

         e1=EXP(-alpha(n)*(R-rmeq(n)))
         e2=e1*e1
         tmpe=e2-e1 
         tmp=2*alpha(n)*D(n)*tmpe/R

         tempV=D(n)*(e2-2*e1)
         e_morseb(n)=tempV

         dpot(1,i)=dpot(1,i)-tmp*tempx
         dpot(1,j)=dpot(1,j)+tmp*tempx

         dpot(2,i)=dpot(2,i)-tmp*tempy
         dpot(2,j)=dpot(2,j)+tmp*tempy

         dpot(3,i)=dpot(3,i)-tmp*tempz
         dpot(3,j)=dpot(3,j)+tmp*tempz

        return
        end
        
