c       reparametrization of the string by linear interpolation  
c       stepwise linear interpolation from Luca Einden paper 

        subroutine reparam(x,nselec,npoint,x1,ll,istp) 

        double precision x(3,nselec*200),x1(3,nselec*200) 
        integer i,j,k ,npoint,nselec,ik,istp
        double precision dist(npoint),leng(npoint),s(npoint),ll  
        logical debug

         debug=.false.
         do i=1,npoint
         leng(i)=0.d0
         dist(i)=0.d0
         enddo 

        ll=0.d0 
        do 150 i=2,npoint
          k = (i-1)*nselec
            dist(i) = 0.d0

           do 140 ik=1,nselec                           
            do n=1,3  
             dist(i) = dist(i) + (x(n,ik+k)-x(n,ik+k-nselec))**2.0 

            x1(n,ik+k)=0.d0 

c            write(72,*) 'first',i,ik,x(n,ik,i)
            enddo 
 140     continue
              
         dist(i)=sqrt(dist(i))
      
        ll=ll+dist(i)
        leng(i)=ll
150     continue 



       do i=2,npoint
         s(i)=(i-1)*leng(npoint)/(npoint-1)
c        write(72,*) i, 'dist',dist(i),"s",s(i),'l',ll       
       enddo 
       

c    two ends are fixed
        k = (npoint-1)*nselec 
        do ik=1,nselec
        do n=1,3
        x1(n,ik)=x(n,ik)
        x1(n,ik+k)=x(n,ik+k) 
        enddo 
        enddo 
       
        


        do 200 i=2,npoint-1
         do ij=2,npoint-1 
        if(leng(ij-1).lt.s(i).and.leng(ij).ge.s(i)) ir=ij
         enddo
        m = (i-1)*nselec
        k= (ir-1)*nselec 
         do 201 ik=1,nselec
            do n=1,3 
            x1(n,ik+m)=x(n,ik+k-nselec)+(s(i)-leng(ir-1))
     $*(x(n,ik+k)-x(n,ik+k-nselec))/dist(ir)

	enddo

 201        continue
200     continue 
      
        
       
         
      
        if(debug) then 
        ll=0.d0 
        do 151 i=2,npoint

            dist(i) = 0.d0

         k=(i-1)*nselec
            do 141 ik=1,nselec                           
            do n=1,3  
             dist(i) = dist(i) + (x1(n,ik+k)-x1(n,ik+k-nselec))**2.0 

      
            enddo 
 141     continue
              
         dist(i)=sqrt(dist(i))
      
        ll=ll+dist(i)
        leng(i)=ll
151     continue 

c        do i=2,npoint-1
c        write(73,*) i,dist(i),dist(i)/(s(i)),'l',ll 
c        enddo 
        endif 
111    format(i8,i8,6f12.5)

        
        return 
        end 
