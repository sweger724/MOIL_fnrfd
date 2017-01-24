	subroutine str_shape(inofrz,nofreez,coor,tofi,lambda,error,
     1		delta,s,stdo,uwtab,uwdelta,uws)
c
c compute shape descriptors based on the radius of inertia
c
	integer inofrz,nofreez(*),error
	double precision coor(3,*),tofi(3,3),lambda(3)
	double precision delta,s
	integer stdo,uwtab,uwdelta,uws

c local
	integer i,j,k,l,i1,j1
	double precision tmp,trace_of_t,trace_of_t2,trace_of_t3
	double precision lambda_average
	double precision e(3)

        do 7 k=1,3      
         do 7 l=k,3
           do 6 i=1,inofrz
            do 6 j=i,inofrz
             i1 = nofreez(i)
             j1 = nofreez(j)
             tmp = (coor(k,i1)-coor(k,j1))*(coor(l,i1)-coor(l,j1))
             tofi(k,l)=tofi(k,l)+tmp
6          continue
          tofi(k,l) = 0.5*tofi(k,l)/(inofrz*inofrz)
          tofi(l,k) = tofi(k,l)
7         continue
         trace_of_t = tofi(1,1) + tofi(2,2) + tofi(3,3)
         call house(tofi,3,3,lambda,E,error)
         if (error.ne.0) then
                write(*,*)' error in diagonalization  '
		write(*,*)' error = ',error
                stop
         end if
         lambda_ave = (lambda(1)+lambda(2)+lambda(3))/3.d0
         delta = 0.d0
         s     = 1.d0
         do 8 i=1,3
	   tmp = (lambda(i)-lambda_ave)
           delta = delta + tmp*tmp
           s     = s * tmp
8        continue
         trace_of_t2 = trace_of_t*trace_of_t
         trace_of_t3 = trace_of_t2*trace_of_t
         delta = (3.d0/2.d0)*delta/trace_of_t2
         s     = 27.d0*s/trace_of_t3
         write(stdo,*)'rg = ',sqrt(trace_of_t)
         write(uwtab,*)'tofi = ',((tofi(i,j),i=1,3),j=1,3)
         write(uwdelta,*)'delta = ',delta
         write(uws,*)'s = ',s

	return
	end
