        double precision function vscalar(x,y,n)
        integer n
        double precision x(*),y(*)
c
c scalar product of two vectors x and y
c
        integer i
        
        vscalar = 0.d0
        do 1 i=1,n
         vscalar = vscalar + x(i)*y(i)
1       continue
        return 
        end 

        double precision function vscalar_prll(x,y,from,to)

        integer from,to
        double precision x(*),y(*)
c
c scalar product of two vectors x and y
c
        integer i
        
        vscalar_prll = 0.d0
        do 1 i=from,to
         vscalar_prll = vscalar_prll + x(i)*y(i)
1       continue
        return 
        end 
