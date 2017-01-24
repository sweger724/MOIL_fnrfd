c
c   This is a set of function to provide the coefficients and variance for
c   BD dyna.exe program !!!!
c   The theoretical derivation see Allen & Tildesley (chepter 9)
c
c   Author V. Zaloj
c   October 02, 2000
c

      Real*8 Function bd_c0(x)
c
c     For reference see Allen & Tildesley (p. 261)
c
c     c0=exp(-x)
c
      Implicit None
      Real*8 x,z
      z=abs(x)
      bd_c0=exp(-z)
      return
      end

      Real*8 Function bd_c1(x)
c
c     For reference see Allen & Tildesley (p. 261)
c
c     c1=(1-c0)/x
c
      Implicit None
      Real*8 x,z,bd_c0
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_c1=1.0d0-(z/2.0d0)*(1.0d0-(z/3.0d0)*(1.0d0-(z/4.0d0)))
      Else
         bd_c1=(1.0-bd_c0(z))/z
      Endif
      return
      end

      Real*8 Function bd_c2(x)
c
c     For reference see Allen & Tildesley (p. 261)
c
c     c2=(1-c1)/x
c
      Implicit None
      Real*8 x,z,bd_c1
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_c2=0.5d0*(1.0d0-(z/3.0d0)*
     >              (1.0d0-(z/4.0d0)*(1.0d0-(z/5.0d0))))
      Else
         bd_c2=(1.0-bd_c1(z))/z
      Endif
      return
      end

      Real*8 Function bd_sgrr(x)
c
c     For reference see Allen & Tildesley (p. 262)
c
c     (sgm_rr/dt)**2 = (k_B * T/m) bd_sgrr(x)
c
      Implicit None
      Real*8 x,z,u
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_sgrr=2.0d0*z/3.0d0
      Else
         u=exp(-z)
         bd_sgrr=(2.0d0-(3.0d0-4*u+u*u)/z)/z
      Endif
      return
      end

      Real*8 Function bd_sgvv(x)
c
c     For reference see Allen & Tildesley (p. 262)
c
c     (sgm_vv)**2 = (k_B * T/m) bd_sgvv(x)
c
      Implicit None
      Real*8 x,z,u
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_sgvv=2.0d0*z*(1.0d0-z*(1.0d0-2.0d0*z/3.0d0))
      Else
         u=exp(-z)
         bd_sgvv=1.0d0-u*u
      Endif
      return
      end

      Real*8 Function bd_sgrv(x)
c
c     For reference see Allen & Tildesley (p. 262)
c
c     sgm_rv/dt=c_rv*(sgm_rr/dt)*sgm_vv = (k_B * T/m) bd_sgrv(x)
c
      Implicit None
      Real*8 x,z,u
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_sgrv=z*(3.0d0-z*7.0d0*z/3.0d0)
      Else
         u=exp(-z)
         bd_sgrv=(1.0d0-u)**2/z
      Endif
      return
      end

      Real*8 Function bd_crv(x)
c
c     For reference see Allen & Tildesley (p. 262)
c
c     c_rv = sgm_rv/(sgm_rr*sgm_vv) = bd_crv(x)
c
      Implicit None
      Real*8 x,z,bd_sgrv,bd_sgvv,bd_sgrr
      z=abs(x)
      If(z.lt.1.0d-5) then
         bd_crv=sqrt(3.0d0)/2.0d0
      Else
         bd_crv=bd_sgrv(x)/sqrt(bd_sgvv(x)*bd_sgrr(x))
      Endif
      return
      end
c     
      Real*8 Function erfcinv(x) 
c
c     This is a real*8 function for computation for inverse             
c     complementary error function
c     The variance = 1/2 (or sigma=1/sqrt(2))
c     For reference see: G. Lamm and K. Schulten, JCP, 78 (5), p.2731
c                                      
      Implicit None  
      Real*8 one,b1,b2,b3,b4,b5,b6,b7,x,y,z
      Parameter  (one =  1.000000000d0)  
      Parameter  (b1 = 1.778739000d0)
      Parameter  (b2 = 0.802853000d0) 
      Parameter  (b3 = 0.014606000d0)
      Parameter  (b4 = 2.026268000d0) 
      Parameter  (b5 = 0.378538000d0) 
      Parameter  (b6 = 0.003699600d0)
      Parameter  (b7 = 0.693147181d0)
      z=max(x,1.d-10)
      y=Sqrt(b7-Log(z)) 
      erfcinv=y - (((b3*y+b2)*y+b1)/(((b6*y+b5)*y+b4)*y+one))          
      Return 
      End 
