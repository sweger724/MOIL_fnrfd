      subroutine vec2matav
c
czva: Subroutine vec2matav should compute sum of all forces that act on slow
c     subsystem and sum of all second derivatives (numerical computation) !!!
c==================================================
c Warning: eforce call by cdye and csym
c
c  virtual shift of coordinates to compute the second derivatives !!!
c
      double precision dlts2d
      parameter (dlts2d=1.0d-5)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/SSBP.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/EWALD.BLOCK'
      include 'COMMON/METAL.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c    
      character*9 name
      integer namel
      integer i,j,k,l,m,n
      integer level
c
      double precision d2x(4),d2y(4),d2z(4)
      double precision wktmp
c
      data d2x/1.0d-5,-1.0d-5,0.0d0,0.0d0/
      data d2y/0.0d0,1.0d-5,-1.0d-5,0.0d0/
      data d2z/0.0d0,0.0d0,1.0d-5,-1.0d-5/
c
      data name/'vec2matav'/
      data namel/9/
c
c     average forces and second derivatives
c     
c the structure of the second order derivatives are
c   
c    dFx/dx     dFx/dy     dFx/dz    
c
c    dFy/dx     dFy/dy     dFy/dz
c
c    dFz/dx     dFz/dy     dFz/dz
c
c  the matrix should be diagonal but because of numerical computations
c it may not obey the symmetry !!!!
c
      do 300 m=1,3
c
c translate all slow part for m=1,3. The slow part comes to initial
c position for m=4
c
         do 100 l=1,npts
            coor(1,l) = coor(1,l)+d2x(m)
            coor(2,l) = coor(2,l)+d2y(m)
            coor(3,l) = coor(3,l)+d2z(m)
 100     continue
c
      e_total = 0.d0
      e_bond= 0.d0
      e_theta= 0.d0
      e_tors= 0.d0
      e_imp= 0.d0
      e_cnst= 0.d0
      e_el14=0.0d0 
      e_vdw14=0.0d0
c     
      e_vdw = 0.d0
      e_el  = 0.d0
      e_corr = 0.d0
      e_dir =0.d0 
      e_vsym = 0.d0
      e_lsym = 0.d0   
c  
      do 103 i=1,npt
	 dpot(1,i) = 0.d0
	 dpot(2,i) = 0.d0
	 dpot(3,i) = 0.d0
 103    continue

        if (evdyes.or.eelyes) then
           if (ctrue) then
              if (shift) then
                 level = 1
                 call alert(name,namel,'Shift not supported',19,level)
              else if (.not.efyes) then
                 if (.not.arith) then
                    if (.not. 3*nwaters.eq.npt) then
                       call cdie()
                    end if
                 else
                    call cdie_arith()
                 end if
              else if (efyes) then
                 call nbfinit()
              end if
           else
              call rdie()
           end if          
        end if

        if((evdyes.or.eelyes).and.totspe.gt.0.and.(.not.efyes)) then
           if (e14el_yes.or.e14v_yes) then
              call ener14()
              e_total = e_total + e_el14 + e_vdw14
           endif
        end if
c
c     only cdie is supported for symmetry related operations
c     
        if (esymyes) then
           e_vsym = 0.d0
           e_lsym = 0.d0
c
           if (iblock1(symnum).gt.0 .or.
     1          iblock2(symnum).gt.0 .or. iblock3(symnum).gt.0)
     2          call symcdie()
        end if

         if (lcent.and.icenter.gt.0) then
            call ecent()
         end if


c     ---------------------------------
         if (eteth_yes) then
            call etether()
          end if

         if (ehyes) then
            call ehydro()
         end if
c
         k=3*(m-1)
c
         do 200 l=1,npts          
            d2pt_ave(k+1,l)=d2pt_ave(k+1,l)+(dpot(1,l)/dlts2d)        
            d2pt_ave(k+2,l)=d2pt_ave(k+2,l)+(dpot(2,l)/dlts2d)        
            d2pt_ave(k+3,l)=d2pt_ave(k+3,l)+(dpot(3,l)/dlts2d)
 200     continue
c
 300  continue

      do 350 l=1,npts
         coor(1,l) = coor(1,l)+d2x(4)
         coor(2,l) = coor(2,l)+d2y(4)
         coor(3,l) = coor(3,l)+d2z(4)
 350  continue
c     
      call eforce()
c
      do 500 l=1,npts
         do 400 j=1,3
c
c    average forces
c
            dpot_ave(j,l) = dpot_ave(j,l) + dpot(j,l) 
c
c   substract the base position part for second derivatives
c 
            wktmp=dpot(j,l)/dlts2d  
            d2pt_ave(j,l)=d2pt_ave(j,l)-wktmp
            d2pt_ave(j+3,l)=d2pt_ave(j+3,l)-wktmp
            d2pt_ave(j+6,l)=d2pt_ave(j+6,l)-wktmp
 400     continue
 500  continue
c
      return
c
      end







