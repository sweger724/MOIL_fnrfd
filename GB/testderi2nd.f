      SUBROUTINE testderi2nd(natom)

      implicit none
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SGB.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/SGB2DRV.BLOCK'
c args       
      integer natom

c local variables
      integer  i, j,indx, overlappair_idx,k
      double precision e_gbnp, e_gbnp1, e_gbnp2, e_gbnp3, e_gbnp4,
     &   diag_xx, diag_xy, diag_xz,
     &   diag_yy, diag_yz, diag_zz
      double precision tempx, tempy, tempz, dx, dy, dz 
      double precision anal_1st1, anal_1st2, num_2nd, d1iff
      double precision gbnp1, gbnp2, num_1st, anal_2nd
c begin
      dx = 0.00000001d0
      dy = 0.00000001d0
      dz = 0.00000001d0


      do i=1, natom
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xx for atom", i  
         tempx = coor(1,i)
         coor(1,i) = tempx + dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(1,i)
         coor(1,i) = tempx - dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(1,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dx) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dx)

         coor(1,i) = tempx 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 1   
         diag_xx = diag(indx) 
         write(*,1000) num_1st, dpot(1,i),num_1st-dpot(1,i)
1000	 format(1x,'xx> 1st> num anal diff',3(1x,e10.4))
         write(*,1001)num_2nd,diag_xx,num_2nd-diag_xx
1001	 format(1x,'xx> 2nd> num anal diff ',3(1x,e10.4))
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xy for atom", i  
         tempy = coor(2,i)
         coor(2,i) = tempy + dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(1,i)
         coor(2,i) = tempy - dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(1,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dy) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dy)

         coor(2,i) = tempy 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 2   
         diag_xy = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,1002)num_2nd,diag_xy,num_2nd-diag_xy
1002	 format(1x,'xy> 2nd> num anal diff ',3(1x,e10.4))
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_xz for atom", i  
         tempx = coor(1,i)
         coor(1,i) = tempx + dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(3,i)
         coor(1,i) = tempx - dx
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(3,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dx) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dx)

         coor(1,i) = tempx 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 3   
         diag_xz = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,1003)num_2nd,diag_xz,num_2nd-diag_xz
1003	 format(1x,'xz> 2nd> num anal diff ',3(1x,e10.4))
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_yy for atom", i  
         tempy = coor(2,i)

         coor(2,i) = tempy + dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(2,i)

         coor(2,i) = tempy - dy
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(2,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dy) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dy)

         coor(2,i) = tempy 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 4   
         diag_yy = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,1004)num_2nd,diag_yy,num_2nd-diag_yy
1004	 format(1x,'yy> 2nd> num anal diff ',3(1x,e10.4))
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_yz for atom", i  
         tempz = coor(3,i)
         coor(3,i) = tempz + dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(2,i)
         coor(3,i) = tempz - dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(2,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dz) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dz)

         coor(3,i) = tempz 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 5   
         diag_yz = diag(indx) 
c         write(*,*) "xx> 1st> num, anal", num_1st, dpot(1,i)
c         write(*,*) "xx> 1st> ", num_1st - dpot(1,i)
         write(*,1045)num_2nd,diag_yz,num_2nd-diag_yz
1045	 format(1x,'yz> 2nd> num anal diff ',3(1x,e10.4))
c --------------------------------------------------------------------- 
         write(*,*) " ************************"  
         write(*,*) "atom> check diag_zz for atom", i  
         tempz = coor(3,i)

         coor(3,i) = tempz + dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom,gbnp1)
         anal_1st1 = dpot(3,i)

         coor(3,i) = tempz - dz
         do k=1, natom
            dpot(1,k) = 0.0d0
            dpot(2,k) = 0.0d0
            dpot(3,k) = 0.0d0
         end do
         call egb_nonpol2(natom, gbnp2)
         anal_1st2 = dpot(3,i)
         num_2nd = ( anal_1st1 - anal_1st2 ) / (2.0d0 * dz) 
	 num_1st = ( gbnp1 - gbnp2) / ( 2.0d0 * dz)

         coor(3,i) = tempz 
         do k = 1, 6 * natom
            diag(k) = 0.0d0
         end do
	 call egb_np2nd(natom)
         indx = 6*(i-1) + 6   
         diag_zz = diag(indx) 
         write(*,*) "zz> 1st> num, anal", num_1st, dpot(3,i)
         write(*,*) "zz> 1st> ", num_1st - dpot(3,i)
	 write(*,1006)num_1st,dpot(3,i),num_1st-dpot(3,i)
1006	 format(1x,'zz> 1st> num anal diff ',3(1x,e10.4))
         write(*,1005)num_2nd,diag_zz,num_2nd-diag_zz
1005	 format(1x,'zz> 2nd> num anal diff ',3(1x,e10.4))

      end do
c ---------------------------------------------------------------------
      write(*,*) "** testing off-diagonal term" 
      do i = 1, natom-1
         do j = i+1, natom
             indx = overlappair_idx(i,j)
c@		write(*,*) 'natom i j indx ',natom,i,j,indx
             if ( indx .gt. 0 ) then
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xx term"
                tempx = coor(1, j)
                coor(1,j) = tempx + dx
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(1,j) = tempx - dx
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dx)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dx) 
                
                coor(1,j) = tempx
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+1)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1007)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1007		format(1x,'xxo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xy term"
                tempy = coor(2, j)
                coor(2,j) = tempy + dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(2,j) = tempy - dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(2,j) = tempy
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
c		write(*,*)' before calling analytical off d'
c		write(*,*) ' anal = ',offdiaggbnp(6*(indx-1)+2)
c		write(*,*) ' 6*(indx-1) 180*npt ',6*(indx-1),180*natom
c		offdiaggbnp(6*(indx-1)+2) = 0.d0
                call egb_np2nd(natom)
c		write(*,*)' after calling analytical off d'
c		write(*,*) ' anal = ',offdiaggbnp(6*(indx-1)+2)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+2)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1008)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1008		format(1x,'xyo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
c --------------------------------------------------------------------- 	
	        write(*,*) "* off-diagonal xz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(1,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(1,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dz)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dz) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+3)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1009)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1009		format(1x,'xzo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal yy term"
                tempy = coor(2, j)
                coor(2,j) = tempy + dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(2,i)
                 
                coor(2,j) = tempy - dy
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(2,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(2,j) = tempy
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+4)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1010)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1010		format(1x,'yyo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal yz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(2,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(2,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dz)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dz) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+5)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1011)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1011		format(1x,'yzo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
c --------------------------------------------------------------------- 
	        write(*,*) "* off-diagonal zz term"
                tempz = coor(3, j)
                coor(3,j) = tempz + dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp1)
                anal_1st1 = dpot(3,i)
                 
                coor(3,j) = tempz - dz
	        do k=1, natom
		   dpot(1,k) = 0.0d0
		   dpot(2,k) = 0.0d0
		   dpot(3,k) = 0.0d0
	        end do
                call egb_nonpol2(natom, gbnp2)
                anal_1st2 = dpot(3,i)
                
                num_1st = ( gbnp1 - gbnp2 ) / (2.0d0 * dy)
                num_2nd = (anal_1st1 - anal_1st2)/ (2.0d0 * dy) 
                
                coor(3,j) = tempz
                do k = 1, 6 * natom
                   diag(k) = 0.0d0
                end do
                call egb_np2nd(natom)
c                call egb_nonpol2(natom, gbnp2)
                anal_2nd = offdiaggbnp(6*(indx-1)+6)
c                write(*,*) "xxo 1st>i,j, num_1st, anal_1st",i,j,
c     &                      num_1st, dpot(1,j) 
c                write(*,*) "xxo 1st> diff:", num_1st - dpot(1,j) 
		write(*,1012)i,j,num_2nd,anal_2nd,num_2nd-anal_2nd
1012		format(1x,'zzo 2nd> i j num anal diff',
     1            2(1x,i7),3(1x,e10.4))
	     end if 
	 end do 
      end do
c --------------------------------------------------------------------- 

c      write(*,*) "DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
c      call egb_np2nd(natom)

999   return
      end
