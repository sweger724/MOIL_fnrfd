c *****************************************************************
c parallel multicplication of two vectors. 

      subroutine mult_vec_vec(v1,v2,res,rstart,rend)

      integer rstart,rend,i
      double precision v1(*),v2(*),res
      
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      res = 0.d0

      do 10 i=rstart,rend
         res = res + v1(i)*v2(i)
 10   continue

c sum the results from all the processors.
      if (prll_on_off) call reduce_1(res)

      return
      end

c ***********************************************************************
c multoplication of sparce matrix and vector (in parallel).

      subroutine mult_mat_vec_sparse(mat_val,mat_idx,v,res,size
     1     ,rstart,rend)
      
      double precision mat_val(*),res(*),v(*)
      integer size,mat_idx(*),rstart,rend

      integer i,k

      do 10 i=rstart,rend
         res(i) = mat_val(i)*v(i)
         do 20 k=mat_idx(i),mat_idx(i+1) - 1
            res(i) = res(i) + mat_val(k)*v(mat_idx(k))
 20      continue
 10   continue
      
      return
      end

c***********************************************************************
c check if all the bonds constrians are satisfied. conv = 0 if all
c the constrains are setisfied and 1 if not.

      subroutine check_converge (velo,cooref,ishak1,ishak2,
     *     nshak,dissq,shak_diff,sk_start,sk_end,conv)

      double precision cooref(3,*)
      double precision velo(3,*),dissq(*),shak_diff(*)
      integer ishak1(*),ishak2(*),nshak,conv,sk_start,sk_end

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'

      integer iat1,iat2,k,all(maxpe)
      double precision rx,ry,rz,r2

      do 200 k=sk_start,sk_end
         iat1 = ishak1(k)
         iat2 = ishak2(k)         
         rx = cooref(1,k) + velo(1,iat1)-velo(1,iat2)
         ry = cooref(2,k) + velo(2,iat1)-velo(2,iat2)
         rz = cooref(3,k) + velo(3,iat1)-velo(3,iat2)
         r2   = rx*rx + ry*ry + rz*rz - dissq(k)
         if (dabs(r2) .gt. shak_diff(k)) then
            conv = 1
            goto 300
         endif
 200  continue

      conv = 0

 300  continue

c if it is a parallel run, gather the data about convergance from
c all the processors.

      if (prll_on_off) then
         call gather_int(conv,all)
         do 10 k=1,num_pes
            if (all(k).gt.0) then
               conv = 1
               goto 400
            endif
 10      continue
      endif

 400  continue
      
      return
      end
           
c **********************************************************************
      subroutine conjugate_grad_shakept (mat_idx,mat_val)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/SHAKE.BLOCK'

      double precision sigma(maxshak),lam(maxshak)
      double precision rx,ry,rz,t1,t2,mat_val(maxshak*10)
      integer mat_idx(maxshak*10)
      integer iat1,iat2,k,i,j,conv

c             do 999 i=1,nshak
c         write(*,*)'matval,idx',i,mat_val(i),mat_idx(i)
c 999     continue

      do 77 i=1,nshak
         lam(i) = 0.d0
 77   continue

      do 70 k=1,50
c sigma is the constrain vector

         do 50 i=sk_start,sk_end
            iat1 = ishak1(i)
            iat2 = ishak2(i)
            rx = cooref(1,i) + velo(1,iat1) - velo(1,iat2)
            ry = cooref(2,i) + velo(2,iat1) - velo(2,iat2)
            rz = cooref(3,i) + velo(3,iat1) - velo(3,iat2)          
            sigma(i) = (rx*rx + ry*ry + rz*rz - dissq(i))*sqri(i)
 50      continue

c                     if (my_pe .eq. 0) then
c            do 9999 i=1,3
c               write(*,*)'lam(i),sigma(i)'
c     *                   ,lam(i),sigma(i)
c 9999       continue
c               end if

         call conjugate_grad (mat_val,mat_idx,nshak,maxshak,sigma,
     1        lam,cg_shak_it)

 

         do 60 i=sk_start,sk_end
            iat1 = ishak1(i)
            iat2 = ishak2(i)
            t1 = lam(i)*sqri(i)*invms(iat1)
            t2 = lam(i)*sqri(i)*invms(iat2)
            velo(1,iat1) = velo(1,iat1) + t1*cooref(1,i)
            velo(2,iat1) = velo(2,iat1) + t1*cooref(2,i)
            velo(3,iat1) = velo(3,iat1) + t1*cooref(3,i)
            velo(1,iat2) = velo(1,iat2) - t2*cooref(1,i)
            velo(2,iat2) = velo(2,iat2) - t2*cooref(2,i)
            velo(3,iat2) = velo(3,iat2) - t2*cooref(3,i)
 60      continue

         
         if (prll_on_off) then
            do 80 j=1,updates_num
               i = updates(j)
               iat1 = ishak1(i)
               iat2 = ishak2(i)
               t1 = lam(i)*sqri(i)*invms(iat1)
               t2 = lam(i)*sqri(i)*invms(iat2)
               velo(1,iat1) = velo(1,iat1) + t1*cooref(1,i)
               velo(2,iat1) = velo(2,iat1) + t1*cooref(2,i)
               velo(3,iat1) = velo(3,iat1) + t1*cooref(3,i)
               velo(1,iat2) = velo(1,iat2) - t2*cooref(1,i)
               velo(2,iat2) = velo(2,iat2) - t2*cooref(2,i)
               velo(3,iat2) = velo(3,iat2) - t2*cooref(3,i)
 80         continue
         endif

c          if (my_pe .eq. 0) then
c             do 999 i=1,30
c                write(*,*)'velo(1,i)',velo(1,i)
c 999            continue
c             end if

 

         call check_converge (velo,cooref,ishak1,ishak2,nshak,
     1        dissq,shak_diff,sk_start,sk_end,conv)
         
         if (conv.eq.0) goto 444 

 70   continue

 444  if (my_pe.eq.0) then
c         write(*,*)'cg_shake converged after ',k,' steps'
       endif

c            stop
c             return
c      k = pt_start-1 
c      do 45 i=pt_start,pt_end
c         velo(1,i-k)  = velo(1,i) 
c         velo(2,i-k)  = velo(2,i) 
c         velo(3,i-k)  = velo(3,i)
c 45   continue

c      call gather_velo(velo,disp,npts)

      return 
      end

c ********************************************************************
      subroutine conjugate_grad (mat_val,mat_idx,size,nmax,b,x,maxiter)

      integer maxiter,size,nmax,mat_idx(*)
      double precision b(*),x(*),mat_val(*)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/SHAKE.BLOCK'

      integer i,iter
      double precision r(maxshak),p(maxshak),temp(maxshak)
      double precision p1,p2,w,fact1,itol

c get the initial gradient 
      itol = 1.d-9   


      call mult_mat_vec_sparse(mat_val,mat_idx,x,temp,size
     1     ,sk_start,sk_end)
   
      do 11 i=sk_start,sk_end
         r(i) = temp(i) + b(i)
         p(i) = -r(i)
 11   continue

c      do 999 i=1,30
c         write(*,*)'p(i)',p(i),i
c 999     continue


      call mult_vec_vec(r,r,p1,sk_start,sk_end)

      do 100 iter=1,maxiter

         if (prll_on_off) then
            call update_vector(p)
         endif

         call mult_mat_vec_sparse(mat_val,mat_idx,p,temp,size
     1        ,sk_start,sk_end)

         
         call mult_vec_vec(p,temp,w,sk_start,sk_end)     
    

         fact1 = p1/w
         do 20 i=sk_start,sk_end            
            x(i) = x(i) + fact1*p(i)
            r(i) = r(i) + fact1*temp(i)
 20      continue


c         do 9 i=1,5
c            write(*,*)'x(i)',x(i),i
c 9         continue
c       do 8 i=nshak-5,nshak
c           write(*,*)'x(i)',x(i),i
c8         continue


         call mult_vec_vec(r,r,p2,sk_start,sk_end)

          fact1 = p2/p1
         do 30 i=sk_start,sk_end
            p(i) = -r(i) + fact1*p(i)
 30      continue

         p1 = p2 

         
 100  continue
 111  continue

      if (prll_on_off) then
         call update_vector(x)
      endif

      return
      end

c **********************************************************************
      subroutine build_matrix(mat_idx,ishak1,ishak2,nshak)

      integer ishak1(*),ishak2(*),mat_idx(*),nshak

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'

      integer idx,iat1,iat2,iat3,iat4,i,j

      idx = 1 + nshak
      mat_idx(1) = idx + 1
c claculate the symmetric matrix (shake_mat)
      
      do 20 i=1,nshak
        iat1 = ishak1(i)
        iat2 = ishak2(i)         
        do 30 j=1,nshak
           if (i.eq.j) goto 30
           iat3 = ishak1(j)
           iat4 = ishak2(j)    
           if ((iat2.eq.iat3).or.(iat2.eq.iat4).or.
     1          (iat1.eq.iat3).or.(iat1.eq.iat4)) then 
              idx = idx + 1
              mat_idx(idx) = j
           endif
 30     continue
        mat_idx(i+1) = idx+1
 20   continue

      write(*,*)'size !!!',idx

      return
      end

c ***********************************************************************
      subroutine build_matrix_from_idx (mat_val,mat_idx)

      double precision mat_val(*)
      integer mat_idx(*)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      
      integer i,iat1,iat2,iat3,iat4,k,j
      double precision temp,mi

      do 5 i=1,nshak
         iat1 = ishak1(i)
         iat2 = ishak2(i)         
         temp = cooref(1,i)*cooref(1,i) +  cooref(2,i)*cooref(2,i) +
     1         cooref(3,i)*cooref(3,i)
        mat_val(i) = 2*temp*(invms(iat1) + invms(iat2))
        sqri(i) = 1/(dsqrt(mat_val(i)))         
 5    continue


      do 10 i=sk_start,sk_end
         iat1 = ishak1(i)
         iat2 = ishak2(i)         
         mat_val(i) = 1.d0
c         temp = cooref(1,i)*cooref(1,i) +  cooref(2,i)*cooref(2,i) +
c     1         cooref(3,i)*cooref(3,i)
c         mat_val(i) = 2*temp*(invms(iat1) + invms(iat2))
c         sqri(i) = 1/(dsqrt(mat_val(i)))

c build the other elements of the i'th line
         do 20 k=mat_idx(i),mat_idx(i+1) - 1
            j= mat_idx(k)
            iat3 = ishak1(j)
            iat4 = ishak2(j)    
            if ((iat2.eq.iat3).or.(iat2.eq.iat4)) then
               mi = invms(iat2)
            else 
               if ((iat1.eq.iat3).or.(iat1.eq.iat4)) then
                  mi = invms(iat1)
               else                   
                  write(*,*)'error!',i,j
                  write(*,*)iat1,iat2,iat3,iat4
                  goto 20
               endif
            endif      
            mat_val(k) = -2*mi*(cooref(1,i)*cooref(1,j) + 
     1           cooref(2,i)*cooref(2,j) + cooref(3,i)*cooref(3,j))
 
            if ((iat1.eq.iat3).or.(iat2.eq.iat4)) then
               mat_val(k) = -1*mat_val(k)              
            endif            
            mat_val(k) = mat_val(k)*sqri(i)*sqri(j)
 20      continue
 10   continue

c      do 110 i=sk_start,sk_end
c         mat_val(i) = 1.d0
c         do 120 k=mat_idx(i),mat_idx(i+1) - 1
c            j = mat_idx(k)
c            mat_val(k) = mat_val(k)*sqri(i)*sqri(j)
c 120      continue
c 110   continue


      return
      end

c **********************************************************************
      subroutine conjugate_grad_shakevl (mat_idx,mat_val,epsilon)

      double precision epsilon

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/SHAKE.BLOCK'

      double precision sigma(maxshak),lam(maxshak)
      double precision rx,ry,rz,t1,t2,mat_val(10*maxshak)
      integer iat1,iat2,k,i,j,mat_idx(10*maxshak),conv
      
      do 68 i=1,nshak
         lam(i) = 0.d0
 68   continue

      do 70 k=1,50
c sigma is the constrain vector
         do 50 i=sk_start,sk_end
            iat1 = ishak1(i)
            iat2 = ishak2(i)
            rx = velo(1,iat1) - velo(1,iat2)
            ry = velo(2,iat1) - velo(2,iat2)
            rz = velo(3,iat1) - velo(3,iat2)
            sigma(i) = (rx*cooref(1,i) + ry*cooref(2,i) + 
     1           rz*cooref(3,i))*sqri(i)
 50      continue    




         call conjugate_grad (mat_val,mat_idx,nshak,maxshak,sigma,lam,
     1        cg_shakv_it)
       

         do 60 i=sk_start,sk_end
            iat1 = ishak1(i)
            iat2 = ishak2(i)
            t1 = lam(i)*2*sqri(i)*invms(iat1)
            t2 = lam(i)*2*sqri(i)*invms(iat2)
            velo(1,iat1) = velo(1,iat1) + t1*cooref(1,i)
            velo(2,iat1) = velo(2,iat1) + t1*cooref(2,i)
            velo(3,iat1) = velo(3,iat1) + t1*cooref(3,i)
            velo(1,iat2) = velo(1,iat2) - t2*cooref(1,i)
            velo(2,iat2) = velo(2,iat2) - t2*cooref(2,i)
            velo(3,iat2) = velo(3,iat2) - t2*cooref(3,i)
 60      continue

         if (prll_on_off) then
            do 80 j=1,updates_num
               i = updates(j)
               iat1 = ishak1(i)
               iat2 = ishak2(i)
               t1 = lam(i)*2*sqri(i)*invms(iat1)
               t2 = lam(i)*2*sqri(i)*invms(iat2)
               velo(1,iat1) = velo(1,iat1) + t1*cooref(1,i)
               velo(2,iat1) = velo(2,iat1) + t1*cooref(2,i)
               velo(3,iat1) = velo(3,iat1) + t1*cooref(3,i)
               velo(1,iat2) = velo(1,iat2) - t2*cooref(1,i)
               velo(2,iat2) = velo(2,iat2) - t2*cooref(2,i)
               velo(3,iat2) = velo(3,iat2) - t2*cooref(3,i)
 80         continue
         endif 


       
         call check_vl_converge (velo,cooref,ishak1,ishak2,
     *        nshak,1.d-7,conv,sk_start,sk_end)
	              
         if (conv.eq.0) goto 444 
           
 70   continue
      
 444  if (my_pe.eq.0) then
c         write(*,*)'cg_shakevl converged after ',k,' steps' 
      endif
c      k = pt_start-1 
c      do 45 i=pt_start,pt_end
c         velo(1,i-k)  = velo(1,i) 
c         velo(2,i-k)  = velo(2,i) 
c         velo(3,i-k)  = velo(3,i)
c 45   continue

c      call gather_velo(velo,disp,npts)
      return 
      end

c ***********************************************************************

      subroutine check_vl_converge (velo,cooref,ishak1,ishak2,
     *     nshak,epsilon,conv,sk_start,sk_end)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'


      double precision cooref(3,*),epsilon
      double precision velo(3,*)
      integer ishak1(*),ishak2(*),nshak,conv,all(maxpe)
      
      integer iat1,iat2,k,sk_start,sk_end
      double precision rx,ry,rz,r2
     
 

      do 200 k=sk_start,sk_end
         iat1 = ishak1(k)
         iat2 = ishak2(k)
         rx = velo(1,iat1)-velo(1,iat2)
         ry = velo(2,iat1)-velo(2,iat2)
         rz = velo(3,iat1)-velo(3,iat2)
         r2   = rx*cooref(1,k) + ry*cooref(2,k) + rz*cooref(3,k)
         if (dabs(r2) .gt.epsilon) then
            conv = 1
            goto 300
         endif

 200  continue

      conv = 0

 300  continue
c if it is a parallel run, gather the data about convergance from
c all the processors.

      if (prll_on_off) then
         call gather_int(conv,all)
         do 10 k=1,num_pes
            if (all(k).gt.0) then
               conv = 1
               goto 400
            endif
 10      continue
      endif

 400  continue

      return
      end

c ***********************************************************************      
      subroutine get_shared_const (mat_idx)

      integer mat_idx(*)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/SHAKE.BLOCK'

      integer i,j,k,count,left,right
      integer temp(maxpe)
      logical inlist

      count = 0
      if (my_pe.eq.0) then
         left = nshak
      else
         left = shake_bal(my_pe-1)
      endif

      if (my_pe.eq.num_pes-1) then
         right = 0
      else
         right = shake_bal(my_pe+2)-1
      endif

      updates_num = 0
      left_num = 0
      right_num = 0
      total_shared = 0
      left_num_rcv = 0
      right_num_rcv = 0
      do 10 i=sk_start,sk_end
         do 20 k=mat_idx(i),mat_idx(i+1) - 1
            j= mat_idx(k) 
            if (j.lt.sk_start) then
               if (j.ge.left) then
                  if (inlist(left_nbrs_rcv,left_num_rcv,j)) goto 20
                  left_num_rcv = left_num_rcv + 1
                  left_nbrs_rcv(left_num_rcv) = j
c                  write(*,*)'L :',j,ishak1(j),ishak2(j)
               else
                  if (inlist(shared_const,total_shared,i)) goto 20
                  total_shared = total_shared + 1
                  shared_const(total_shared) = i
c                  write(*,*)'G :',i,ishak1(i),ishak2(i)
               endif
            else
               if (j.gt.sk_end) then
                  if (j.le.right) then
                     if (inlist(right_nbrs_rcv,left_num_rcv,j)) goto 20
                     right_num_rcv = right_num_rcv + 1
                     right_nbrs_rcv(right_num_rcv) = j
c                     write(*,*)'R :',j,ishak1(j),ishak2(j)
                  else
                     if (inlist(shared_const,total_shared,i)) goto 20
                     total_shared = total_shared + 1
                     shared_const(total_shared) = i
c                     write(*,*)'G :',i,ishak1(i),ishak2(i)
                  endif
               endif
            endif
 20      continue
 10   continue

      call gather_int(total_shared,shared_num)
      total_shared = 0
      disp1(0) = 0
      do 30 i=1,num_pes
         disp1(i) = disp1(i-1) + shared_num(i-1)
         total_shared = total_shared + shared_num(i-1)
 30   continue
  
      if (total_shared.gt.0) then
         call gather_int1(shared_const,disp1,shared_num)
      endif

      call gather_int(left_num_rcv,temp)
      call getmax(temp,num_pes,mxleft)
      if (my_pe.ne.num_pes-1) right_num = temp(my_pe+2)

      call gather_int(right_num_rcv,temp)
      call getmax(temp,num_pes,mxright)
      if (my_pe.ne.0) left_num = temp(my_pe)

      if (mxright.gt.0) then
         call send_right_int(right_nbrs_rcv,left_nbrs,mxright)
      endif
      if (mxleft.gt.0) then
         call send_left_int(left_nbrs_rcv,right_nbrs,mxleft)
      endif

      do 40 i=1,right_num_rcv
         updates_num = updates_num + 1
         updates(updates_num) = right_nbrs_rcv(i)
 40   continue
      do 50 i=1,left_num_rcv
         updates_num = updates_num + 1
         updates(updates_num) = left_nbrs_rcv(i)
 50   continue
      do 60 i=1,total_shared
         j = shared_const(i)
         if ((j.gt.sk_end).or.(j.lt.sk_start)) then
            updates_num = updates_num + 1
            updates(updates_num) = shared_const(i)
         endif
 60   continue

c      write(*,*)'global left right',total_shared,left_num,right_num
c      write(*,*)'left_rcv right_rcv',left_num_rcv,right_num_rcv 

      return
      end

c logical function : returns true if integer elemant (el) is member
c in the list (list) (of size size).

      logical function inlist(list,size,el)      

      integer list(*),size,el,i
      
      do 10 i=1,size
         if (list(i).eq.el) then 
            inlist = .true.
            return
         endif
 10   continue

      inlist = .false.
      return 
      end

c get the maximum element in the list .      
      subroutine getmax(list,size,max)

      integer list(*),size,max,i

      max = 0
      do 10 i=1,size
         if (list(i).gt.max) max = list(i)
 10   continue

      return
      end

      subroutine update_vector(vec)
      
      double precision vec(*)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      
      integer i,j,idx
      double precision temp(10),temp1(10)

c      write(*,*)right_num,left_num,right_num_rcv,left_num_rcv
c      write(*,*)mxleft,mxright

      if (mxright.gt.0) then
         do 10 i=1,right_num
            j = right_nbrs(i)
            temp(i) = vec(j)
 10      continue

         call send_right_double(temp,temp1,mxleft)

         do 20 i=1,left_num_rcv
            j = left_nbrs_rcv(i)
            vec(j) = temp1(i)          
 20      continue
      endif
      
      if (mxleft.gt.0) then
         do 30 i=1,left_num
            j = left_nbrs(i)
            temp(i) = vec(j)
 30      continue

         call send_left_double(temp,temp1,mxright)

         do 40 i=1,right_num_rcv
            j = right_nbrs_rcv(i)
            vec(j) = temp1(i)    
 40      continue

      endif

      if (total_shared.gt.0) then
         idx = 1
         do 50 i=disp1(my_pe)+1,disp1(my_pe+1)
            j = shared_const(i)
            temp(idx) = vec(j)
            idx = idx + 1
 50      continue

         call gather_double1(temp,disp1,shared_num)
         do 60 i=1,total_shared
            j = shared_const(i) 
            vec(j) = temp(i)            
 60      continue
      endif

      return
      end






