      program test_alabar

c     for alanine dipeptide

      implicit none

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/DEBUG.BLOCK'


      character*10 name
      integer namel, i, j, k, ntest
      integer ucon,urcrd,of,geti,uwgrd
      logical find, fopen
      double precision e0, tors_grad(3,maxpt), d_e, e_i, d_g
      integer phii(4), psii(4)
      double precision getd, stepsize

      name = 'test_alabar'
      stdi   = 5
      stdo   = 6
      stderr = 0
      ntest = 0

c     an unfortunate necessity
      jnkf = 25
      open (unit=jnkf,status='scratch')

c     begin input loop
 1    continue
      call rline(name,namel,stdi)
      if (find('file')) then
         if (find('conn')) then
            ucon = of()
            call rconn(ucon)
         endif

         if (find('rcrd')) urcrd = of() 
         if (find('wgrd')) uwgrd = of() 
      end if
      ntest   = geti('#tes',ntest)
      if (find('acti')) go to 2
      go to 1
 2    continue
c     end input loop


      write(*,*) 'ntest = ', ntest


      read(urcrd) e0, ((coor(j,i),i=1,npt),j=1,3)


c     atoms for psi
      psii(1) = 4
      psii(2) = 6
      psii(3) = 8
      psii(4) = 10


c     first zero the energy/force
      e_total = 0.d0
      do 20 j = 1,npt
         do 200 k = 1,3
            dpot(k,j) = 0.d0
 200     continue
 20   continue


c     compute energy and gradient of initial point
      call alabarrier()

c     save energy and gradient in e_i and tors_grad, resp.
      e_i = e_total
      write(*,*) 'Energy of initial structure = ', e_i
      do 21 j = 1,npt
         do 201 k = 1,3
            tors_grad(k,j) = dpot(k,j)
 201     continue
 21   continue


c     now perturb initial position and compare
      do 10 i = 1,4
         do 30 k = 1,3
            do 100 j = 1,ntest
               stepsize = 0.00001 * j

c              take a step, get new energy, then return
               e_total = 0.d0
               coor(k,psii(i)) = coor(k,psii(i)) + stepsize
               call alabarrier()
               coor(k,psii(i)) = coor(k,psii(i)) - stepsize
               d_e = e_total - e_i
               d_g = stepsize * tors_grad(k,psii(i))
c              write(*,*) 'd_e, d_g: ', d_e, d_g
               write(uwgrd,*) stepsize, d_e - d_g

 100        continue
 30      continue
 10   continue


      end
