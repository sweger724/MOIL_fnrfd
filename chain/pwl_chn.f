      subroutine pwl_chn(ndegf,minstep,cgrd2,estred,ncalls,
     1  rall,r,dsall,derivs,s,pointr,ipick,igrid,nselec,d0,e0,e1,
     2  grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,scalar,
     3  sigma,sigmag,gamma,rho,lambda,fixend,finished)
        implicit none
c 
c ndegf   - the total number of degrees of freedom for the chain
c minstep - maxmum number of minimization steps
c cgrd2   - convergence criteria. if the current square of the gradient
c           is smaller than cgrd2, we are done.
c estred  - ESTimated REDuction in energy during first step. Employed only
c           in guessing the norm of gradient in the first optimization step
c ncalls  - number of energy/force calls
c rall    - vector of all coordinates including the end points.
c r       - vector of all coordinates EXCLUDING the end points.
c           Note that r and rall overlap. for npt3 number of degrees of
c           freedom r(1)==rall(npt3+1)
c dsall   - chain forces for all degrees of freedom including end points
c derivs  - S derivatives EXCLUDING end points. Similarly to r/rall
c           derivs(i)==dsall(npt3+1)
c s       - The energy of the chain
c pointr  - integer pointer to selected particles
c igrid   - number of polymer segments
c nselec  - number of selected particles
c d0 e0 e1- buffer vectors to caculate certain derivatives, used to store
c           distances and energies of individual monomers
c grdcmx grdcmy grdcmz grdlx grdly grdlz - gradients of the rigid body constraints
c scalar&sigma - where the current and initial values for rigid body constraints
c                are stored
c gamma, rho&lambda   - parameters for the chain energy
c fixend  - logical variable if true the end points of the chain are fixed
c           if false they are free
c 
      logical fixend
      integer ndegf,minstep,ncalls,igrid,nselec
      integer pointr(*),ipick(*)
      integer finished
      double precision cgrd2,estred,s,gamma,rho,lambda
      double precision rall(*),r(*),dsall(*),derivs(*)
      double precision d0(*),e0(*),e1(*),grdcmx(*),grdcmy(*),grdcmz(*)
      double precision grdlx(*),grdly(*),grdlz(*),scalar(*),sigma(*)
      double precision sigmag(*)

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'

c local

c
c INTEGERS
c *** The two variables below intends for internal consistency checks
c mxln   - maximum number of steps in line search
c mxfcn  - maximum number of steps without a reduction in the function value
c isteps  - current iteration number
c ielow - last iteration number that yields a reduction in energy value
c iagain - the most recent restart iteration
c ncalls - internal counter for the number of optimization steps
c nstart  - *** still unknown ***
c nbest  - used for an internal test against ncall. If the Function CHange
c          in the last step (fch) was positive or if fch was zero AND
c          the norm of the new gradient was larger than the norm of the
c          previous gradient DO NOT SET nbest to ncalls
c
      character*6 name
      integer i,j,namel,level,nall,npt3
      integer mxln,mxfcn,isteps,ielow,iagain
      integer nstart,nbest
      integer imore

c DOUBLE PRECISION
c eps    - lowest acceptable step size
c mix    - degree of mixinig required before reinitializing
c beta   - the powell magic paremeter for mixing old and new directions
c        i.e. the way the conujugation is made
c dspl - the second derivative of the spline interpolation
c redene - an estimate for the expected reduction in the energy at the next
c          step
c ediff    - calculated difference between the previously lowest
c               and current function
c einit  - the value of the energy before entering the conjugation loop
c emin   - the lowest value of function found till now
c gamden - unknown yet given by GMIN-GINIT where GINIT is the scalar
c          product between -dx1,-dy1,-dz1 and dx,dz,dz at the beginning
c          the current search direction and GMIN is update to GNEW
c gama   - another of Powell magic parameters
c grad0  - Initial gradient overlap in the main loop
c gnew   - is the overlap of the new gradient with the old one
c drvspl  - estimated nor of gradient for the spline function
c gsup2  - norm of the current gradient
c steplow - lower bound on the step size
c mindr  - the recent value of step size than minus change of step gave the
c          lowest value
c dr - the change in the current step length
c step   - the step size up to the step modification dr
c sum    - the gradient squared
c
      double precision epsilon,mix,beta,dspl,redene,ediff,einit,
     1  emin,gamden,gama,grad0,glow,gnew,drvspl,gsup2,
     2  steplow,step,dr,mindr,sum,temp
c temporary storage array for the history of the optimization
c
      double precision dx1(lgrid*3*maxpt)
      double precision dx1rest(lgrid*3*maxpt)
      double precision dx2rest(lgrid*3*maxpt)
      double precision dxinit(lgrid*3*maxpt)
      double precision xbest(lgrid*3*maxpt)
      double precision dxbest(lgrid*3*maxpt)
      double precision divms(maxpt)

c@
c       double precision sold,dr_test,dstest(lgrid*3*maxpt)
c       double precision delta
c       integer k
c@

        name    = 'powell'
        epsilon= 1.d-12
        namel   = 6
        isteps  = 0                                                     
        ielow   = 0
        iagain  = 0
        mxln    = 100
        mxfcn   = 40
        npt3    = 3*npt
        nall    = ndegf + npt3 + npt3
        mix     = 0.2d0

c------------------------------------------------------------------
c
c check that minstep is not too small
c
        if (minstep.lt.3) then
         level = 0
         call alert(name,namel,'# of steps must be > 2',22,level)
         return
        end if
c------------------------------------------------------------------
c
c initialize inverted mass
c
        do 20 i=1,nselec
         j = pointr(i)
         divms(i) = 1.d0/ptms(j)
20      continue

c initialize the optimization
c First cylce
c
        ncalls = 1
c
c call sds to get energy and forces of the chain. Note that since dsall
c and derivs overlap (see opening remarks) then the calculation of dsall
c also compute derivs (derivatives of S). Furthermore, changes in r modify
c the coordinates sent to ds (rall).
c

c@
c compare derivatives numerically and analytically
c
c@
c        delta = 1.d-6
c        do i = 1,ndegf+6*npt
c          dstest(i) = 0.d0
c        end do
c       do i =1,igrid
c       k = 0
c        do j =1,npt
c               k = k + 1
c               rall((i-1)*npt3+k) = 0
c               k = k + 1
c               if (i.eq.2) then
c                        rall((i-1)*npt3+k) = j
c               else
c                        rall((i-1)*npt3+k) = 0
c               end if 
c               k = k + 1
c               rall((i-1)*npt3+k) = 0
c       end do
c       end do

c       do i=3*npt+1,ndegf+3*npt
c       rall(i) = rall(i) + delta
c       call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,ipick,
c     2         npt,nall,gamma
c     1         ,rho,lambda,fixend,debug)
c       sold = s
c       rall(i) = rall(i) - 2*delta
c       call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,ipick,
c     2         npt,nall,gamma
c     1         ,rho,lambda,fixend,debug)
c       dstest(i) = (sold-s)/(2.d0*delta)
c       write(*,*)'delta i sold s',delta,i,sold,s
c       rall(i) = rall(i) + delta
c       end do

        call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,ipick,
     2          npt,nall,gamma
     1          ,rho,lambda,fixend,debug)
c
c       write(*,*)' i dstest ds '
c       do i = 3*npt+1,ndegf+3*npt
c        write(*,*) i,dstest(i),dsall(i)
c       end do
c       stop
c
c project out from the gradient the components corresponding to
c rigid body motions
c
        j = -npt3+1
        do i=1,igrid-2
         j = j + npt3 
         call crbm(derivs(j),scalar,sigmag,grdcmx,grdcmy,grdcmz,
     1          grdlx,grdly,grdlz,divms,npt,igrid,nselec,pointr,
     2          debug,stdo)
        end do

c
c store the negative of the first gradient in temporary storage area
c store also optimum values of coordinates that are now simply the first set
c calculate the square of the gradient(sum)
c
        sum = 0.d0
        do i=1,ndegf
          dx1(i)    = -derivs(i)
          xbest(i)  = r(i)
          dxbest(i) = derivs(i)
          sum       = sum + derivs(i)*derivs(i)
        end do

c-------------------------------------------------------------------
c
c if the norm of the gradient is less than the required accuracy
c return 
c
        if (sum.lt.cgrd2) then
         write(stdo,*)' POWELL: Gradient converged in first step!'
         finished=1
         return
        end if
c-------------------------------------------------------------------

c
c initialize the index of the best call to one
c
        nbest = 1
c
c initialize the lowest function value to the value of the first call
c
        emin = s
c
c initialize the norm of the gradient gsup2 to be sum
c
        gsup2 = sum

c
c initialize the interpolation for energy reduction and
c calculate minimum step size
c
        redene = estred
        mindr  = estred/sum


c*******************************************************************
c
c*** END OF FIRST STEP AND INITIALIZATION PROCEDURE
c*** BEGIN CENTRAL LOOP
c



3       continue

c increment the number of iterations by 1

        isteps = isteps + 1
c
c initialize the line seacrh
c
        einit = s
        grad0 = 0.d0
        do i=1,ndegf
         dxinit(i) = derivs(i)
         grad0     = grad0 + dx1(i)*derivs(i)
        end do

c-----------------------------------------------------------------
        if (grad0.ge.0.0d0) then
         level = 1
         write(stdo,*)' ncalls isteps grad0 ',ncalls,isteps,grad0
         call alert(name,namel,'Minimization is MAXimization',26,level)
         return
        end if
c-----------------------------------------------------------------

        glow   = grad0
        steplow = -999.d0
        nstart  = ncalls
        imore = -1
        dr = dmin1(mindr,dabs(redene/grad0))
        mindr  = 0.d0

5       continue

c check for maximum number of energy calls
        if (ncalls.eq.minstep) then
         write(stdo,*)'POWELL: Number of steps limit'
         return
        end if

c Generate new coordinates

c       write(stdo,*)' step mindr dr ',step,mindr,dr
        step = mindr + dr
        do i=1,ndegf
         r(i) = xbest(i) + dr*dx1(i)
        end do

c@
c@ check gradient
c       dr_test = 1.d-7
c       s = 0.d0
c       do 101 i=1,ndegf
c        rall(npt3+i) = rall(npt3+i) + dr_test
c        call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,npt,nall,gamma
c     1         ,rho,lambda,fixend,debug)
c        s_old = s
c        rall(npt3+i) = rall(npt3+i) - 2.d0*dr_test
c        call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,npt,nall,gamma
c     1         ,rho,lambda,fixend,debug)
c        dxbest(npt3+i) = (s_old - s) / (2.d0*dr_test)
c        rall(npt3+i) = rall(npt3+i) + dr_test
c101    continue
c         call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,npt,nall,gamma
c     1          ,rho,lambda,fixend,debug)
c       do 102 i=1,ndegf
c               j = npt3+i
c               if (dabs(dsall(j)-dxbest(j)).gt.1.d-4) then
c               write(*,*)' npt3+i ds num_ds ',j,dsall(j),dxbest(j)
c               end if
c102    continue
c@@

c
c       Should set project rigid body motion here
c
c Now call sds to get new S value and new forces
c
        ncalls = ncalls + 1
        call sds(s,dsall,rall,d0,e0,e1,nselec,pointr,ipick
     1          ,npt,nall,gamma
     2          ,rho,lambda,fixend,debug)
c
c project out from the gradient the components corresponding for
c rigid body motions
c
        j = -npt3+1
        do 7 i=1,igrid-2
         j = j + npt3
         call crbm(derivs(j),scalar,sigmag,grdcmx,grdcmy,grdcmz,
     1          grdlx,grdly,grdlz,divms,npt,igrid,nselec,pointr,
     2          debug,stdo)
7       continue

c
c gnew: calculate overlap between current gradient and MINUS the previous grad
c sum : gradient square
c
        gnew = 0.d0
        sum  = 0.d0
        do 8 i=1,ndegf
         gnew = gnew + dx1(i)*derivs(i)
         sum  = sum + derivs(i)*derivs(i)
8       continue
c
c compare current energy with the previous minimum value
c
        ediff = s - emin
c       write(stdo,*)' gnew sum ediff ',gnew,sum,ediff
        if (ediff.lt.0 .or.
     1   (ediff.eq.0 .and. gnew/glow.ge.-1.d0 .and. ncalls.gt.1)) then

c---------------------------------------------------------------
c check for convergence
         if (sum.lt.cgrd2) then
          write(stdo,*)' POWELL: Gradient converged at step ',ncalls
          return
         end if
c---------------------------------------------------------------

c save it!
         emin  = s
         gsup2 = sum
         nbest = ncalls
         do 9 i=1,ndegf
          xbest(i)  = r(i)
          dxbest(i) = derivs(i)
9        continue
c setting below nbest to ncalls means that the last call reduce the energy
c i.e. this was the BESTcall

        else if (ediff.eq.0 .and. gnew/glow.lt.-1.d0) then
c---------------------------------------------------------------
c check for convergence
         if (sum.lt.cgrd2) then
          write(stdo,*)' POWELL: Gradient converged at step ',ncalls
          return
         end if
c---------------------------------------------------------------
        end if
c
c calculate the spline coefficient beta
c
        temp  = (ediff+ediff)/dr-gnew-glow
        dspl  = (gnew-glow)/dr
c if the energy was not reduced we need to update steplow (the length of
c step that not reduces the energy
c
        if (ncalls.gt.nbest) then
         steplow = step
        else
         if (glow*gnew.le.0.d0) steplow = mindr
         mindr  = step
         glow   = gnew
         dr = -dr
        end if
        if (ediff.ne.0.d0) dspl = dspl + (temp+temp)/dr
        

c test converegnce of the line search
c convergence of line search is assumed when
c (i) glow=0 OR glow/grad0<mix
c else tested are possibility for new step initialization
c
        if ((glow.ne.0.d0 .and. ncalls.gt.nstart+1 .and.
     1   dabs(glow/grad0).gt.mix .and. ncalls.lt.nbest+mxln)
     2   .or. (glow.ne.0.d0 .and. ncalls.le.nstart+1)) then
          if (steplow.lt.-0.5d0) then
           dr = 9.0d0*mindr
          else
           dr = 0.5d0*(steplow-mindr)
          end if
          drvspl = glow + dr*dspl
          if (glow*drvspl.lt.0.d0) dr = dr*glow/(glow-drvspl)
c retrun to generate a new coordinate set
c
          go to 5
        else if (glow.ne.0 .and. ncalls.le.nstart+1) then
c the line search did not converge yet, calculate new step
          if (steplow.lt.-0.5d0) then
           dr = 9.0d0*mindr
          else
           dr = 0.5d0*(steplow-mindr)
          end if
          drvspl = glow+dr*dspl
          if (glow*drvspl.gt.0.d0) dr = dr*glow/(glow-drvspl)
c retrun to generate a new coordinate set
c
          go to 5

        else if (glow.eq.0. .or. dabs(glow/grad0).le.mix) then
          if (ncalls.ne.nbest) then
c
c set current values to best found so far
c
           s = emin
           do 10 i=1,ndegf
            r(i) = xbest(i)
            derivs(i) = dxbest(i)
10         continue
          end if
          sum = 0.d0
          do 11 i=1,ndegf
           sum = sum + derivs(i)*dxinit(i)
11        continue
          beta = (gsup2-sum)/(glow-grad0)
c
c test that the new search direction can be made downhill (?!)
c
          if (dabs(beta*glow).gt.mix*gsup2) then
c 
c the current search direction CANNOT be made downhill. Make one more try
c
                if (ncalls.lt.nbest+mxln) then
                 if (steplow.lt.-0.5d0) then
                  dr=9.0d0*mindr
                 else
                  dr = 0.5d0*(steplow-mindr)
                 end if
                 drvspl=glow+dr*dspl
                 if (glow*drvspl.lt.0.d0) dr=dr*glow/(glow-drvspl)
c
c back-branch to another trial in generating a new coordinate set
                 go to 5
                else
                 level = 1
                 write(stdo,1000)isteps,beta,glow,mix,gsup2
1000             format(1x,'isteps = ',i5,' beta =',f10.5,' glow = ',
     1                  f10.5,' mix = ',f10.5,' gsup2 = ',f10.5)
                 call alert(name,namel,'Line search fail',16,level)
                end if
          end if
         if (s.lt.einit) ielow = isteps

c----------------------------------------------------------------------
c test if reduction of energy was achieved in allowed number of iterations
c
          if (isteps.gt.ielow+mxfcn) then
                level = 1
                call alert(name,namel,'Mxfcon failed',13,level)
                return
          end if
c----------------------------------------------------------------------
         
c prepare a new estimate for expected reduction in energy
         redene = mindr*grad0
         if (imore.gt.0) then
          level = 1
          call alert(name,namel,'Imore > 0 ',10,level)
         end if
         if ((iagain.eq.0) .or. (isteps-iagain.ge.ndegf)
     1          .or. (dabs(sum).ge.mix*gsup2)) then
                gamden = glow - grad0
                do  13 i=1,ndegf
                 dx1rest(i) = dx1(i)
                 dx2rest(i) = derivs(i) - dxinit(i)
                 dx1(i)     = -derivs(i) + beta*dx1(i)
13              continue
                iagain = isteps
c back-branch to restart the optimization
                go to 3
        else
c
c calculate the gamma condition
c
                gama = 0.d0
                sum = 0.d0
                do 14 i=1,ndegf
                 gama = gama + derivs(i)*dx2rest(i)
                 sum = sum + derivs(i)*dx1rest(i)
14              continue
                gama = gama/gamden
c
c check the new search direction if it is sufficiently down
c
                if (dabs(beta*glow+gama*sum).ge.mix*gsup2) then
c need to restart
                 gamden = glow-grad0
                 do 15 i=1,ndegf
                        dx1rest(i) = dx1(i)
                        dx2rest(i) = derivs(i) - dxinit(i)
                        dx1(i)     = -derivs(i) + beta*dx1(i)
15               continue
                 iagain = isteps
c Do a new iteration
                 go to 3
                else
c modify the search direction
                 do 16 i=1,ndegf
                   dx1(i) = -derivs(i) + beta*dx1(i) + gama*dx1rest(i)
16               continue
c Do a new iteration
                 go to 3
                end if
         endif
        else
         imore = imore + 1
         if (imore.le.0) then
          if (steplow.lt.-0.5d0) then
           dr = 9.0d0*mindr
          else
           dr = 0.5d0*(steplow-mindr)
          end if
          drvspl = glow + dr*dspl
          if (glow*drvspl.lt.0.d0) dr = dr*glow/(glow-drvspl)
          go to 5
         else 
          if (s.lt.einit) ielow = isteps

c----------------------------------------------------------------------
c test if reduction of energy was achieved in allowed number of iterations
c
          if (isteps.gt.ielow+mxfcn) then
                level = 1
                call alert(name,namel,'Mxfcon failed',13,level)
                return
          end if
c----------------------------------------------------------------------
         end if
        end if
        return
        end
