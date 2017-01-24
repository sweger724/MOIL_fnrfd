      subroutine powel(cgrd2,minstep,estred,ncalls)
c 
c minstep - maxmum number of minimization steps
c cgrd2   - convergence criteria. if the current square of the gradient
c           is smaller than cgrd2, we are done.
c estred  - ESTimated REDuction in energy during first step. Employed only
c           in guessing the norm of gradient in the first optimization step
c ncalls  - number of energy/force calls
c 
      integer minstep,ncalls
      double precision cgrd2,estred
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
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
      integer i,j,namel,level
      integer mxln,mxfcn,isteps,ielow,iagain
      integer nstart,nbest
      integer imore

c DOUBLE PRECISION
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
c work   - used as work space
c
      double precision beta,dspl,redene,ediff,einit,emin,gamden,gama,  
     1  grad0,glow,gnew,drvspl,gsup2,steplow,step,dr,
     2  mindr,sum,temp
c temporary storage array for the history of the optimization
c
      double precision dx1(maxpt),dy1(maxpt),dz1(maxpt)
      double precision dx1rest(maxpt),dy1rest(maxpt),dz1rest(maxpt)
      double precision dx2rest(maxpt),dy2rest(maxpt),dz2rest(maxpt)
      double precision dxinit(maxpt),dyinit(maxpt),dzinit(maxpt)
      double precision xbest(maxpt),ybest(maxpt),zbest(maxpt)
      double precision dxbest(maxpt),dybest(maxpt),dzbest(maxpt)

        name    = 'powell'
        namel   = 6
        isteps  = 0                                                           
        ielow   = 0
        iagain  = 0
        mxln    = 10
        mxfcn   = 40

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

c initialize the optimization
c First cylce
c
        ncalls = 1
c
c call eforce to get energy and forces
c
C@
C       if (my_pe.eq.0) then
C       write(*,*)' before calling eforce crd = '
C       do 101 i=1,npt
C               write(31,*)coor(1,i),coor(2,i),coor(3,i)
C101    continue
C       end if
        call eforce()
        if (prll_on_off) call reduce_energies()
c@
C       if (my_pe.eq.0) then
C        write(31,*)' npt = ',npt
C        do 102 i=1,npt
C               write(31,103)dpot(1,i),dpot(2,i),dpot(3,i)
C102     continue
C103     format(1x,' dVxyz ',3f10.5)
C        call wener(31)
C       end if
c
c store the negative of the first gradient in temporary storage area
c store also optimum values of coordinates that are now simply the first set
c calculate the square of the gradient(sum)
c
        sum = 0.d0
        do 1 j=1,npt_par
         i         = prtc_pointer(j)
         dx1(i)    = -dpot(1,i)
         dy1(i)    = -dpot(2,i)
         dz1(i)    = -dpot(3,i)
         xbest(i)  = coor(1,i)
         ybest(i)  = coor(2,i)
         zbest(i)  = coor(3,i)
         dxbest(i) = dpot(1,i)
         dybest(i) = dpot(2,i)
         dzbest(i) = dpot(3,i)
         sum = sum + dpot(1,i)*dpot(1,i) + dpot(2,i)*dpot(2,i)
     1           + dpot(3,i)*dpot(3,i)
1       continue
        if (prll_on_off) call reduce_1(sum)
C@
C       if (my_pe.eq.0) write(31,*)' sum = ',sum
c-------------------------------------------------------------------
c
c if the norm of the gradient is less than the required accuracy
c return 
c
        if (sum.lt.cgrd2) then
         write(stdo,*)' POWELL: Gradient converged in first step!'
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
        emin = e_total
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
        einit = e_total
        grad0 = 0.d0
        do 4 j=1,npt_par
         i         = prtc_pointer(j)
         dxinit(i) = dpot(1,i)
         dyinit(i) = dpot(2,i)
         dzinit(i) = dpot(3,i)
         grad0     = grad0 + dx1(i)*dpot(1,i) + dy1(i)*dpot(2,i)
     1           + dz1(i)*dpot(3,i)
4       continue
        if (prll_on_off) call reduce_1(grad0)
C@
C       if (my_pe.eq.0) write(31,*) ' grad0 = ',grad0
c-----------------------------------------------------------------
        if (grad0.ge.0.0d0) then
         level = 0
         if (my_pe.eq.0) then
           write(stdo,*)' ncalls isteps grad0 ',ncalls,isteps,grad0
           call alert(name,namel,'Search is up',12,level)
         end if
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
         if(ncalls.eq.nbest) return
         e_total=emin
         do 1001 j=1,npt_par
          i         = prtc_pointer(j)
          coor(1,i) = xbest(i)
          coor(2,i) = ybest(i)
          coor(3,i) = zbest(i)
          dpot(1,i) = dxbest(i)
          dpot(2,i) = dybest(i)
          dpot(3,i) = dzbest(i)
1001     continue
         if (prll_on_off) call gather_crd()
         if (prll_on_off) call gather_force()
         write(stdo,*)'POWELL: Number of steps limit'
         return
        end if

c Generate new coordinates
        step = mindr + dr
        do 6 j=1,npt_par
         i         = prtc_pointer(j)
         coor(1,i) = xbest(i) + dr*dx1(i)
         coor(2,i) = ybest(i) + dr*dy1(i)
         coor(3,i) = zbest(i) + dr*dz1(i)
6       continue
        if (prll_on_off) call gather_crd()
c@
C       if (my_pe.eq.0) then
C        write(31,*)' ****************'
C        do 1002 i=1,npt
C         write(31,1003)coor(1,i),coor(2,i),coor(3,i)
C1002    continue
C1003    format(1x,'crd ',3f10.5)
C        write(31,*)' ***************'
C       end if
         
c
c Now call energy to get new value and new forces
c
        ncalls = ncalls + 1
        call eforce()
        if (prll_on_off) call reduce_energies()
c
c gnew: calculate overlap between current gradient and MINUS the previous grad
c sum : gradient square
c
        gnew = 0.d0
        sum  = 0.d0
        do 8 j=1,npt_par
         i    = prtc_pointer(j)
         gnew = gnew + dx1(i)*dpot(1,i) + dy1(i)*dpot(2,i) 
     1          + dz1(i)*dpot(3,i)
         sum  = sum + dpot(1,i)*dpot(1,i) + dpot(2,i)*dpot(2,i)
     1           + dpot(3,i)*dpot(3,i)
8       continue
        if (prll_on_off) call reduce_1(gnew)
        if (prll_on_off) call reduce_1(sum)
c
c compare current energy with the previous minimum value
c
        ediff = e_total - emin
        if (ediff.lt.0) then

c---------------------------------------------------------------
c check for convergence
         if (sum.lt.cgrd2) then
          if (my_pe.eq.0)
     1     write(stdo,*)' POWELL: Gradient converged at step ',ncalls
          return
         end if
c---------------------------------------------------------------

c save it!
         emin  = e_total
         gsup2 = sum
         nbest = ncalls
         do 9 j=1,npt_par
          i        = prtc_pointer(j)
          xbest(i) = coor(1,i)
          ybest(i) = coor(2,i)
          zbest(i) = coor(3,i)
          dxbest(i) = dpot(1,i)
          dybest(i) = dpot(2,i)
          dzbest(i) = dpot(3,i)
9        continue
c setting below nbest to ncalls means that the last call reduce the energy
c i.e. this was the BEST call

        end if
c
c calculate the spline coefficient beta
c
        if (dr.eq.0.d0) then
                level = 0
                call alert(name,namel,'Divide by 0, dr=0',17,level)
        else
                temp   = (ediff+ediff)/dr-gnew-glow
                dspl = (gnew-glow)/dr
        end if
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
c (i) glow=0 OR glow/grad0<0.1
c else tested are possibility for new step initialization
c
        if ((glow.ne.0.d0 .and. ncalls.gt.nstart+1 .and.
     1   dabs(glow/grad0).gt.0.1d0 .and. ncalls.lt.nbest+mxln)
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

        else if (glow.eq.0. .or. dabs(glow/grad0).le.0.1d0) then
          if (ncalls.ne.nbest) then
c
c set current values to best found so far
c
           e_total = emin
           do 10 j=1,npt_par
            i         = prtc_pointer(j)
            coor(1,i) = xbest(i)
            coor(2,i) = ybest(i)
            coor(3,i) = zbest(i)
            dpot(1,i) = dxbest(i)
            dpot(2,i) = dybest(i)
            dpot(3,i) = dzbest(i)
10         continue
          end if
          sum = 0.d0
          do 11 j=1,npt_par
           i   = prtc_pointer(j)
           sum = sum + dpot(1,i)*dxinit(i) + dpot(2,i)*dyinit(i)
     1          + dpot(3,i)*dzinit(i)
11        continue
          if (prll_on_off) call reduce_1(sum)
          beta = (gsup2-sum)/(glow-grad0)
c
c test that the new search direction can be made downhill (?!)
c
          if (dabs(beta*glow).gt.0.1d0*gsup2) then
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
                 level = 0
                 if (my_pe.eq.0)
     1           call alert(name,namel,'Line search fail',16,level)
                end if
          end if
         if (e_total.lt.einit) ielow = isteps

c----------------------------------------------------------------------
c test if reduction of energy was achieved in allowed number of iterations
c
          if (isteps.gt.ielow+mxfcn) then
                level = 0
                call alert(name,namel,'Mxfcon failed',13,level)
                return
          end if
c----------------------------------------------------------------------
         
c prepare a new estimate for expected reduction in energy
         redene = mindr*grad0
c 
c below we test for a restart procedure
c
         if (imore.gt.0) then
          level = 0
          call alert(name,namel,'Imore > 0 ',10,level)
         end if
         if ((iagain.eq.0) .or. (isteps-iagain.ge.3*npt)
     1          .or. (dabs(sum).ge.0.1d0*gsup2)) then
                gamden = glow - grad0
                do  12 j=1,npt_par
                 i          = prtc_pointer(j)
                 dx1rest(i) = dx1(i)
                 dy1rest(i) = dy1(i)
                 dz1rest(i) = dz1(i)
                 dx2rest(i) = dpot(1,i) - dxinit(i)
                 dy2rest(i) = dpot(2,i) - dyinit(i)
                 dz2rest(i) = dpot(3,i) - dzinit(i)
                 dx1(i)     = -dpot(1,i) + beta*dx1(i)
                 dy1(i)     = -dpot(2,i) + beta*dy1(i)
                 dz1(i)     = -dpot(3,i) + beta*dz1(i)
12              continue
                iagain = isteps
c back-branch to restart the optimization
                go to 3
        else
c
c calculate the gamma condition
c
                gama = 0.d0
                sum = 0.d0
                do 13 j=1,npt_par
                 i    = prtc_pointer(j)
                 gama = gama + dpot(1,i)*dx2rest(i) +
     1                  dpot(2,i)*dy2rest(i) + dpot(3,i)*dz2rest(i)
                 sum = sum + dpot(1,i)*dx1rest(i) +
     1                  dpot(2,i)*dy1rest(i) + dpot(3,i)*dz1rest(i)
13              continue
                gama = gama/gamden
                if (prll_on_off) call reduce_1(gama)
                if (prll_on_off) call reduce_1(sum)
c
c check the new search direction if it is sufficiently down
c
                if (dabs(beta*glow+gama*sum).ge.0.1d0*gsup2) then
c need to restart
                 gamden = glow-grad0
                 do 14 j=1,npt_par
                        i          = prtc_pointer(j)
                        dx1rest(i) = dx1(i)
                        dy1rest(i) = dy1(i)
                        dz1rest(i) = dz1(i)
                        dx2rest(i) = dpot(1,i) - dxinit(i)
                        dy2rest(i) = dpot(2,i) - dyinit(i)
                        dz2rest(i) = dpot(3,i) - dzinit(i)
                        dx1(i)     = -dpot(1,i) + beta*dx1(i)
                        dy1(i)     = -dpot(2,i) + beta*dy1(i)
                        dz1(i)     = -dpot(3,i) + beta*dz1(i)
14               continue
                 iagain = isteps
c Do a new iteration
                 go to 3
                else
c modify the search direction
                 do 15 j=1,npt_par
                  i      = prtc_pointer(j)
                  dx1(i) = -dpot(1,i) + beta*dx1(i) + gama*dx1rest(i)
                  dy1(i) = -dpot(2,i) + beta*dy1(i) + gama*dy1rest(i)
                  dz1(i) = -dpot(3,i) + beta*dz1(i) + gama*dz1rest(i)
15               continue
c Do a new iteration
                 go to 3
                end if
         endif
        else
         level = 0
         call alert(name,namel,'Total mess-up, Powel fail',24,level)
         return
        end if
        end
