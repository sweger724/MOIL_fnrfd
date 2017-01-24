      subroutine sds(debug1,upcpueps)
c
c.v1.0 (last changed 24/2/97)
c This routine outputs onsager, the current value of the onsager 
c action, and donsger(*), the vector of derivatives of the onsager 
c action with respect to all degrees of freedom.
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/SYMM.BLOCK'
c
      integer i,j,jjx,l,ll,k
      integer lnow,lminus,lplus
      integer start,sof,npt3s
      double precision upcpueps,elapsedtime
      logical debug1
c
c.....initialize slow/fast separation for action if required
c     NOTE: currently fast modes are those associated with waters
c           and waters must be located after slow part in con
c           COUNTERIONS are also assumed to be fast
c
c.....initialize onsager and the derivative of the Onsager action 
c.....to zero.
      npt3s=3*npts
      onsager = 0.d0
      upcpueps = 0.d0
c
      do 1 i=1,degfmax
         donsger(i) = 0.d0
1     continue

c.....We always start two structures in advance and end two structrues
c.....early. Finaly correction to edge effects (first and last chain
c.....parts) will be added at the end of the main loop.
c
      start = 3
      sof   = pseg + 2
c 
c.....Start MAIN  loop on structures
c
      do 6 j=start,sof
c
c........first coordinate of the structure j (-1)
c
         jjx=j
         call sto_force(jjx,debug1)
c
 6    continue
c
c    echange the stored derivatives on the boundaries
c    and added in the boundaries ones
c
      upcpueps=elapsedtime()
      call excheps_chain()
      upcpueps=elapsedtime()-upcpueps
c
c
c   Scaling to some number stoscale 
c
      onsager=stoscale*onsager
      donsmax=0.0d0
      donsmsq=0.0d0
      do 602 j=start,sof
         jjx=npt3*(j-1)
         do 601 l=1,npt3s
            lnow    = l + jjx
            donsger(lnow) = stoscale*donsger(lnow) 
            donsmax=max(abs(donsger(lnow)),donsmax)
            donsmsq=donsmsq+donsger(lnow)**2
 601     continue
 602  continue
c
      donsmsq=sqrt(donsmsq/(npt3s*(sof-start+1)))
c
      return
      end
c



