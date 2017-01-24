        subroutine set_cutoffs()
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'

c rmax is maintained here for old input (with a single
c cutoff) to work
c
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if


         if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
         if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

         cutvdw2  = cutvdw2*cutvdw2
         cutvbig2 = cutvbig2*cutvbig2
         cutele2  = cutele2*cutele2
         cutebig2 = cutebig2*cutebig2
         if (cutmono2.lt.0) then
                cutmono2 = cutebig2*1.5625d0
         else
                cutmono2 = cutmono2*cutmono2
         end if

        return
        end
