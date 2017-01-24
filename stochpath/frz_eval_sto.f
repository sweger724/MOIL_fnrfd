      subroutine frz_eval_sto()
c
c pick freezing atoms and re-assign energy calculations
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'


      integer i,j

      freeze = .true.

      call pick(ipick,i)
c
      inofrz = 0
c
      do 4 i=1,npt
         if (ipick(i).ne.0) then
	    inofrz               = inofrz + 1
	    nofreez(inofrz)      = i
	    zerofrz(i)           = 1
         else
	    zerofrz(i)           = 0
         end if
 4    continue

      sepfast = .true.
c
      slow_frz=.false.
c
      ebyes=.false. 
      ethyes=.false.
      etoyes=.false.
      eimyes=.false.
      e14el_yes  = .false.
      e14v_yes  = .false.
c
      npts=npt-inofrz

C logical list of frozen monomers

      do 7 i=1,totdmon
         frzM(i) = .false.
         do 5 j=dpoipt(i-1)+1,dpoipt(i)
            if (zerofrz(j).eq.1) go to 6
 5       continue
         frzM(i) = .true.
 6       continue
 7    continue

      return
      end	 
