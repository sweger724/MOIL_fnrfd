      subroutine sto_force(nstruc,debug1)
c     
c     STO mean force calculations
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/FREEZ.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
      include 'COMMON/DYNA.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/SPECL.BLOCK'
      include 'COMMON/CONSTRAN.BLOCK'
      include 'COMMON/TETHER.BLOCK'
      include 'COMMON/VELOC.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c     
      integer nstruc
      double precision dv(maxpt3),fastonsg 
      double precision wktemp,tempstruc   
      double precision tmp,avgforce
c     
      integer istep,k
c     
      integer i,j,istart,iend,n
      integer l,ll,jjx,npt3s,lstart,lplus,lnow,lminus
c
      logical debug1
c
      logical savesymm
      integer nwatsave
c
czva  store  "symm" to restore after STO ! no "symm" for STO
c
      savesymm=esymyes
      nwatsave=nwaters
c
      npt3s=3*npts
      tempstruc=tempi(1)
c     
      jjx=npt3*(nstruc-1)
c     
      ll = -3 + jjx
      do 10 l=1,npt
	 ll = ll + 3
	 coor(1,l) = r(ll+1)
	 coor(2,l) = r(ll+2)
	 coor(3,l) = r(ll+3)
 10   continue
c 
czva initialization for the average forces and second derivatives
c    (should be done anyway even dynamics average is not used !!!)
      do 25 i=1,npt
         do 20 j=1,3
            dpot_ave(j,i) = 0.d0
 20      continue
         do 22 j=1,9
            d2pt_ave(j,i) = 0.d0
 22      continue
 25   continue
c     
c     start fast forces (dynamics simulations) if (freeze=.true.)
c     
      if(freeze) then
c     
c.....set slow_frz= .false. to run DYNA
c     
         slow_frz=.false.
         ebyes=.false. 
         ethyes=.false.
         etoyes=.false.
         eimyes=.false.
         e14el_yes  = .false.
         e14v_yes  = .false.
c     
c     the temperature should be defined in different way( not tempi)
c     
c
         write(stdo,'(/10x,a,i4)') "Structure N=",nstruc
c         write(stdo,*)   " STO_FORCE : slow_frz=",slow_frz
c
         if (esymyes) call squeeze()
         call nbondm()
         if (esymyes) call syminit()
c
         call vec2matav
czva
         if(nstep.lt.1) goto 999
c
         call str_velinit(tempstruc)
c
c
         do 600 istep=1,nstep

c     
c     start one Verlet step !!!
c     
            do 402 k=1,inofrz
               l = nofreez(k)
               velo(1,l) =  velo(1,l)*dt - factor1(l)*dpot(1,l)
               velo(2,l) =  velo(2,l)*dt - factor1(l)*dpot(2,l)
               velo(3,l) =  velo(3,l)*dt - factor1(l)*dpot(3,l)
 402        continue
c     call matrix shake if (shakm)
            if (shakm) then
               call mshakpt(dt,dt2)
            end if
c     
c     Make a coordinate step and prepare the velocity calculation
c     
            tmp = 1.d0/dt
c     
            do 403 k=1,inofrz
               l = nofreez(k)
               coor(1,l) = coor(1,l) + velo(1,l)
               coor(2,l) = coor(2,l) + velo(2,l)
               coor(3,l) = coor(3,l) + velo(3,l)
               velo(1,l) = velo(1,l)*tmp
               velo(2,l) = velo(2,l)*tmp
               velo(3,l) = velo(3,l)*tmp
 403        continue
c
            if(istep.eq.nstep) goto 600            
c   
            if(mod(istep,nlist).eq.0) then
               if (esymyes) call squeeze()
               call nbondm()
               if (esymyes) call syminit()
            endif
c
c            call eforce()
            call vec2matav
c     
c     calculate new velocities
c     
            do 404 k=1,inofrz
               l = nofreez(k)
               velo(1,l) = velo(1,l) - factor2(l)*dpot(1,l)
               velo(2,l) = velo(2,l) - factor2(l)*dpot(2,l)
               velo(3,l) = velo(3,l) - factor2(l)*dpot(3,l)
 404        continue
c     
c     if matrix shake on, correct velocities
            if (shakm) then
               call mshakvl(dt)
            end if
c     
c     Are we still in equilibration period
c     
            if (istep.lt.neqstep) then
               call veleqtemp(tempstruc)
            end if   
c     
c     end one Verlet step !!!
c
 600     continue
c...  end loop over all dynamics
c
c        compute averages from sum !!!
c
         avgforce=0.0d0
	 do 650 l=1,npts
            do 610 j=1,3
               dpot_ave(j,l) = dpot_ave(j,l)/nstep
               avgforce=avgforce+(dpot_ave(j,l)**2)
 610        continue
            do 620 j=1,9
               d2pt_ave(j,l) = d2pt_ave(j,l)/nstep
 620        continue
c            write(stdo,'(i6,a,2x,3E12.3)') 
c     >           l," VEL=",(dpot_ave(j,l),j=1,3)
c            write(stdo,'(10x,3E12.3)') (d2pt_ave(j,l),j=1,9)
 650     continue
c
        avgforce=sqrt(avgforce/(3.0d0*npts))
c
        write(stdo,'(/10x,a,i6,4x,a,f14.3)')
     >  "Dyna end after step : ",nstep,"Av. fast force=",avgforce
c
 999    continue   
        debug=.false.
        if(debug1) then
           call wener(stdo)	
c               do 6030 n=1,npt
c		  write(stdo,'(i4,2x,3f13.8,2x,3E17.9)')
c     >              n,(coor(j,n),j=1,3),(dpot(j,n),j=1,3)
c 6030	       continue 
        endif
c     
c     zva  store back the fast coordinates for the next iteration
c     
	 ll = -3 + npt3s + jjx
	 do 630 l=npts+1,npt
	    ll = ll +3
	    r(ll+1) = coor(1,l)
	    r(ll+2) = coor(2,l)
	    r(ll+3) = coor(3,l)
 630     continue
c     
c.....set slow_frz= .true. to run STO
c     
         slow_frz= .true.
         ebyes = .true. 
         ethyes = .true.
         etoyes = .true.
         eimyes = .true.
         e14el_yes  = .true.
         e14v_yes  = .true.
c no "symm" for STO. flag should be restored after STO
         esymyes = .false.
c no waters in the system if frozen: set to avoid double counting in "eforce"
         nwaters=0
      endif
c     
c     
c     zva   start STO calculations ...
c     
c      write(stdo,'(/10x,a,i4)') "Structure N=",nstruc
c      write(stdo,*)   " STO_FORCE : slow_frz=",slow_frz
c     
      if (esymyes) then
	 call squeeze()
c     
c     zva     previous setting is dangerous if the water is not selected with fast
c     
c     zva	 ll = -3 + npt3s + jjx
c     zva	 do 30 l=npts+1,npt
	 ll = -3  + jjx
	 do 700 l=1,npt
	    ll = ll +3
	    r(ll+1) = coor(1,l)
	    r(ll+2) = coor(2,l)
	    r(ll+3) = coor(3,l)
 700	 continue
      endif
c     
      call nbondm()
      if(esymyes) call syminit()
c     
      call eforce()
c     
      if(debug1) then
         call wener(stdo)	
c         do 6040 n=1,npt
c            write(stdo,'(i8,2(4x,3f12.3))') 
c     1           n,(coor(j,n),j=1,3),(dpot(j,n),j=1,3)
c 6040    continue 
      endif   
c
      ll = -3
c........this loop is over all prts, so dv(wats) is ready except for -1/kT
c     which is anyway removed from the final formula
      do 810 l=1,npt
         ll = ll + 3
         dv(ll+1) = dpot(1,l)+dpot_ave(1,l)
         dv(ll+2) = dpot(2,l)+dpot_ave(2,l)
         dv(ll+3) = dpot(3,l)+dpot_ave(3,l)
 810  continue
c     
      do 860 l=1,npt3s
         lnow     = l + jjx
         lplus    = lnow + npt3
         lminus   = lnow - npt3
         dv(l)    = mdt2(l)*(r(lplus)+r(lminus)-2.d0*r(lnow))
     >        + gamma*(r(lplus)-r(lminus)) + dv(l)
         onsager  = onsager + dv(l)*dv(l)
         donsger(lminus)=donsger(lminus)+2.0d0*(mdt2(l)-gamma)*dv(l)
         donsger(lnow)  =donsger(lnow)  -4.d0*mdt2(l)*dv(l)
         donsger(lplus) =donsger(lplus) +2.0d0*(mdt2(l)+gamma)*dv(l)
 860  continue
c     
c     
      call getvec(dv,sepfast)
c     
      do 880 l=1,npt3s
         lnow    = l + jjx
         donsger(lnow) = donsger(lnow) + 2.0d0*dv(l)
 880  continue
      
c   
      esymyes=savesymm 
      nwaters=nwatsave
c
      if(freeze) then
c     
c.....set slow_frz= .false. to run DYNA
c     
         slow_frz=.false.
         ebyes=.false. 
         ethyes=.false.
         etoyes=.false.
         eimyes=.false.
         e14el_yes  = .false.
         e14v_yes  = .false.
      endif
c
      return
      end	 
