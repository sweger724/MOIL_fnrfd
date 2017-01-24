	subroutine chndyn(temp,tmpr,step,tfac,igrid,npri,nstep,nsvel,
     1		pointr,nselec,dmass,divms,d0,e0,e1,grdcmx,grdcmy,grdcmz,
     2		grdlx,grdly,grdlz,r,v,dv,dvold,fact1,fact2,udata,ucrd,
     3		nwcrd,ntest,icol,newv,irand,debug,scalar,sigma,
     4		gamma,rho,lambda,annl,fixend,nlist)
c
	double precision temp,tmpr,step,tfac
	double precision gamma,rho,lambda,annl
	integer igrid,npri,nstep,nsvel,nselec,nlist
	integer ucrd,udata,nwcrd,ntest,icol,newv,irand
	logical debug,fixend
c
c temp   -  assigned temperature
c tmpr   -  actually calculated temperature (from kinetic energy)
c		using function hotchaf
c annl   - scaling factor for velocity annealing
c step   -  time integration step (in AMKA units)
c tfac   -  conversion factor for time from PS to AMKA units
c igrid  -  number of chain monomers
c npri   -  print some useful(?) data each NPRI steps
c nstep  -  number of integration steps 
c nsvel  -  number of steps between scaling of the velocities
c nselec -  number of selected particles, on which the chain constraints
c		are imposed.
c ucrd   -  unit number of file on which the coordinates (path format)
c		are written.
c udata  -  write info on the run on unit UDATA
c nwcrd  -  write coordinates on ucrd each NWCRD steps
c ntest  -  test the constraints each NTEST steps
c icol   -  use ICOL steps for equilibratio (before collecting data)
c newv   -  assign new velocities each NEWV steps
c irand  - a seed for a random nmber generator
c debug  - if .true. print a LOT of debugging info
c
c
	double precision dmass(*),grdcmx(*),grdcmy(*),grdcmz(*)
	double precision divms(*),grdlx(*),grdly(*),grdlz(*)
	double precision r(*),v(*),dv(*),dvold(*),fact1(*),fact2(*)
	double precision d0(*),e0(*),e1(*)
	double precision scalar(*),sigma(*)
c
c dmass - double precision mass vector (m(i=1,npt)(j),j=1,3*igrid)
c divms - double precision 1/mass vector(1/m(i=1,npt)(j),j=1,3*igrid)
c
c d0         -  vecotor of distances between i,i+1 & i,i+2 pairs
c e0         -  vecotor of the individual energies of the monomers
c e1	     -  work vector for the polymer energy routine (sd)
c grdcm[x-z] -  gradient of center of mass constraints
c grdl[x-z]  -  gradient of infitesimal rotation constraints
c r - coordinates (rx(i=1,npt),ry(i=1,npt),rz(i=1,npt)(j=1,igrid))
c v - velocities  (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c dv- forces      (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c dvold - old forces 
c           (_x(i=1,npt),_y(i=1,npt),_z(i=1,npt)(j=1,igrid))
c fact[1-2]  - constant vectors useful for integration (see chmain)
c
	integer pointr(*)
c
c pointr - a pointer to the selected particles
c          (which are subject to chain const.
c

c And now a LOT of neccessary common blocks
c

C
C  next common blocks are used in the energy calculation.
C  list, debugging, coordinate transfer.. etc
C
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/VELOC.BLOCK'
	include 'COMMON/LINE.BLOCK'
c
c local
	integer i,j,k,npt3,ndegf,smlp,idyna,istp,lend,middle
	integer icons,jcons,kcons
	integer one
	double precision enetot,sener,hotchaf,kinet,vfac,sigmav(6)
	double precision tmp1,tmp2
	double precision tpo(1)

	npt3  = npt*3
	ndegf = igrid*npt3
	middle= igrid/2*npt3
	smlp  = 1
c
c Data for velocities constraints.
c sigmav for the velocities is zero (in contrast to the coordinates)
c
		do 1 i=1,6
		 sigmav(i)=0.d0
1		continue
c
c Make a first energy call. Note that at first a call to nbondm is made
c 
	  j = -2
	  do 15 i=1,npt
		j = j + 3
		coor(1,i) = r(middle+j)
		coor(2,i) = r(middle+j+1)
		coor(3,i) = r(middle+j+2)
15	  continue
          call nbondm()
c
c
c Get the energy of the chain = S and the derivatives of S
c
		call sds(sener,dv,r,d0,e0,e1,nselec,pointr,
     1			npt,ndegf,gamma,rho,lambda,fixend,debug)
c
c orie all the structures according to the middle one
c
	if (nlist .ne. 0 ) smlp=nlist
c
c Start integration loop
c

	do  100  idyna = 1,nstep,smlp
c
c do internal update each inbfrq for the non-bonded AND the images
c (if applicable). A single non-bonded list is generated according
c to the middle structure.
c
	 if (nlist.ne.0) then
	  j = -2
	  do 2 i=1,npt
		j = j + 3
		coor(1,i) = r(middle+j)
		coor(2,i) = r(middle+j+1)
		coor(3,i) = r(middle+j+2)
2	  continue
          call nbondm()
	 end if
	 lend=min(nstep,smlp-1+idyna)
	 do 100 istp = idyna,lend
c
c each npri steps printout some infor
c
		if (istp/npri*npri .eq. istp) then
      		 write(udata,*)' chain energy at step ',istp,' = ',sener
		 write(udata,*)' Stru #    dist(i,i+1) dist(i,i+2)   e '
		 do 3 i=1,igrid
		  if (i.le.igrid-2) then
		   tmp1 = d0(i)/dsqrt(dfloat(npt))
		   tmp2 = dsqrt(d0(igrid-1+i)/dfloat(npt))
		  else if (i.le.igrid-1) then
		   tmp1 = d0(i)/dsqrt(dfloat(npt))
		   tmp2 = 0.d0
		  else
		   tmp1 = 0.d0
		   tmp2 = 0.d0
		  end if
		  write(udata,101)i,tmp1,tmp2,e0(i)
101		  format(1x,i6,1x,3(f10.5,1x))
3		 continue
		end if
c
c each nwcrd step write coordinates in path format on ucrd
c and individual monomer energies on udata
c put down also the monomer potential energy
c
		if (debug) write(*,*)' istp nwcrd icol ',istp,nwcrd,icol
		if (istp/nwcrd*nwcrd .eq. istp 
     1				.and. istp.gt.icol) then
c
c save in unit udata the energies of individual monomers
c
		 write(udata,*)' Writing coordinate at step ',istp
                 k = - npt3
		 do 41 i=1,igrid
			k = k + npt3
                write(ucrd)e0(i),(r(j),j=1,npt3-2,3),
     1                  (r(j),j=2,npt3-1,3),(r(j),j=3,npt3,3)
41		 continue
		end if
c
c copy old potential derivatives to temporary data
c
		do 5 i=1,ndegf
			dvold(i)=dv(i)
5		continue
c
c calculate a "free" (unconstrained) step
c
		if (.not.fixend) then
		 do 6 i=1,ndegf
		  r(i) = r(i) + v(i)*step - dvold(i)*fact1(i)
6		 continue
		else
		 do 61 i=npt3+1,ndegf-npt3
		  r(i) = r(i) + v(i)*step - dvold(i)*fact1(i)
61		 continue
		end if
c
c correct the new position in order to satisfy the constraints
c call Correct Rigid Body Motion (CRBM)
c
		if (.not.fixend) then
		 jcons = -npt3 + 1
		 kcons = -5
		 do 62 icons=1,igrid
		  jcons = jcons + npt3
		  kcons = kcons + 6
		  call crbm(r(jcons),scalar,sigma(kcons),grdcmx(jcons),
     1		   grdcmy(jcons),grdcmz(jcons),grdlx(jcons),
     2		   grdly(jcons),grdlz(jcons),
     2		   divms,npt,igrid,nselec,pointr,debug,udata)
62		 continue
		else
		 jcons = 1
		 kcons = 1
		 do 63 icons=1,igrid-1
		  jcons = jcons + npt3
		  kcons = kcons + 6
		  call crbm(r(jcons),scalar,sigma(kcons),grdcmx(jcons),
     1		   grdcmy(jcons),grdcmz(jcons),grdlx(jcons),
     2		   grdly(jcons),grdlz(jcons),
     2		   divms,npt,igrid,nselec,pointr,debug,udata)
63		 continue
		end if
c
c Get the energy of the chain = S and the forces = -grad S
c at the newly generated coordinates.
c
		call sds(sener,dv,r,d0,e0,e1,nselec,pointr,
     1		 npt,ndegf,gamma,rho,lambda,fixend,debug)

c
c calculate "free" velocities
c
		if (.not.fixend) then
		 do 7 i=1,ndegf
		  v(i) = v(i) - fact2(i)*(dv(i)+dvold(i))
7		 continue
		else
		 do 71 i=1+npt3,ndegf-npt3
		  v(i) = v(i) - fact2(i)*(dv(i)+dvold(i))
71		 continue
		end if
c
c correct the velocities to satisfy the constraints
c
		if (.not.fixend) then
	call crbm(v,scalar,sigmav,grdcmx,grdcmy,grdcmz,grdlx,
     1	  grdly,grdlz,divms,npt,igrid,nselec,pointr,debug,udata)
		else
	 call crbm(v(1+npt3),scalar,sigmav,grdcmx,grdcmy,grdcmz,grdlx,
     1	  grdly,grdlz,divms,npt,igrid-2,nselec,pointr,debug,udata)
		end if
c
c check conservation of energy
c
	if (istp/ntest*ntest .eq. istp ) then
		kinet = 0.d0
		do 8 i = 1,ndegf
			kinet = kinet + dmass(i)*v(i)*v(i)
8		continue
		kinet = kinet/2.d0
		enetot = sener + kinet
		write(udata,*)' at step ',istp,' total energy is ',enetot	
		write(udata,*)' tot. pot. ',sener,' total kine. ',kinet	
	end if
c
c Two options:
c if the run is normal MD, do the following:
c if the temperature is more than 20 degrees different from the
c assigned temperature: scale the velocities
c if the run is simulated annealing, do the following:
c multiply current velocities by the annealing factor
c
	 if (nsvel.ne.0) then
	  if  (istp/nsvel*nsvel.eq.istp) then
	   if (dabs(1.d0-annl).lt.1.d-10) then
		tmpr = hotchaf(v(1),v(1+ndegf/3),v(1+2*ndegf/3),
     1			dmass,ndegf/3,ndegf)
	 	if ((dabs(tmpr-temp).gt.20.d0))then
		  vfac=dsqrt(temp/tmpr)
		  write(udata,*)' *** scaling velocities at step ',istp
		  write(udata,*)' current temperature is ',tmpr
		  write(udata,*)' desired temperature is ',temp
		  write(udata,*)' scaling factor = ',vfac
		  do 9 i=1,ndegf
			v(i)=v(i)*vfac
9		  continue
	 	end if
	   else
		tmpr = hotchaf(v(1),v(1+ndegf/3),v(1+2*ndegf/3),
     1			dmass,ndegf/3,ndegf)
		write(udata,*)' *** Annealing at step ',istp
		write(udata,*)' current temperature is ',tmpr
		write(udata,*)' Annealing factor = ',annl
		do 10 i=1,ndegf
			v(i)=v(i)*annl
10		continue
	   end if
	 end if
	end if
c
c the time to reassign velocities?!
c
        if (istp/newv*newv.eq.istp) then

	 write(udata,*)' *** assigning new velocities at step ',istp
	 write(udata,*)' irand = ',irand
c
c set velocities - initial conditions
c
	do 12 i=1,igrid
	 call velinit(irand,temp,one,tpo)
	 k = (i-1)*npt - 3
	 do 11 j=1,npt
	  k = k + 3
	  v(k+1) = velo(1,j)
	  v(k+2) = velo(2,j)
	  v(k+3) = velo(3,j)
11	 continue
12	continue

c
c correct the newly assigned velocities to satisfy the constraints
c
	 if (.not.fixend) then
	  call crbm(v,scalar,sigmav,grdcmx,grdcmy,grdcmz,grdlx,
     1	 grdly,grdlz,divms,npt,igrid,nselec,pointr,debug,udata)
	 else
       call crbm(v(1+npt3),scalar,sigmav,grdcmx,grdcmy,grdcmz,grdlx,
     1	 grdly,grdlz,divms,npt,igrid-2,nselec,pointr,debug,udata)
	 end if 

	 if (debug) then
       write(udata,*)' reassigning velocities at ',istp,' irand '
     1  ,irand
	  write(udata,103)(v(i),i=1,ndegf)
103	  format(1x,8(f9.4,1x),/)
         end if
        end if
100	continue
	return
	end

