       program tempmesh
c 
c Extracting the local temperature of a subset of atoms
c 
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
        double precision getd
        integer mesh
	parameter(mesh=30)
        character*4 styl1,styl2
        character*6 name
        integer nofreez(maxpt),ipick(maxpt),jpick(maxpt)
        integer namel,inofrz,npick,npic1
        integer rbin
        integer ibox(mesh,mesh,mesh)
        integer ix,iy,iz,ncount
        integer ngrp
        real x,y,z
        double precision tcut
        real grid_size
	real tgrid(mesh,mesh,mesh)
        real v2
        double precision xmin,xmax,ymin,ymax,zmin,zmax
        real dx,dy,dz
        integer ucon,urcrd,urvel,uwcrd,utemp
        integer geti,of
        logical find,fopen
        integer i,j,k,l,level
        integer lpend

        
        data ucon,urcrd,urvel,uwcrd,utemp/5*99/
        data tgrid/27000*0./
	data ibox/27000*0/
c
c General initialization
c
        stdi   = 5
        stdo   = 6
        stderr = 0
        rbin = 1

        totmon = 0
        npt    = 0
        nb     = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        nbulk  = 0
        inofrz = 0
        debug  = .false.
        norew  = .true. 
        name   = 't-grid'
        namel  = 6
        xmin = 0.d0
        xmax = 0.d0
        ymin = 0.d0
        ymax = 0.d0
        zmin = 0.d0
        zmax = 0.d0
        tcut = 0.d0
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c 
c
        lpstr = 1
        lpend = 1
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) debug = .true.
        if (find('file')) then
         if (find('conn')) then
                 ucon  = of()
                 call rconn(ucon)
                 inofrz = npt
		 do i=1,npt
                  nofreez(i) = i
                 end do
         end if
         if (find('rcrd')) urcrd = of()
         if (find('rvel')) urvel = of()
	 if (find('wcrd')) uwcrd = of()
         if (find('utem')) utemp = of()
        end if
        if (find('pick')) then
                call pick(ipick,ngrp)
        end if  
        lpstr  = geti('lpst',lpstr)
        lpend = geti('lpen',lpend)
        xmax  = getd('xmax',xmax)
        xmin  = getd('xmin',xmin)
        ymax  = getd('ymax',ymax)
        ymin  = getd('ymin',ymin)
        zmax  = getd('zmax',zmax)
        zmin  = getd('zmin',zmin)
        tcut  = getd('tcut',tcut)
        if (find('acti')) go to 2
        go to 1
2       continue
        dx = (xmax-xmin)/(mesh)
        dy = (ymax-ymin)/(mesh)
        dz = (zmax-zmin)/(mesh)
c
c check that required files were opened
c
        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (.not.fopen(urvel) ) then
         level = 1
         call alert(name,namel,'urvel not opened',16,level)
        end if
          do l=lpstr,lpend
           norew = .true.
           call rdyncrd(urvel,l,inofrz,nofreez,rbin)
           call vdcopy(coor,velo,3*npt)
           call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
           do j =1,npt
            x = coor(1,j)
            y = coor(2,j)
            z = coor(3,j)
            ix = (x-xmin)/dx+1
            iy = (y-ymin)/dy+1
            iz = (z-zmin)/dz+1

            if (ix.le.0) ix=1
            if (ix.gt.mesh) ix=mesh
            if (iy.le.0) iy=1
            if (iy.gt.mesh) iy=mesh
            if (iz.le.0) iz=1
            if (iz.gt.mesh) iz=mesh
            ibox(ix,iy,iz) = ibox(ix,iy,iz) + 3
            if (moname(poimon(j)).eq.'TIP3'.and.ptnm(j).eq.'OH2') then
             ibox(ix,iy,iz) = ibox(ix,iy,iz) - 1
            else if (moname(poimon(j)).eq.'LYS'
     1           .and.ptnm(j).eq.'NZ') then
             ibox(ix,iy,iz) = ibox(ix,iy,iz) - 3
            else if (moname(poimon(j)).eq.'ARG'
     1           .and.ptnm(j).eq.'NH1') then
             ibox(ix,iy,iz) = ibox(ix,iy,iz) - 2
            else if (ptnm(j)(1:1).eq.'H') then
             ibox(ix,iy,iz) = ibox(ix,iy,iz) - 1
            end if
            v2 = velo(1,j)*velo(1,j)
            v2 = v2 + velo(2,j)*velo(2,j)
            v2 = v2 + velo(3,j)*velo(3,j)
            tgrid(ix,iy,iz) = tgrid(ix,iy,iz) + ptms(j)*v2
           end do
          end do
          ncount = 0
          do i=1,mesh
           do j=1,mesh
            do k=1,mesh
             if (ibox(i,j,k).eq.0) then
               write(*,*)' warning: box(',i,j,k,') is zero '
               tgrid(i,j,k) = -1.
             else
              tgrid(i,j,k) = tgrid(i,j,k)/(0.002*ibox(i,j,k))
             end if
             write(utemp,*)i,j,k,tgrid(i,j,k)
             if (tgrid(i,j,k).gt.tcut) then
              ncount = ncount + 1
              coor(1,ncount) = xmin + (i-1)*dx
              coor(2,ncount) = ymin + (j-1)*dy
              coor(3,ncount) = zmin + (k-1)*dz
              ptnm(ncount) = 'XE'
              moname(ncount) = 'XE'
             end if
            end do
           end do
          end do
         npt = ncount
         call putcrd(uwcrd,'CHARM')
        stop
        end
