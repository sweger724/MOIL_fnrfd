       program xtemp
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
        include 'COMMON/OVERLAP.BLOCK'
        character*4 styl1,styl2
        character*5 name
        integer nofreez(maxpt),ipick(maxpt),jpick(maxpt)
        integer nofrez1(maxpt),ipic1(maxpt),jpic1(maxpt)
        integer namel,inofrz,npick,npic1,k
        integer utemp(10),lgrp(10)
        integer ngrp
        double precision rms, tmp_mass(maxpt)
        double precision v2
        double precision kinetic(10)
        integer ucon,urcrd,urvel,ucmbn
        integer geti,of
        logical find,fopen,lcmbn,wsub,ovlp
        integer i,ich,level
        integer lpend,ncmbn,ntotal
        integer rbin

        integer l,j
        
        data ucon,urcrd,urvel/3*99/
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
        norew  = .false. 
        name   = 'ccrd'
        namel  = 4
        wsub = .false.
        ovlp = .false.
        npic1 = 0
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
         end if
         if (find('rcrd')) urcrd = of()
         if (find('rvel')) urvel = of()
        end if
        if (find('pick')) then
                call pick(ipick,ngrp)
		if (ngrp.gt.10) then
                  write(*,*)' max no. of temp group is 10!'
                 stop
                end if
                do i=1,ngrp
                utemp(i) = 70+i
                open(unit=utemp(i))
                end do
        end if  
        lpstr  = geti('lpst',lpstr)
        lpend = geti('lpen',lpend)
        if (find('ovlp')) ovlp = .true.
        if (find('acti')) go to 2
        go to 1
2       continue
c
c check that required files were opened
c
        if (.not.fopen(ucon)) then
         level = 1
         call alert(name,namel,'ucon not opened',15,level)
        else if (.not.fopen(urvel) .and. (.not.lcmbn)) then
         level = 1
         call alert(name,namel,'urvel not opened',16,level)
        end if
          do i=lpstr,lpend
           norew = .true.
           call rdyncrd(urvel,i,inofrz,nofreez,rbin)
	    do j=1,ngrp
             kinetic(j) = 0
             lgrp(j) =0
            end do
            do k=1,npt
             do j=1,ngrp
              if (ipick(k).eq.j) then
               lgrp(j) = lgrp(j) + 1
               v2 = coor(1,k)*coor(1,k)
               v2 = v2 + coor(2,k)*coor(2,k)
               v2 = v2 + coor(3,k)*coor(3,k)
               kinetic(j) = kinetic(j) + ptms(k)*v2*0.5
              end if
             end do
            end do
                do j=1,ngrp
                 write(utemp(j),*)kinetic(j)/lgrp(j)/0.002
                end do
           end do
        stop
        end
