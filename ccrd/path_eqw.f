       program path_eqw
c 
c Creating a triplet for a fp calculation.
c 
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        character*4 styl1,styl2
        character*8 name
        integer nofreez(maxpt),ipick(maxpt),jpick(maxpt)
        integer nofrez1(maxpt),ipic1(maxpt),jpic1(maxpt)
        integer namel,inofrz,npick,npic1
        double precision rms
        integer ucon,urcrd,uwcrd,ucmbn
        integer geti,of
        logical find,fopen,lcmbn,wsub,ovlp
        integer i,j,ich,level
        integer lpend,ncmbn,ntotal

        integer str_indx,all_struc,desired_water

        data ucon,urcrd,uwcrd/3*99/
c
c General initialization
c
        stdi   = 5
        stdo   = 6
        stderr = 0
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
        name   = 'path_eqw'
        namel  = 8
        wsub = .false.
        ovlp = .false.
c
c open scratch file for line manipulation
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c 
c ccrd default
c
        styl1 = 'unkw'
        styl2 = 'unkw'
        lpstr = 1
        lpend = 1
        desired_water = 0
        all_struc = 0

        call rline(name,namel,stdi)
        ucon  = of()
        call rconn(ucon)
        desired_water = nwaters
        close(ucon)

        call rline(name,namel,stdi)
        all_struc = geti('#str',all_struc)
        call rline(name,namel,stdi)
        uwcrd = of()

        do 10 str_indx=1,all_struc
        call rline(name,namel,stdi)
                 ucon  = of()
                 call rconn(ucon)
                 close (ucon)
                 call freeunit(ucon)

                call rline(name,namel,stdi)
                urcrd = of()
                call getcrd(urcrd,'CHARM')
                close (urcrd)
                 call freeunit(urcrd)
c@ 
        write(*,*) ' coordinate set # ',str_indx
        write(*,*)(coor(j,1),j=1,3)
        write(*,*)' desired_water nwaters ',desired_water,nwaters
        if (nwaters.eq.desired_water) then
                call wpath(uwcrd)
        else if (nwaters.gt.desired_water) then
                npt = npt - 3*(nwaters-desired_water)
                call wpath(uwcrd)
        else if (nwaters.lt.desired_water) then
                do 2 i=npt+1,npt+3*(desired_water-nwaters)
                        coor(1,i) = 9999.d0
                        coor(2,i) = 9999.d0
                        coor(3,i) = 9999.d0
2               continue
                npt = npt + 3*(desired_water-nwaters)
                call wpath(uwcrd)
        end if
10      continue
        stop
        end
