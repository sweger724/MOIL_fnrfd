       program put_missing 
c
c Reading pdb files, identifying missing hydrogens and placing them in
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/DEBUG.BLOCK'

        character*4 styl1,styl2
        character*4 name
        integer namel

        integer ucon,urcrd,uwcrd

        integer of
        logical find
        character*6 getchar

        integer level,lc
        integer lpst,lpend
        character*5 ctyp

        data ucon,urcrd,uwcrd/3*99/
c
c General initialization
c
        lc     = 5
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

        debug  = .false.
        name   = 'puth'
        namel  = 4
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
        lpst  = 1
        lpend = 1

1       continue
        call rline(name,namel,stdi)
        if (find('acti')) stop
        if (find('debu')) debug = .true.
        if (find('file')) then
         if (find('conn')) then
          ucon = of()
          call rconn(ucon)
         end if
         if (find('rcrd')) then
          if (find('bina')) then
           level = 1
           call alert(name,namel,'Binary crd not supported',24,level)
          end if
          ctyp = getchar('ctyp','pdb  ',lc)
          urcrd = of()
          call getcrd(urcrd,ctyp)
         end if
         call put_mis1()
         if (find('wcrd'))then
          if (find('bina')) then
           level = 1
           call alert(name,namel,'Binary crd not supported',24,level)
          end if
          ctyp = getchar('ctyp','CHARM',lc)
          uwcrd = of()
          call putcrd(uwcrd,'CHARM')
         end if
        end if
        go to 1
        end
