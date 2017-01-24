       program puth 
       implicit none
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

        integer i,j,k,ires_old
        real xtmp,ytmp,ztmp,x,y,z
        real CACAdist
        real threshold

        data ucon,urcrd,uwcrd/3*99/
        data threshold/1.0/
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
         if (find('wcrd'))then
          if (find('bina')) then
           level = 1
           call alert(name,namel,'Binary crd not supported',24,level)
          end if
          ctyp = getchar('ctyp','CHARM',lc)
          uwcrd = of()
         end if
        end if
        if (.not.find('action')) go to 1

         call puth1()

c
c fixing chain id
c
	do i=1,totmon
         if (moname(i).eq.'NTER'
     1      .or. moname(i).eq.'NTRG'
     2      .or. moname(i).eq.'NTRP') then
             do j=poipt(i-1)+1,poipt(i)
              chain(j) = chain(poipt(i)+1)
             end do
         else if (moname(i).eq.'CTER'
     1     .or. moname(i).eq.'CTRG'
     2     .or. moname(i).eq.'CTRP') then
           chain(poipt(i))=chain(poipt(i-2)+1)
         else 
          do j=poipt(i-1)+1,poipt(i)
           if (chain(j).eq.'*') then
            do k=poipt(i-1)+1,poipt(i)
             if (chain(k).ne.'*') then
              chain(j)=chain(k)
              go to 2
             end if
            end do
2           continue
           end if
          end do
         endif
        end do

c 
c do final check on distances between the CA
c
        xtmp = -9999.
        ytmp = -9999.
        ztmp = -9999.
	do i=1,totmon-1
         do j=poipt(i-1)+1,poipt(i)
          if (ptnm(j).eq.'CA') then
            if (xtmp.lt.-9998) then
             xtmp = coor(1,j)
             ytmp = coor(2,j)
             ztmp = coor(3,j)
             ires_old = i
            else
             x = xtmp - coor(1,j)
             y = ytmp - coor(2,j)
             z = ztmp - coor(3,j)
             CACAdist = x*x+y*y+z*z
             CACAdist = sqrt(CACAdist) 
             if (abs(CACAdist-3.8).gt.threshold) then
              if(i-ires_old.eq.1) then
                write(*,*)
                write(*,*)' CA-CA distance ',CACAdist,
     1          ' deviate from 3.8A. Res: ',ires_old,i,
     2          ' ',moname(ires_old),' ',moname(i)
                write(*,*)
              end if
             end if
             xtmp = coor(1,j)
             ytmp = coor(2,j)
             ztmp = coor(3,j)
             ires_old = i
            end if
           end if
          end do
         end do 

        call putcrd(uwcrd,'CHARM')
	stop
        end
