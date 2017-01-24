         program rot
c special program written in order to rotate gramicidin along the X
c axis by 180 degrees and change back the identity of the amino acids to
c make the calculation symmetric
c

	 include 'COMMON/LENGTH.BLOCK'
	 include 'COMMON/COORD.BLOCK'
	 include 'COMMON/CONNECT.BLOCK'
	 include 'COMMON/LINE.BLOCK'
	 include 'COMMON/DEBUG.BLOCK'
	 include 'COMMON/FREEZ.BLOCK'
	 include 'COMMON/CONVERT.BLOCK'
	 include 'COMMON/UNITS.BLOCK'
	   
c Rotate structure according to predetermined angle
c
         integer ucon,urcrd,u1crd,u2crd
         integer of,i,namel
	 double precision xcm,ycm,zcm,tcm,xtmp,ytmp,ztmp
	 character*3 name
	 logical find

	 stdi=5
	 stdo=6
	 npt=0
	 name='rot'
	 namel=3
c  open junk file for rline
c
	 jnkf=25
	 open(unit=jnkf,status='scratch')
c    defalt parameters
1         continue
	  call rline(name,namel,stdi)
	  if (find('acti')) go to 2
	  if (find('debu')) debug = .true.
	  if (find ('file')) then
	    if (find ('conn')) then
	     ucon=of()
c  get connectivity
             call rconn(ucon)
            end if
	    if (find('w1cr')) u1crd = of()
	    if (find('w2cr')) u2crd = of()
	    if (find ('rcrd')) then
	     urcrd=of()
c   read referance structure
             call getcrd(urcrd,'CHARM')
c calculate center of mass
	     xcm = 0.d0
	     ycm = 0.d0
	     zcm = 0.d0
	     tcm = 0.d0
	     do 81 i=1,npt
	     xcm = xcm + ptms(i)*coor(1,i)
	     ycm = ycm + ptms(i)*coor(2,i)
	     zcm = zcm + ptms(i)*coor(3,i)
	     tcm = tcm + ptms(i)
81	     continue
	     xcm = xcm/tcm
	     ycm = ycm/tcm
	     zcm = zcm/tcm
		write(*,*)' Center of mass position is '
		write(*,*)' Xcm Ycm Zcm ',xcm,ycm,zcm
c correct the center of mass for the x1,y1,z1 coordinate set
	     do 82 i=1,npt
	     coor(1,i) = coor(1,i) - xcm
	     coor(2,i) = coor(2,i) - ycm
	     coor(3,i) = coor(3,i) - zcm
82	     continue
c Write coordinates with zeroed center of mass position
c
		call putcrd(u1crd,'CHARM')
c
c Do now the SPECIFIC rotation for gramicidin
c rotate first the protein
c
	     do 83 i=1,157
	     xtmp     = coor(1,i)
	     coor(1,i)     = coor(1,i+157)
	     coor(1,i+157) = xtmp
	     ytmp     = coor(2,i)
	     coor(2,i)     = -coor(2,i+157)
	     coor(2,i+157) = -ytmp
	     ztmp     = coor(3,i)
	     coor(3,i)     = -coor(3,i+157)
	     coor(3,i+157) = -ztmp
83	     continue
	     do 84 i=316,333
	     xtmp     = coor(1,i)
	     coor(1,i)     = coor(1,i+18)
	     coor(1,i+18)  = xtmp
	     ytmp     = coor(2,i)
	     coor(2,i)     = -coor(2,i+18)
	     coor(2,i+18) = -ytmp
	     ztmp     = coor(3,i)
	     coor(3,i)     = -coor(3,i+18)
	     coor(3,i+18) = -ztmp
84	     continue
	     coor(2,315)   = -coor(2,315)
	     coor(3,315)   = -coor(3,315)
	     call putcrd(u2crd,'CHARM')
            end if
          endif
          goto 1
2	  continue
	  stop
	  end
