	subroutine rcon_specl1(urcon)
c
c Read connectivity data, the one stored in
c connectivity common block: CONNECT.BLCOK
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/CONSPECL1.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/MUTA.BLOCK'
	integer urcon,namel,level
	integer i,j,k,k1,l

	character*5 name
	character*4 version
	character*1 TILDA
	data name/'rcon1'/
	data namel/5/

	rewind urcon

	call init_var()

	read(urcon,99)version
99	format(a4)
	if (version.ne.'12.1') then
		write(*,*)' Wrong version number '
		write(*,*)' Expected 12.1  found ',version
		call alert(name,namel,'ERROR IN CON FILE',17,0)
	end if
	read(urcon,100)TILDA
100	format(a1)
	read(urcon,100)TILDA
c totmon npt nb nangl ntors nimp totex totspe lestyp NBULK
	read(urcon,101)stmon1,snpt1,snb0,snang0,stotex1,
     1   stotspe1,lestyp,NBULK
101	format(1x,8i8)
	read(urcon,100)TILDA
        read(urcon,555)(snb1(i),i=1,nmb),(snang1(i),i=1,nmb)
555     format(1x,8i8)
	read(urcon,100)TILDA
	read(urcon,556)(lz14_1(i),i=0,nmb)
556	format(5(1x,i7))
	read(urcon,100)TILDA
	read(urcon,405)(poitype1(i),i=1,nmb)
405	format(10(1x,i7))
c ---------------------------------------
c	print*,'poitype1=',(poitype1(i),i=1,nmb)
c	print*,'snb0=',snb0
c	print*,'nmb=',nmb
c	print*,'snang0=',snang0
c	print*,'NBULK=',NBULK
c -----------------------------------------
	if (NBULK.EQ.0) then
	 level = 0
	 call alert(name,namel,'No molecules assembled',22,level)
	 return
	else if (stmon1.eq.0) then
	 level = 0
	 call alert(name,namel,'No monomers ?!',14,level)
	 return
	else if (snpt1.eq.0) then
	 level = 0
	 call alert(name,namel,'No particles ?!',15,level)
	 return
	end if

	do 1 i=1,NBULK
	 read(urcon,102)BULK(I)
102	 format(a)
1	continue
	read(urcon,100)TILDA
	read(urcon,103)( pbulk(i),i=1,nbulk )
103	format(10(1x,i7))
	read(urcon,100)TILDA
	read(urcon,104)(smoname1(i),i=1,stmon1)
104	format(10(1x,a4))
	read(urcon,100)TILDA
	read(urcon,105)(spoipt1(i),i=1,stmon1)
105	format(10(1x,i7))
	read(urcon,100)TILDA
	read(urcon,100)TILDA

	do 2 i=1,snpt1
       read(urcon,107)j,newpoit(i),spoimon1(i),styp1(i),sptid1(i),
     1slesid1(i),smutaid1(i),
     1sptnm1(i),sptms1(i),sptchg1(i),sepsgm16(i),sepsgm112(i),
     1sptwei1(i)
107	format(1x,i3,1x,i7,1x,i3,1x,i5,1x,i3,1x,i3,1x,i3,1x,a4,1x,f7.3,
     12x,f9.5,1x,e12.5,e12.5,e12.5)
	pt_to_spcl(newpoit(i)) = i
	if (debug) then
	 write(stdo,*) ' at 107'
	 write(stdo,107)j,newpoit(i),spoimon1(i),styp1(i),sptid1(i),
     1slesid1(i),smutaid1(i),
     1sptnm1(i),
     1sptms1(i),sptchg1(i),sepsgm16(i),sepsgm112(i),sptwei1(i)
	end if
2	continue

	if (snb0.gt.0) then
	 read(urcon,100)TILDA
	 read(urcon,100)TILDA
	 do 3 i=1,snb0
	  read(urcon,109)sib11(i),sib12(i),skbond1(i),sreq1(i)
109	  format(2(1x,i7),2(1x,f10.4))
3	 continue
	end if
c==================================================	 
c ==================================================
	 if (debug) write(stdo,*) ' at 107'
	if (snang0.gt.0) then
	  read(urcon,100)TILDA
	  read(urcon,100)TILDA
	  do 4 i=1,snang0
	   read(urcon,111)siangl11(i),siangl12(i),siangl13(i),
     1skangl1(i),sangleq1(i)
111	   format(3(1x,i7),2(1x,f10.5))
	   sangleq1(i) = sangleq1(i)/pi180
4	  continue
	end if
	if(stotex1 .gt.0) then
           read(urcon,100)TILDA
           read(urcon,100)TILDA
           k = 0
           i = 0
           sexc11(0) = 0
7        continue
           read(urcon,117)j,k1
71       continue
           i = i + 1
           if (debug) write(stdo,*)' j k1 = ',j,k1
117        format(2(1x,i7))
             if (i.lt.j) then
              sexc11(i) = k
              go to 71
           else if (i.gt.j) then
              level = 1
              call alert(name,namel,'Fishy exclusion list!',21,level)
           else
              sexc11(j) = k + k1
              read(urcon,118)(sexc12(l),l=k+1,sexc11(j))
118         format(1x,10(i7,1x))
              k = k + k1
           end if
           if (k.ne.stotex1) go to 7
          end if

        if (stotspe1.gt.0) then
         read(urcon,100)TILDA
         read(urcon,100)TILDA
	 do 8 j=1,stotspe1
	  read(urcon,120)sspec11(j),sspec12(j),s1p14(1,j)
     1,s1p14(2,j),s1p14(3,j)
120       format(2(1x,i7),1x,f30.5,2(1x,f15.5))
8	 continue
        end if

	return
	end
