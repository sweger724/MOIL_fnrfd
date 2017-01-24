	subroutine rcon_specl2(urcon)
c
c Read connectivity data, the one stored in
c connectivity common block: CONNECT.BLCOK
c
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/CONSPECL2.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	integer urcon,namel,level
	integer i,j,k,k1,l

	character*5 name
	character*4 version
	character*1 TILDA
	data name/'rcon2'/
	data namel/5/

	rewind urcon

	call init_var()

        read(urcon,99)version
99      format(a4)
        if (version.ne.'12.1') then
                write(*,*)' Wrong version number '
                write(*,*)' Expected 12.1  found ',version
                call alert(name,namel,'ERROR IN CON FILE',17,0)
        end if

	read(urcon,100)TILDA
100	format(a1)
	read(urcon,100)TILDA
c totmon npt nb nangl ntors nimp totex totspe lestyp NBULK
	read(urcon,101)stmon2,snpt2,snb0,snang0,stotex2,
     1   stotspe2,lestyp,NBULK
101	format(1x,8i8)
        read(urcon,100)TILDA
        read(urcon,555)(snb2(i),i=1,nmb),(snang2(i),i=1,nmb)
555     format(1x,8i8)
        read(urcon,100)TILDA
        read(urcon,556)(lz14_2(i),i=0,nmb)
556     format(5(1x,i7))
	read(urcon,100)TILDA
        read(urcon,405)(poitype2(i),i=1,nmb)
405     format(10(1x,i7))
c ---------------------------------------
c	print*,'snpt1=2',snpt2
c	print*,'snb2=',snb2
c	print*,'nmb=',nmb
c	print*,'snangl2=',snangl2
c	print*,'NBULK=',NBULK
c -----------------------------------------
	if (NBULK.EQ.0) then
	 level = 0
	 call alert(name,namel,'No molecules assembled',22,level)
	 return
	else if (stmon2.eq.0) then
	 level = 0
	 call alert(name,namel,'No monomers ?!',14,level)
	 return
	else if (snpt2.eq.0) then
	 level = 0
	 call alert(name,namel,'No particles ?!',15,level)
	 return
	end if

	do 1 i=1,NBULK
	 read(urcon,102)BULK(I)
102	 format(a)
1	continue
	read(urcon,100)TILDA
	read(urcon,103)(pbulk(i),i=1,nbulk)
103	format(10(1x,i7))
	read(urcon,100)TILDA
	read(urcon,104)(smoname2(i),i=1,stmon2)
104	format(10(1x,a4))
	read(urcon,100)TILDA
	read(urcon,105)(spoipt2(i),i=1,stmon2)
105	format(10(1x,i7))
	read(urcon,100)TILDA
	read(urcon,100)TILDA
	do 2 i=1,snpt2
	read(urcon,107)j,newpoit(i),spoimon2(i),styp2(i),sptid2(i),
     1	slesid2(i),smutaid2(i),
     1	sptnm2(i),sptms2(i),sptchg2(i),sepsgm26(i),sepsgm212(i),
     1  sptwei2(i)
107	format(1x,i3,1x,i7,1x,i3,1x,i5,1x,i3,1x,i3,1x,i3,1x,a4,1x,f7.3,
     1   2x,f9.5,1x,e12.5,e12.5,e12.5)
	if (debug) then
	 write(stdo,*) ' at 107'
	 write(stdo,107)j,newpoit(i),spoimon2(i),styp2(i),sptid2(i),
     1	slesid2(i),smutaid2(i),
     1	 sptnm2(i),
     1   sptms2(i),sptchg2(i),sepsgm26(i),sepsgm212(i),sptwei2(i)
	end if
2	continue

	if (snb0.gt.0) then
	 read(urcon,100)TILDA
	 read(urcon,100)TILDA
	 do 3 i=1,snb0
	  read(urcon,109)sib21(i),sib22(i),skbond2(i),sreq2(i)
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
	   read(urcon,111)siangl21(i),siangl22(i),siangl23(i),
     1	   	skangl2(i),sangleq2(i)
111	   format(3(1x,i7),2(1x,f10.5))
	   sangleq2(i) = sangleq2(i)/pi180
4	  continue
	end if
        if(stotex2 .gt.0) then
           read(urcon,100)TILDA
           read(urcon,100)TILDA
           k = 0
           i = 0
           sexc21(0) = 0
7        continue
           read(urcon,117)j,k1
71       continue
           i = i + 1
           if (debug) write(stdo,*)' j k1 = ',j,k1
117        format(2(1x,i7))
             if (i.lt.j) then
              sexc21(i) = k
              go to 71
           else if (i.gt.j) then
              level = 1
              call alert(name,namel,'Fishy exclusion list!',21,level)
           else
              sexc21(j) = k + k1
              read(urcon,118)(sexc22(l),l=k+1,sexc21(j))
118         format(1x,10(i7,1x))
              k = k + k1
           end if
           if (k.ne.stotex2) go to 7
          end if

        if (stotspe2.gt.0) then
         read(urcon,100)TILDA
         read(urcon,100)TILDA
         do 8 j=1,stotspe2
           read(urcon,121)sspec21(j),sspec22(j),s2p14(1,j)
     1		,s2p14(2,j),s2p14(3,j)
121        format(2(1x,i7),1x,f30.5,2(1x,f15.5))
8	 continue
        end if
	return
	end
