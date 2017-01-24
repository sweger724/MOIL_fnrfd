        subroutine read_CG_parameters(uCGpar)

        implicit none
     

c
c local
c
        character*18 name
        character*5 getchar
        character*2 keyword
        character*4 mono1
        character*4 mono2
        integer namel,coun,lc,level,id1,id2,getCGid,i,j,IDi,IDj
        logical find
        integer uCGpar 

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision getd, Tmp(Namino,12), Tmp2(Namino*Namino,10)

        name = 'read_CG_parameters'
        namel = 18
        lc    = 4

        echoNO=.TRUE.

        call rline(name,namel,uCGpar)
        
        if (.not.find('LJP1')) then
          level = 1
          call alert(name,namel,'Missing LJP1 header',19,level)
        end if

c loop for reading LJP1 parameters
c
        do coun=1,Namino * (Namino+1)/2
          call rline(name,namel,uCGpar)

          if (find('DONE')) then
             level = 1
             call alert(name,namel,'Missing LJP1 parameters',20,level)
          end if
          
          mono1  = getchar('MON1','????',lc)
          id1 = getCGid(mono1)
          mono2  = getchar('MON2','????',lc)
          id2 = getCGid(mono2)
          IDi = Namino*(id1-1)+id2
          IDj = Namino*(id2-1)+id1
C         write(6,*)coun,id1,id2,IDi,IDj
          do i =0,9
            write(keyword,'(a1,i1.1)') 'a',i
C           write(stdo,*) "Read: ",keyword
            LJa(IDi,i) = getd(keyword,999.d0)
            LJa(IDj,i) = LJa(IDi,i)
C           write(stdo,*) "Read: ",LJa(IDi,i)
          end do
   
          LJr(IDi,1) = getd('r1',999.d0)
          LJr(IDi,2) = getd('r2',999.d0)
          LJr(IDi,3) = getd('r3',0.d0)
C         write(stdo,*) "Read: r123",LJr(IDi,1),LJr(IDi,2),LJr(IDi,3)
          LJr(IDj,1) = LJr(IDi,1)
          LJr(IDj,2) = LJr(IDi,2)
          LJr(IDj,3) = LJr(IDi,3)

          LJa(IDi,10) = getd('cutt',4.2d0)
          LJa(IDi,11) = getd('midl',8.55d0)
          LJa(IDj,10) = LJa(IDi,10)
          LJa(IDj,11) = LJa(IDi,11)
          
        end do

        call rline(name,namel,uCGpar)
        
        if (.not.find('DONE')) then
          level = 1
          call alert(name,namel,
     &    'Missing DONE after LJP1 parameters',23,level)
        end if
        write(stdo,*)' *** DONE WITH LJP1 POTENTIAL PARAMETERSS *** '
        
c DONE with reading LJP1 parameters
        
        call rline(name,namel,uCGpar)

        if (.not.find('LJP2')) then
          level = 1
          call alert(name,namel,'Missing LJP2 header',19,level)
        end if

c loop for reading LJP2 parameters
c
        do coun=1,Namino * Namino
          call rline(name,namel,uCGpar)

          if (find('DONE')) then
             level = 1
             call alert(name,namel,'Missing LJP2 parameters',20,level)
          end if

          mono1  = getchar('MON1','????',lc)
          id1 = getCGid(mono1)
          mono2  = getchar('MON2','????',lc)
          id2 = getCGid(mono2)
          IDi = Namino*(id1-1)+id2 + Namino**2
          do i =0,9
            write(keyword,'(a1,i1.1)') 'a',i
            LJa(IDi,i) = getd(keyword,999.d0)
          end do

          LJr(IDi,1) = getd('r1',999.d0)
          LJr(IDi,2) = getd('r2',999.d0)
          LJr(IDi,3) = getd('r3',0.d0)
C          write(stdo,*) "Read: r123",coun,LJr(IDi,1),
C     &                         LJr(IDi,2),LJr(IDi,3)
          LJa(IDi,10) = getd('cutt',4.2d0)
          LJa(IDi,11) = getd('midl',8.55d0)

        end do

        call rline(name,namel,uCGpar)

        if (.not.find('DONE')) then
          level = 1
          call alert(name,namel,
     &    'Missing DONE after LJP2 parameters',23,level)
        end if
        write(stdo,*)' *** DONE WITH LJP2 POTENTIAL PARAMETERSS ***'

c DONE with reading LJP2 parameters

        call rline(name,namel,uCGpar)

        if (.not.find('LJP3')) then
          level = 1
          call alert(name,namel,'Missing LJP3 header',19,level)
        end if

c loop for reading LJP3 parameters
c
        do coun=1,Namino * (Namino+1)/2
          call rline(name,namel,uCGpar)

          if (find('DONE')) then
             level = 1
             call alert(name,namel,'Missing LJP3 parameters',20,level)
          end if

          mono1  = getchar('MON1','????',lc)
          id1 = getCGid(mono1)
          mono2  = getchar('MON2','????',lc)
          id2 = getCGid(mono2)
          IDi = Namino*(id1-1)+id2 + 2*Namino**2
          IDj = Namino*(id2-1)+id1 + 2*Namino**2
C         write(6,*)coun,id1,id2,IDi,IDj
          do i =0,9
            write(keyword,'(a1,i1.1)') 'a',i
C           write(stdo,*) "Read: ",keyword
            LJa(IDi,i) = getd(keyword,999.d0)
            LJa(IDj,i) = LJa(IDi,i)
C           write(stdo,*) "Read: ",LJa(IDi,i)
          end do

          LJr(IDi,1) = getd('r1',999.d0)
          LJr(IDi,2) = getd('r2',999.d0)
          LJr(IDi,3) = getd('r3',0.d0)
C         write(stdo,*) "Read: r123",LJr(IDi,1),LJr(IDi,2),LJr(IDi,3)
          LJr(IDj,1) = LJr(IDi,1)
          LJr(IDj,2) = LJr(IDi,2)
          LJr(IDj,3) = LJr(IDi,3)

          LJa(IDi,10) = getd('cutt',4.2d0)
          LJa(IDi,11) = getd('midl',8.55d0)
          LJa(IDj,10) = LJa(IDi,10)
          LJa(IDj,11) = LJa(IDi,11)

        end do

        call rline(name,namel,uCGpar)

        if (.not.find('DONE')) then
          level = 1
          call alert(name,namel,
     &    'Missing DONE after LJP3 parameters',23,level)
        end if
        write(stdo,*)' *** DONE WITH LJP3 POTENTIAL PARAMETERSS *** '

c DONE with reading LJP3 parameters

c Reading multiple- well  energy parameters
c
C  set CA - CA bond parameters
        CGBond(0,1) = 1
        CGBond(0,2) = 3.804d0
        CGBond(0,3) = 0.d0
        CGBond(0,4) = 200.d0

C set cis CA - CA bond parameters
        CGBond(21,1) = 1
        CGBond(21,2) = 2.965
        CGBond(21,3) = 0.d0
        CGBond(21,4) = 100.d0

C set CYS - bridge parameters
        CGBond(22,1) = 1
        CGBond(22,2) = 2.365
        CGBond(22,3) = 0.d0
        CGBond(22,4) = 300.d0

        call ReadWellParameters(name,namel,uCGpar,Namino,Tmp,"BOND")
        do i = 1, Namino
          do j = 1,12
            CGBond(i,j)=Tmp(i,j)
          end do
        end do

        call ReadWellParameters(name,namel,uCGpar,Namino,Tmp,"ANG1")
        do i = 1, Namino
          do j = 1,12
             CGAngle(i,j)=Tmp(i,j)
          end do
        end do

        call ReadWellParameters(name,namel,uCGpar,Namino,Tmp,"ANG2")
        do i = 1, Namino
          do j = 1,12
             CGAngle(Namino+i,j)=Tmp(i,j)
          end do
        end do 

        call ReadWellParameters(name,namel,uCGpar,Namino,Tmp,"ANG3")
        do i = 1, Namino
          do j = 1,12
             CGAngle(2*Namino+i,j)=Tmp(i,j)
          end do
        end do

C Reading parameters for torsional angles
C
      call ReadTorsionParameters(name,namel,uCGpar,Namino,Tmp2,"TOR1")
        do i = 1, Namino*Namino
          do j = 1,10
             CGTor(i,j)=Tmp2(i,j)
          end do
        end do

      call ReadTorsionParameters(name,namel,uCGpar,Namino,Tmp2,"TOR2")
        do i = 1, Namino*Namino
          do j = 1,10
             CGTor(Namino**2+i,j)=Tmp2(i,j)
          end do
        end do

      call ReadTorsionParameters(name,namel,uCGpar,Namino,Tmp2,"TOR3")
        do i = 1, Namino*Namino
          do j = 1,10
             CGTor(2*Namino**2+i,j)=Tmp2(i,j)
          end do
        end do

      call ReadTorsionParameters(name,namel,uCGpar,Namino,Tmp2,"TOR4")
        do i = 1, Namino*Namino
          do j = 1,10
             CGTor(3*Namino**2+i,j)=Tmp2(i,j)
          end do
        end do

c check if end of file
        call rline(name,namel,uCGpar)
        if (.not. find('*EOD')) then
          level = 0
          call alert(name,namel,'Missing the end of file',23,level)
        end if

        echoNO=.FALSE.

        return
        end

C*************************************************************************
        subroutine CheckEndOfFile(name,namel,uCGpar,message,length)
           implicit none
           character*(*) message
           character*(*) name
           integer namel, uCGpar, level, length
           logical find

           call rline(name,namel,uCGpar)
           if (find('*EOD')) then
             level = 1
             call alert(name,namel,message,length,level)
            end if
        end
        
C*************************************************************************
        subroutine ReadWellParameters(name,namel,uCGpar,N,
     &             DataStore,KEYWORD)
        implicit none
        
        include 'COMMON/UNITS.BLOCK'

        integer N, namel, uCGpar
        character*4 KEYWORD
        character*(*) name
        double precision DataStore(1:N,12)
        
        character*30 message
        character*4 mono 
        character*5 getchar
        integer id, type1, coun, lc, getCGid, level
        integer geti
        double precision getd
        logical find

        lc = 4
        
c Check if end of file
        write(message,'(a,a4,a)')'No ',KEYWORD,' and bellow'
        call CheckEndOfFile(name,namel,uCGpar,message,20)

        if (.not.find(KEYWORD)) then
                level = 1
                write(message,'(a4,a)') KEYWORD,' header missing'
                call alert(name,namel,message,19,level)
        end if

        do coun=1,N
          call rline(name,namel,uCGpar)

          if (find('DONE')) then
             level = 1
      write(message,'(a,a4,a)')'Missing some ',KEYWORD,' parameters'
             call alert(name,namel,message,22,level)
          end if

          mono  = getchar('MONO','????',lc)
          id = getCGid(mono)
          type1 = geti('type',-1)
          if (id .eq. -1) write(6,*) "ERROR!"
          DataStore(id,1) = type1
          if (type1 .eq. 1) then
             DataStore(id,2) = getd('r',999.d0)
             DataStore(id,3) = getd('e',999.d0)
             DataSTore(id,4) = getd('k',999.d0)
          endif

          if (type1 .eq. 2) then
             DataStore(id,2) = getd('r1',999.d0)
             DataStore(id,3) = getd('e1',999.d0)
             DataStore(id,4) = getd('k1',999.d0)
             DataStore(id,5) = getd('r2',999.d0)
             DataStore(id,6) = getd('e2',999.d0)
             DataStore(id,7) = getd('k2',999.d0)
             DataStore(id,8) = getd('beta',999.d0)
C            write(6,*)"Read: ",DataStore(id,2), DataStore(id,3),
C     & DataStore(id,4), DataStore(id,5), DataStore(id,6)            
          endif

          if (type1 .eq. 3) then
             DataStore(id,2) = getd('r1',999.d0)
             DataStore(id,3) = getd('e1',999.d0)
             DataStore(id,4) = getd('k1',999.d0)
             DataStore(id,5) = getd('r2',999.d0)
             DataStore(id,6) = getd('e2',999.d0)
             DataStore(id,7) = getd('k2',999.d0)
             DataStore(id,8) = getd('bet1',999.d0)
             DataStore(id,9) = getd('r3',999.d0)
             DataStore(id,10) = getd('e3',999.d0)
             DataStore(id,11) = getd('k3',999.d0)
             DataStore(id,12) = getd('bet2',999.d0)
          endif 

        end do

        call rline(name,namel,uCGpar)
        if (find('DONE')) then
          write(stdo,*)' *** DONE WITH ',KEYWORD,' PROPERTIES *** '
        else
          level = 1
          write(message,'(a,a4,a)')'NO DONE AFTER '
     &          ,KEYWORD,' parameters'
          call alert(name,namel,message,30,level)
        end if

        return
        end


C*************************************************************************
        subroutine ReadTorsionParameters(name,namel,uCGpar,N,
     &             DataStore,KEYWORD)
        implicit none
        
        include 'COMMON/UNITS.BLOCK'
        integer N, namel, uCGpar
        character*4 KEYWORD
        character*(*) name
        character*5 getchar
        double precision DataStore(N*N,*)

        character*30 message
        character*4 mono1
        character*4 mono2
        integer id1,id2,type1,type2,coun, getCGid, level, geti,lc,IDi
        double precision getd
        logical find

        lc=4

c Check if end of file
        write(message,'(a,a4,a)')'No ',KEYWORD,' and bellow'
        call CheckEndOfFile(name,namel,uCGpar,message,20)

        if (.not.find(KEYWORD)) then
                level = 1
                write(message,'(a4,a)') KEYWORD,' header missing'
                call alert(name,namel,message,19,level)
        end if

        do coun=1,N*N
          call rline(name,namel,uCGpar)

          if (find('DONE')) then
             level = 1
             write(message,'(a,a4,a)')'Missing some '
     &             ,KEYWORD,' parameters'
             call alert(name,namel,message,22,level)
          end if

          mono1  = getchar('MON1','????',lc)
          mono2  = getchar('MON2','????',lc)
          id1 = getCGid(mono1)
          id2 = getCGid(mono2)
          IDi = N*(id1-1)+id2
             DataStore(IDi,1) = getd('c1',999.d0)
             DataStore(IDi,2) = getd('c2',999.d0)
             DataSTore(IDi,3) = getd('c3',999.d0)
             DataSTore(IDi,4) = getd('c4',0.d0)
             DataSTore(IDi,5) = getd('c5',0.d0)
             DataStore(IDi,6) = getd('s1',999.d0)
             DataStore(IDi,7) = getd('s2',999.d0)
             DataStore(IDi,8) = getd('s3',999.d0)
             DataStore(IDi,9) = getd('s4',0.d0)
             DataStore(IDi,10) =getd('s5',0.d0)
        end do

        call rline(name,namel,uCGpar)
        if (find('DONE')) then
          write(stdo,*)' *** DONE WITH ',KEYWORD,' PROPERTIES *** '
        else
          level = 1
          write(message,'(a,a4,a)')'NO DONE AFTER '
     &          ,KEYWORD,' parameters'
          call alert(name,namel,message,30,level)
        end if

        return
        end     
