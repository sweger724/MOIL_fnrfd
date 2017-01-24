         program structural_measures

	 include 'COMMON/LENGTH.BLOCK'
	 include 'COMMON/COORD.BLOCK'
	 include 'COMMON/CONNECT.BLOCK'
	 include 'COMMON/LINE.BLOCK'
	 include 'COMMON/DEBUG.BLOCK'
	 include 'COMMON/FREEZ.BLOCK'
	 include 'COMMON/CONVERT.BLOCK'
	 include 'COMMON/UNITS.BLOCK'
	 include 'COMMON/CCRD.BLOCK'
	   
c exmine structural measures from a dcd or a path file
c t(alpha,beta): tensor of inertia
c delta: measure of sphericity -- (3/2)*(SUM (lambda(i) - <lambda>)^2 ) / trace(t)^2
c S: shape measure -- 27* (PI (lambda(i) - <lambda>) / trace(t)^3
c
         integer ucon,urdcrd,uwflu,urcrd,urpcrd,uwtab,uwdelta,uws
         integer geti,of,nstep,i,j,ncoor,namel,n,k,l,i1,j1,error
         integer rbin
	 integer npick,level
         integer ipick(maxpt)
	 double precision tofi(3,3),delta,s,lambda(3)
	 character*13 name
	 logical find,fopen,ltab,ldelta,ls
	 logical lchr,ldyn,lpth,pickpt

	 stdi=5
	 stdo=6
         rbin = 1

	 npt=0
	 do 1 i=1,3
	  do 1 j=1,3
	   tofi(i,j)=0.d0
1	 continue
	 name='stru_measures'
	 namel=13
c  open junk file for rline
c
	 jnkf=25
	 open(unit=jnkf,status='scratch')
c    defalt parameters
	 pickpt=.false.
	 ltab=.false.
	 ldelta=.false.
	 ls=.false.
	 lchr=.false.
	 ldyn=.false.
	 lpth=.false.
	  ncoor=1
         lpstr = 1
2         continue
	  call rline(name,namel,stdi)
	  if (find ('norw')) norew=.true.
	  if (find ('file')) then
	    if (find ('conn')) then
	     ucon=of()
c  get connectivity
             call rconn(ucon)
            end if
	    if (find ('rcrd')) then
	     urcrd=of()
c   read reference structure
	     lchr = .true.
             call getcrd(urcrd,'CHARM')
	     do 80 i=1,npt
	     coor2(1,i)=coor(1,i)
	     coor2(2,i)=coor(2,i)
	     coor2(3,i)=coor(3,i)
80           continue
            end if
	    if (find ('rdyc')) then
		urdcrd=of()
		norew = .true.
		ldyn = .true.
	    end if
	    if (find ('rpth')) then
		 urpcrd=of()
		 norew = .true.
		 lpth = .true.
	    end if
	    if (find('wtab')) then 
	     uwtab=of()
	     ltab=.true.
            end if
	    if (find('wdel')) then
	      uwdelta=of()
	      ldelta=.true.
            end if
	    if (find('wris')) then
	      uws=of()
	      ls=.true.
            end if
          else
	    ncoor=geti('#crd',ncoor)
            if (find('pick')) then
	     call pick(ipick,npick)
	     pickpt=.true.
	    end if
	    if (find('action')) goto 3
          endif
          goto 2
3         continue
          if (.not. fopen(ucon)) then
	    level=1
	    call alert(name,namel,'ucon not opened',15,level)
          end if
c  initialize the vector nofreez and coordinates
	  inofrz=npt
	  do 4 i=1,npt
	   nofreez(i)=i
4         continue
c ============================================
        if (pickpt) then
	 inofrz = 0
         do 5 i=1,npt
          if(ipick(i) .ne. 0) then
		inofrz = inofrz + 1
		nofreez(inofrz)=i
          end if
5     continue
        end if
c ============================================
       if(lchr) then
	call str_shape(inofrz,nofreez,coor,tofi,lambda,error,
     1           delta,s,stdo,uwtab,uwdelta,uws)
       else if (ldyn) then
	do 6 k=1,ncoor
	 call rdyncrd(urdcrd,k,inofrz,nofreez,rbin)
	 call str_shape(inofrz,nofreez,coor,tofi,lambda,error,
     1           delta,s,stdo,uwtab,uwdelta,uws)
6	continue
       else if (lpth) then
	 do 7 k=1,ncoor
	  call rpath(urpcrd,k)
	  call str_shape(inofrz,nofreez,coor,tofi,lambda,error,
     1           delta,s,stdo,uwtab,uwdelta,uws)
7	 continue
	end if

          stop
	  end
