	subroutine add_les()
c Subroutine for multiplying les particles.

	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/PROPERT.BLOCK'
	include 'COMMON/MONOMERS.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/TMP_CONNECT.BLOCK'

	character*7 name
	integer namel,i,j,k,l,jbeg,kbeg,level,lid,cpyn,m,n,temp
	integer ipick(maxpt),n_cpy(maxpt),poicpy(maxpt)
	double precision weight
	logical find,noles
	integer geti
	integer igroup

	integer i_save

	name  = 'add_les'
	namel = 7
        igroup = 1
	noles = .true.
	if (lesflag) then
	     	level = 1
	  	call alert(name,namel,
     *			   'Old fashion les exists',20,level)
	end if
	
	do 101 i = 1,npt
		n_cpy(i) = 1
		lesid(i) = 0
101 	continue
	lestyp = 0

			
	lid = 100
1	continue
		call rline(name,namel,stdi)
		if (find('*EOD')) then
			if (noles) return
			 go to 2
		end if
		noles = .false.
		if (find('MULT')) then
		   lestyp = lestyp + 1
		   lesflag = .true.
	   	   call pick(ipick,igroup)
	   	   cpyn = 0
	   	   cpyn = geti('#cpy',cpyn)
		   if (cpyn .eq. 0) then
			    level = 1
	  		    write(stdo,*) moname(poimon(i)),' ',
     *			        poimon(i),' ', ptnm(i)
     	  		    call alert(name,namel,
     *				'Number of copies in les',
     *			        23,level)
	       	    end if
	      	    if (cpyn .gt. maxcpy) then
	      		level = 1
	  		write(stdo,*) moname(poimon(i)),' ',
     *				     poimon(i),' ', ptnm(i)
	  		call alert(name,namel,
     *			   'Particle with too many copies',29,level)
	       	     end if
		     temp = geti('lesi',-1)
		     if (temp .gt. 100) then
	      		  	level = 1
	  		  	call alert(name,namel,
     *			  	'Explicit les id more then 100',
     *					28,level)				
		        else
		   	    	if (temp .eq. -1) then
			     		l = lid
		   	     		lid = lid + 1
			            else 
					l = temp
				end if

		      end if
				
	   	      do 100 i = 1,npt
			if (ipick(i) .eq. 1) then
				lesid(i) = l
		 		n_cpy(i) = cpyn
			end if
100	   	       continue
		end if
	go to 1

2 	continue

	do 290 i=1,totmon
		do  291 j=1,mdivlist(0,i)
			i_save = mdivlist_tmp(j,i)
			do 292 k=mdivlist(j,i)+1,mdivlist(j+1,i)
			i_save = i_save + n_cpy(k)
292		continue
		mdivlist_tmp(j+1,i) = i_save
291		continue
		mdivlist_tmp(1,i+1) = mdivlist_tmp(mdivlist(0,i)+1,i)
290	continue

	do 293 i=1,totmon
		do 293 j=1,mdivlist(0,i)+1
			mdivlist(j,i) = mdivlist_tmp(j,i)
293	continue

c Create temporary particle list.
	tnpt = 1
	do 300 i = 1,npt
		weight = 1.d0/n_cpy(i)
		do 310 j = 1,n_cpy(i)
			tpoimon(tnpt) 	= poimon(i)
			if (mdivyes) then
				tpoidmon(tnpt) = poidmon(i)
			end if
			tptid(tnpt)   	= ptid(i)
			tptnm(tnpt) 	= ptnm(i)			    
			tptms(tnpt) 	= ptms(i) * weight		
			tptchg(tnpt) 	= ptchg(i) * weight		    
			if (.not.arith) then
			 tepsgm6(tnpt) 	= epsgm6(i) * weight			
			else
			 tepsgm6(tnpt)  = epsgm6(i)
			end if
			tepsgm12(tnpt) 	= epsgm12(i) * weight		    
			tptwei(tnpt) 	= weight			    
			tlesid(tnpt) 	= lesid(i)			    
			tflagchr(tnpt) 	= flagchr(i)
			if (lesid(i) .eq. 0) then
				cplbl(tnpt) = 0
			    else 
				cplbl(tnpt) = j
			end if
			tnpt = tnpt + 1
310		continue
		poicpy(i) = tnpt - 1
300	continue
	tnpt = tnpt - 1
c Creat tempotrary bond list 
	tnb = 1
	do 400 i = 1,nb
		m = ib1(i)
		n = ib2(i)
		if (m .eq. 1) then 
			jbeg = 1
		   else 
			jbeg = poicpy(m - 1) + 1
		end if 
		if (n .eq. 1) then 
			kbeg = 1
		   else 
			kbeg = poicpy(n - 1) + 1
		end if 
		if (((lesid(m) .eq. 0) .and. (lesid(n) .eq. 0)) .or.
     *	              (lesid(m) .ne. lesid(n))) then
			weight = 1.d0/n_cpy(n) * 1.d0/n_cpy(m)
			do 410 j = jbeg, poicpy(m)
			   do 420 k = kbeg, poicpy(n)
				tib1(tnb)  = j
				tib2(tnb)  = k
				tkbond(tnb) = kbond(i) * weight		    
				treq(tnb)   = req(i)
				tnb = tnb + 1
420			   continue
410			continue
		   else
			k = kbeg
			weight = 1.d0/n_cpy(n)   
			do 430 j = jbeg, poicpy(m)
				tib1(tnb)  = j
				tib2(tnb)  = k
				tkbond(tnb) = kbond(i) * weight		    
				treq(tnb)   = req(i)
				tnb = tnb + 1
				k = k + 1
430			continue
		end if
400 	continue
	tnb = tnb - 1
c Creat tempotrary morse bond list 
	tnmb = 1
	do 450 i = 1,nmb
		m = imb1(i)
		n = imb2(i)
		if (m .eq. 1) then 
			jbeg = 1
		   else 
			jbeg = poicpy(m - 1) + 1
		end if 
		if (n .eq. 1) then 
			kbeg = 1
		   else 
			kbeg = poicpy(n - 1) + 1
		end if 
		if (((lesid(m) .eq. 0) .and. (lesid(n) .eq. 0)) .or.
     *	             (lesid(m) .ne. lesid(n))) then
			weight = 1.d0/n_cpy(n) * 1.d0/n_cpy(m)
			do 460 j = jbeg, poicpy(m)
			   do 470 k = kbeg, poicpy(n)
				timb1(tnb)  = j
				timb2(tnb)  = k
				trmeq(tnb)  = rmeq(i)
				tD(tnb)     = D(i) * weight
				talpha(tnb) = alpha(i)
				tbeta1(tnb) = beta1(i)
				tnmb = tnmb + 1
470			   continue
460			continue
		   else
			k = kbeg
			weight = 1.d0/n_cpy(n)   
			do 480 j = jbeg, poicpy(m)
				timb1(tnb)  = j
				timb2(tnb)  = k
				trmeq(tnb)  = rmeq(i)
				tD(tnb)     = D(i) * weight
				talpha(tnb) = alpha(i)
				tbeta1(tnb) = beta1(i)
				tnmb = tnmb + 1
				k = k + 1
480			continue
		end if
450 	continue
	tnmb = tnmb - 1
c Update poipt

	do 500 i = 1, totmon
		poipt(i) = poicpy(poipt(i))
500 	continue

c	if (mdivyes) then
c		i_save = tpoidmon(1)
c		do 510 i=1,tnpt
c			j = 0
c			if (tpoidmon(i).ne.i_save) then
c				j              = j + 1
c				dpoipt(j) = i - 1
c				i_save         = tpoidmon(i)
c			end if
c510		continue
c		i_save = -1
c		do 511 i=1,totmon
c		 do 511 j=1,mdivlist(0,i)
c			i_save = i_save + 1
c			mdivlist(j,i) = dpoipt(i_save)
c511		continue
c	end if

c Copy final molecule

	npt = tnpt
	nb  = tnb
	do 600 i = 1,npt
		poimon(i)   = tpoimon(i)
		ptid(i)     = tptid(i)			    
		ptnm(i)     = tptnm(i)			    
		ptms(i)     = tptms(i)			    
		ptchg(i)    = tptchg(i)			    
		epsgm6(i)   = tepsgm6(i)			    
		epsgm12(i)  = tepsgm12(i)			    
		ptwei(i)    = tptwei(i)			    
		lesid(i)    = tlesid(i)			    
		flagchr(i)  = tflagchr(i)
600 	continue

	do 700 i = 1,nb
		ib1(i)   = tib1(i) 		    
		ib2(i)   = tib2(i) 		    
		kbond(i) = tkbond(i) 		    
		req(i)   = treq(i)
700	continue
	do 800 i = 1,nmb
		imb1(i)  = timb1(i) 		    
		imb2(i)  = timb2(i) 		    
		rmeq(i)  = trmeq(i)
		D(i)     = tD(i)
		alpha(i) = talpha(i)
		beta1(i) = tbeta1(i)
800	continue
	pbulk(1) = poipt(totmon)
	return
	end

