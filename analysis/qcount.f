           program qcount
c
c count the occurance of q along free energy path to do umbrella
c sampling calculations
c
	   integer uqfi,urho
	   double precision getd

	   integer steps
	   double precision xfc,q,q0,qmin,qmax
	   real rho(100)
	   character*6 name
	   logical find
	   data rho/100*0/

           stdi=5
	   stdo=6
	    open(unit=jnkf,status='scratch')
	    nstru= 45000
	    steps = 10
	    xfc  = 0.5d0
	do 100 istep=1,steps
		rewind uqfi
		do 4 i=1,istep-1
		 do 4 j=1,nstru+1
		 read(uqfi,*) q0
4		continue
		read (uqfi,*)q0
		qmin=1.d6
		qmax=-1.d6
		do 5 i=1,nstru
		 read(uqfi,*)q
		 if (q.lt.qmin) qmin = q
		 if (q.gt.qmax) qmax = q
5		continue
		write(*,*)' q0 qmin qmax ',q0,qmin,qmax
		dq = (qmax-qmin)/99.
		rewind uqfi
		do 6 i=1,istep-1
		 do 6 j=1,nstru+1
		 read(uqfi,*) q0
6		continue
		 read(uqfi,*)q0
		do 7 i=1,nstru
		 read(uqfi,*)q
		 j = (q-qmin)/dq+1
		 rho(j) = rho(j)+1
7		continue
		do 8 j=0,99
		 q = qmin + j*dq
		 u = 0.6*log(rho(j)) + xfc*(q-q0)**2
		 write(urho,*)q,u
8		continue
100	continue
	stop
	end
