	Subroutine cdie3()

C	cdie 

      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'
	include 'COMMON/GETVEC.BLOCK'
c	include 'COMMON/SCNDRV.BLOCK'

	double precision pick 
	double precision epstmp
	double precision rx,ry,rz,r2,s2,a,b,e1,e2,q,df,df1,df2,ai,bi
	double precision s,s6,tmp
	double precision qi,b6,b12
	double precision xi,yi,zi,dxi,dyi,dzi
c add for second derivatives
	double precision s8,s3,abq1,abq2,add

	integer i,j,k,jbeg1,jend1,jbeg2,jend2,jbeg3,jend3
	integer ptbeg,ptend
c adds for second derivatives
	integer index,jndex,kndex

	if (eelyes) then
		epstmp = 332.0716d0/eps
	else
		epstmp =  0.0d0
	end if
	if (nocut) then
		cutvdw2 = 10000.d0
		cutele2 = 10000.d0
	end if

	if (evdyes) then
		pick = 1.0d0
	else
		pick = 0.0d0
	end if
	e_vdw = 0.d0
	e_el  = 0.d0


c -------------------------------------------------------------
c       yael

	if (prll_on_off) then	
	   ptbeg = dpoipt(monp(my_pe))+1
	   if (my_pe.eq.(num_pes-1)) then
	      ptend = npt-1
	   else
	      ptend = dpoipt(monp(my_pe+1))
	   endif   
	   if (my_pe.eq.0) then
	      jbeg1 = 1
	      jbeg2 = 1
	      jbeg3 = 1
           else
	      jbeg1=point1(ptbeg-1)+1
	      jbeg2=point2(ptbeg-1)+1
	      jbeg3=point3(ptbeg-1)+1
	   endif
	else
	   ptbeg = 1
	   ptend = npt-1
	   jbeg1=1
	   jbeg2=1
	   jbeg3=1
	end if	
c
c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
c


	do 400 i=ptbeg,ptend


		jend1 = point1(i)
		tmp = 1.d0/ptwei(i)
		ai  = epsgm12(i)*pick
                if (arith) then
		bi  = epsgm6(i)
                else
                bi  = epsgm6(i)*pick
                endif
		qi  = ptchg(i)*epstmp
		xi  = coor(1,i)
		yi  = coor(2,i)
		zi  = coor(3,i)
		dxi = 0.d0
		dyi = 0.d0
		dzi = 0.d0

		if (jbeg1.le.jend1) then

		do 200 k=jbeg1,jend1
			j=list1(k)
			rx = xi - coor(1,j)
			ry = yi - coor(2,j)
			rz = zi - coor(3,j)
			r2=rx*rx+ry*ry+rz*rz
			s2=1.0d0/r2

			q = qi*ptchg(j)
			s = dsqrt(s2)

			if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
				df1 = 0.d0
				if ((lesid(i).ne.0) .and.
     *			       (lesid(i) .eq.lesid(j))) then
				 q = q*tmp
				end if
				e2  = q*s
				df2 = -e2*s2
				s3  = s2*s
				abq1 = 3.d0*q*s3*s2
				abq2 = -q*s3
				go to 100
			end if

			a= ai*epsgm12(j)
                        if (arith) then
			b= bi+epsgm6(j)
                        else
                        b= bi*epsgm6(j)
                        endif
			if ((lesid(i).ne.0) .and.
     *			    (lesid(i) .eq.lesid(j)))  then
			   a = a*tmp
			   b = b*tmp
			   q = q*tmp
			end if

			e2 = q*s
			df2 = -e2*s2
			s3  = s2*s
			s6  = s3*s3
			s8  = s2*s6
                        if (arith) then
                        b6 = b*b*b*b*b*b
                        b12 = b6*b6
			abq1 = 3.d0*q*s3*s2 + 168.d0*a*b12*s8*s8
     1                         -48.d0*b6*a*s8*s2
			abq2 =-q*s3-12.d0*a*b12*s6*s8+6.d0*b6*a*s8
                        b   = s2*b*b
                        s6  = b*b*b
			e1 = a*(s6*s6 - s6) 
			df1 = -6.0d0*s2*a*(2.d0*s6*s6-s6)
                        else
               abq1 = 3.d0*q*s3*s2 + 168.d0*a*s8*s8-48.d0*b*s8*s2
                        abq2 = -q*s3 - 12.d0*a*s6*s8 + 6.d0*b*s8
                        a = a*s6*s6
                        b = b*s6
                        e1 = a - b
                        df1 = -6.0d0*s2*(a+e1)
                        endif
			e_vdw = e_vdw + e1 

c the electrostatic part should be done without buffering
100			continue
			df = df1 + df2
c
c start second derivative part
c
                        index = 6*(i-1)+1
                        jndex = 6*(j-1)+1
                        kndex = 6*(k-1)+1
                        add         = abq1*rx*rx + abq2
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
                        add         = abq1*ry*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
                        add         = abq1*rz*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
                        add        = abq1*ry*ry + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
                        add        = abq1*rz*ry
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
                        add        = abq1*rz*rz + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_1(kndex) = -add
c
c end second derivative part
c
			rx = df*rx
			ry = df*ry
			rz = df*rz
			dxi = dxi + rx
			dyi = dyi + ry
			dzi = dzi + rz
			dpot(1,j) = dpot(1,j) - rx
			dpot(2,j) = dpot(2,j) - ry
			dpot(3,j) = dpot(3,j) - rz
			e_el = e_el + e2


200		continue
		end if
		jbeg1 = jend1 + 1
  
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c


              jend2 = point2(i)


              if (jbeg2.le.jend2) then
              do 205 k=jbeg2,jend2
                      j=list2(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutvdw2) go to 205
                        a = ai*epsgm12(j)
                        if (arith) then
                        b = bi+epsgm6(j)
                        else
                        b = bi*epsgm6(j)
                        endif
			if ((lesid(i).ne.0) .and.
     *			    (lesid(i) .eq.lesid(j)))  then
                         a = a*tmp
                         b = b*tmp
                      end if
                      s2=1.0d0/r2
                      s6=s2*s2*s2
                        s8=s2*s6
                        if (arith) then
                        b6 = b*b*b*b*b*b
                        b12 = b6*b6
                        abq1=168.d0*a*b12*s8*s8-48.d0*b6*a*s8*s2
                        abq2=-12.d0*a*b12*s6*s8+6.d0*b6*a*s8
                        else
                        abq1=168.d0*a*s8*s8-48.d0*b*s8*s2
                        abq2=-12.d0*a*s6*s8+6.d0*b*s8
                        endif
                        index = 6*(i-1)+1
                        jndex = 6*(j-1)+1
                        kndex = 6*(k-1)+1
                        add         = abq1*rx*rx + abq2
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add
                        add         = abq1*ry*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add
                        add         = abq1*rz*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add
                        add         = abq1*ry*ry + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add
                        add         = abq1*rz*ry
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add
                        add         = abq1*rz*rz + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_2(kndex) = -add

                      if (arith) then
                      b   = s2*b*b
                      s6  = b*b*b
                      e1  = a*(s6*s6 - s6)
                      df = -6.0d0*s2*a*(2.d0*s6*s6-s6)
                      else
                      a = a*s6*s6
                      b = b*s6
                      e1 = a - b
                      df = -6.0d0*s2*(a+e1)
                      endif

                      rx = df*rx
                      ry = df*ry
                      rz = df*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz

                      e_vdw = e_vdw + e1 
205           continue
              end if
              jbeg2 = jend2 + 1



c start THIRD loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c


              jend3 = point3(i)


              if (jbeg3.le.jend3) then

              do 210 k=jbeg3,jend3
                      j=list3(k)
                      rx = xi - coor(1,j)
                      ry = yi - coor(2,j)
                      rz = zi - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutele2) go to 210

                      q=qi*ptchg(j)

c chen
			if ((lesid(i).ne.0) .and.
     *			    (lesid(i) .eq.lesid(j)))  q = q * tmp
                      s2=1.0d0/r2
                      s = dsqrt(s2)
                        s3 = s2*s
                        abq1=3.d0*q*s3*s2
                        abq2=-q*s3
c                       write(*,*)'#3 abq1 abq2 ',abq1,abq2
                        index = 6*(i-1)+1
                        jndex = 6*(j-1)+1
                        kndex = 6*(k-1)+1
                        add         = abq1*rx*rx + abq2
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add
                        add         = abq1*ry*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add
                        add         = abq1*rz*rx
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add
                        add         = abq1*ry*ry + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add
                        add         = abq1*rz*ry
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add
                        add         = abq1*rz*rz + abq2
                        index = index + 1
                        jndex = jndex + 1
                        kndex = kndex + 1
                        diag(index) = diag(index) + add
                        diag(jndex) = diag(jndex) + add
c                        d2nb_3(kndex) = -add

                      e2 = q*s
                      df2 = -e2*s2

                      rx = df2*rx
                      ry = df2*ry
                      rz = df2*rz
                      dxi = dxi + rx
                      dyi = dyi + ry
                      dzi = dzi + rz
                      dpot(1,j) = dpot(1,j) - rx
                      dpot(2,j) = dpot(2,j) - ry
                      dpot(3,j) = dpot(3,j) - rz
                      e_el = e_el + e2


210           continue
              end if
              jbeg3 = jend3 + 1



		dpot(1,i) = dpot(1,i) + dxi
		dpot(2,i) = dpot(2,i) + dyi
		dpot(3,i) = dpot(3,i) + dzi
400	continue

c remove electrostatic contributions of the type: VP - its related real prts
        if (vp_flag) call vp_cdie()


	return 
	end
