	Subroutine cdie_spcl(n)

C	cdie 

      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
	include 'COMMON/CONVERT.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/CONSPECL1.BLOCK'
      	include 'COMMON/CONSPECL2.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SPECL.BLOCK'

	double precision pick
	double precision epstmp
	double precision rx,ry,rz,r2,s2,a,b,e1,e2,q,df,df1,df2,ai,bi
	double precision s,s6,tmp
	double precision qi
	double precision xi,yi,zi,dxi,dyi,dzi
        double precision dfx(maxmorsb),dfy(maxmorsb),dfz(maxmorsb)
        double precision tempdx1(maxspcl),tempdy1(maxspcl)
        double precision tempdz1(maxspcl)
        double precision xtemp1(maxspcl),ytemp1(maxspcl),ztemp1(maxspcl)
        double precision tempdx2(maxspcl),tempdy2(maxspcl)
        double precision tempdz2(maxspcl)
        double precision xtemp2(maxspcl),ytemp2(maxspcl),ztemp2(maxspcl)
	integer i,j,k,n,initial,factr(maxspcl)
	integer j1,k1
	integer jbeg11(0:maxmorsb)
	integer jbeg21(0:maxmorsb)
        integer jend1


	if (eelyes) then
		epstmp = kofdie/eps
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

	se_vdw1 = 0.0d0
	se_vdw2 = 0.0d0
	se_el1=0.d0
	se_el2=0.d0
	poitype1(0) = 0
	poitype2(0) = 0

        spoint11(0) = 0
        initial = poitype1(n-1)+1
	jbeg11(n)=spoint11(poitype1(n-1))+1
c
c  For "special" all nonbonded interactions of one electronic
c curve are in a single list.
c
	do 400 i=initial,poitype1(n)-1

                   if (newpoit(i) .eq. imb1(n)) then
                    factr(i)=1
                    else if (newpoit(i) .eq. imb2(n)) then
                     factr(i)=-1
                    else
                     factr(i)=0
                   end if
		jend1 = spoint11(i)
		tmp = 1.d0/sptwei1(i)

c       ileana

		if(.not.arith)then
		ai  = sepsgm112(i)*pick
		bi  = sepsgm16(i)*pick
		else
		ai  = sepsgm112(i)*pick
		bi  = sepsgm16(i)

		endif

		qi  = sptchg1(i)*epstmp
		xi  = coor(1,newpoit(i))
		yi  = coor(2,newpoit(i))
		zi  = coor(3,newpoit(i))
		dxi = 0.d0
		dyi = 0.d0
		dzi = 0.d0

		if (jbeg11(n).le.jend1) then

		do 200 k=jbeg11(n),jend1
			j=slist11(k)
			rx = xi - coor(1,newpoit(j))
			ry = yi - coor(2,newpoit(j))
			rz = zi - coor(3,newpoit(j))
			r2=rx*rx+ry*ry+rz*rz
			s2=1.0d0/r2

			q = qi*sptchg1(j)
			s = dsqrt(s2)
			e2 = q*s
			df2 = -e2*s2
			if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
				df1 = 0.d0
				if ((slesid1(i).ne.0) .and.
     *			       (slesid1(i) .eq.slesid1(j))) q = q*tmp
				go to 100
			end if

c       ileana

			if(.not.arith)then

			a= ai*sepsgm112(j)
			b= bi*sepsgm16(j)

			else
			a= ai*sepsgm112(j)
			b= 0.5d0*(bi+sepsgm16(j))

			endif
			
			if ((slesid1(i).ne.0) .and.
     *			    (slesid1(i) .eq.slesid1(j)))  then
			   a = a*tmp
			   q = q*tmp

c       ileana
			   if(.not.arith)then
			   b = b*tmp

			   endif

			end if

			s6=s2*s2*s2
			a = a*s6*s6
			b = b*s6
			e1 = a - b
			df1 = -6.0d0*s2*(a+e1)
			se_vdw1 = se_vdw1 + e1 

100			continue
			df = df1 + df2

			rx = df*rx
			ry = df*ry
			rz = df*rz
			dxi = dxi + rx
			dyi = dyi + ry
			dzi = dzi + rz
			dx1(j) = dx1(j) - rx
			dy1(j) = dy1(j) - ry
			dz1(j) = dz1(j) - rz
			se_el1 = se_el1 + e2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                   if (newpoit(j) .eq. imb1(n)) then
                    factr(j)=1
                    else if (newpoit(j) .eq. imb2(n)) then
                     factr(j)=-1
                    else
                     factr(j)=0
                   end if
		   xtemp1(j) = -f(n)*rx
		   ytemp1(j) = -f(n)*ry
		   ztemp1(j) = -f(n)*rz
		   xtemp1(i) = f(n)*rx
		   ytemp1(i) = f(n)*ry
		   ztemp1(i) = f(n)*rz
                   dfx(n)=tempdfx(n)*(e1+e2)
                   dfy(n)=tempdfy(n)*(e1+e2)
                   dfz(n)=tempdfz(n)*(e1+e2)
                   tempdx1(j) = xtemp1(j) + factr(j) * dfx(n)
                   tempdy1(j) = ytemp1(j) + factr(j) * dfy(n)
                   tempdz1(j) = ztemp1(j) + factr(j) * dfz(n)
                   tempdx1(i) = xtemp1(i) + factr(i) * dfx(n)
                   tempdy1(i) = ytemp1(i) + factr(i) * dfy(n)
                   tempdz1(i) = ztemp1(i) + factr(i) * dfz(n)
                   sdx(j)=sdx(j)+tempdx1(j)
                   sdy(j)=sdy(j)+tempdy1(j)
                   sdz(j)=sdz(j)+tempdz1(j)
                   sdx(i)=sdx(i)+tempdx1(i)
                   sdy(i)=sdy(i)+tempdy1(i)
                   sdz(i)=sdz(i)+tempdz1(i)

200		continue
		end if
		jbeg11(n) = jend1 + 1
                dx1(i) = dx1(i) + dxi
                dy1(i) = dy1(i) + dyi
                dz1(i) = dz1(i) + dzi
400		continue
  

C       now calculate for the 1-4 special interactions

                 do 300 i=lz14_1(n-1)+1,lz14_1(n)
c@                 do 300 i=lz14_1(n-1)+1,-1
                        j = sspec11(i)
			k = sspec12(i)
			j1 = pt_to_spcl(j)
			k1 = pt_to_spcl(k)
                        rx = coor(1,j) - coor(1,k)
                        ry = coor(2,j) - coor(2,k)
                        rz = coor(3,j) - coor(3,k)
                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2
                        s6=s2*s2*s2
                        a = s1p14(1,i)*s6*s6
                        b = s1p14(2,i)*s6
                        e1 = (a - b)
                        df1 = -6.0d0*s2*(a+e1)
			if (dabs(s1p14(3,i)).gt.1.d-12) then
                        	s   = dsqrt(s2)
                        	e2  = (s1p14(3,i)*s)
                        	df2 = -e2*s2
			else
				e2  = 0.d0
				df2 = 0.d0
			end if
                        df = df1 + df2
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        se_vdw1 = se_vdw1 + e1
                        se_el1 = se_el1 + e2
                   	if (j .eq. imb1(n)) then
                    		factr(j1)=1
                    	else if (j .eq. imb2(n)) then
                     		factr(j1)=-1
                    	else
                     		factr(j1)=0
                   	end if
		   	if (k .eq. imb1(n)) then
		    		factr(k1)=1
		   	else if (k .eq. imb2(n)) then
		    		factr(k1)=-1
		   	else
		    		factr(k1)=0
		   	end if
                        dx1(j1) = dx1(j1) + rx
                        dy1(j1) = dy1(j1) + ry
                        dz1(j1) = dz1(j1) + rz
                        dx1(k1) = dx1(k1) - rx
                        dy1(k1) = dy1(k1) - ry
                        dz1(k1) = dz1(k1) - rz
		    xtemp1(k1) = -f(n)*rx
		    ytemp1(k1) = -f(n)*ry
		    ztemp1(k1) = -f(n)*rz
		    xtemp1(j1) = f(n)*rx
		    ytemp1(j1) = f(n)*ry
		    ztemp1(j1) = f(n)*rz
                   dfx(n)=tempdfx(n)*(e1+e2)
                   dfy(n)=tempdfy(n)*(e1+e2)
                   dfz(n)=tempdfz(n)*(e1+e2)
                    tempdx1(j1) = xtemp1(j1) + factr(j1) * dfx(n)
                    tempdy1(j1) = ytemp1(j1) + factr(j1) * dfy(n)
                    tempdz1(j1) = ztemp1(j1) + factr(j1) * dfz(n)
                   tempdx1(k1) = xtemp1(k1) + factr(k1) * dfx(n)
                   tempdy1(k1) = ytemp1(k1) + factr(k1) * dfy(n)
                   tempdz1(k1) = ztemp1(k1) + factr(k1) * dfz(n)
                   sdx(j1)=sdx(j1)+tempdx1(j1)
                   sdy(j1)=sdy(j1)+tempdy1(j1)
                   sdz(j1)=sdz(j1)+tempdz1(j1)
                   sdx(k1)=sdx(k1)+tempdx1(k1)
                   sdy(k1)=sdy(k1)+tempdy1(k1)
                   sdz(k1)=sdz(k1)+tempdz1(k1)
300              continue

        spoint21(0) = 0
        initial = poitype2(n-1)+1
	jbeg21(n)=spoint21(poitype2(n-1))+1
c
c All nonbonded interactions for the second curve in the next loops
c
	do 1400 i=initial,poitype2(n)-1

                   if (newpoit(i) .eq. imb1(n)) then
                    factr(i)=1
                    else if (newpoit(i) .eq. imb2(n)) then
                     factr(i)=-1
                    else
                     factr(i)=0
                   end if

		jend1 = spoint21(i)
		tmp = 1.d0/sptwei2(i)
		ai  = sepsgm212(i)*pick

c       ileana

		if(.not.arith)then
		bi  = sepsgm26(i)*pick
		else
		bi  = sepsgm26(j)

		endif

		qi  = sptchg2(i)*epstmp
		xi  = coor(1,newpoit(i))
		yi  = coor(2,newpoit(i))
		zi  = coor(3,newpoit(i))
		dxi = 0.d0
		dyi = 0.d0
		dzi = 0.d0

		if (jbeg21(n).le.jend1) then

		do 1200 k=jbeg21(n),jend1
			j=slist21(k)
			rx = xi - coor(1,newpoit(j))
			ry = yi - coor(2,newpoit(j))
			rz = zi - coor(3,newpoit(j))
			r2=rx*rx+ry*ry+rz*rz
			s2=1.0d0/r2

			q = qi*sptchg2(j)
			s = dsqrt(s2)
			e2 = q*s
			df2 = -e2*s2
			if (r2.gt.cutvdw2) then
c then only electrostatic should be calculated
				df1 = 0.d0
				if ((slesid2(i).ne.0) .and.
     *			       (slesid2(i) .eq.slesid2(j))) q = q*tmp
				go to 1100
			end if


c       ileana

			if(.not.arith)then
			a= ai*sepsgm212(j)
			b= bi*sepsgm26(j)

			else

			a= ai*sepsgm212(j)
			b= 0.5d0*(bi+sepsgm26(j))

			endif

			if ((slesid2(i).ne.0) .and.
     *			    (slesid2(i) .eq.slesid2(j)))  then
			   a = a*tmp
			   q = q*tmp
			   
			   if(.not.arith)then
			      b= b*tmp
			      endif
			   
			end if

			s6=s2*s2*s2
			a = a*s6*s6
			b = b*s6
			e1 = a - b
			df1 = -6.0d0*s2*(a+e1)
			se_vdw2 = se_vdw2 + e1 

c the electrostatic part should be done without buffering
1100			continue
			df = df1 + df2

			rx = df*rx
			ry = df*ry
			rz = df*rz
			dxi = dxi + rx
			dyi = dyi + ry
			dzi = dzi + rz
			se_el2 = se_el2 + e2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                   if (newpoit(j) .eq. imb1(n)) then
                    factr(j)=1
                    else if (newpoit(j) .eq. imb2(n)) then
                     factr(j)=-1
                    else
                     factr(j)=0
                   end if
		   xtemp2(j) = -(1-f(n))*rx
		   ytemp2(j) = -(1-f(n))*ry
		   ztemp2(j) = -(1-f(n))*rz
		   xtemp2(i) = (1-f(n))*rx
		   ytemp2(i) = (1-f(n))*ry
		   ztemp2(i) = (1-f(n))*rz
                   dfx(n)=tempdfx(n)*(e1+e2)
                   dfy(n)=tempdfy(n)*(e1+e2)
                   dfz(n)=tempdfz(n)*(e1+e2)
                   tempdx2(j) = xtemp2(j) - factr(j) * dfx(n)
                   tempdy2(j) = ytemp2(j) - factr(j) * dfy(n)
                   tempdz2(j) = ztemp2(j) - factr(j) * dfz(n)
                   tempdx2(i) = xtemp2(i) - factr(i) * dfx(n)
                   tempdy2(i) = ytemp2(i) - factr(i) * dfy(n)
                   tempdz2(i) = ztemp2(i) - factr(i) * dfz(n)
                   sdx(j)=sdx(j)+tempdx2(j)
                   sdy(j)=sdy(j)+tempdy2(j)
                   sdz(j)=sdz(j)+tempdz2(j)
                   sdx(i)=sdx(i)+tempdx2(i)
                   sdy(i)=sdy(i)+tempdy2(i)
                   sdz(i)=sdz(i)+tempdz2(i)

1200		continue
		end if
		jbeg21(n) = jend1 + 1
		dx1(i) = dx1(i) + dxi
		dy1(i) = dy1(i) + dyi
		dz1(i) = dz1(i) + dzi
1400	continue
  
                 do 1300 i=lz14_2(n-1)+1,lz14_2(n)
                        j = sspec21(i)
			k = sspec22(i)
			j1 = pt_to_spcl(j)
			k1 = pt_to_spcl(k)
                        rx = coor(1,j) - coor(1,k)
                        ry = coor(2,j) - coor(2,k)
                        rz = coor(3,j) - coor(3,k)
                        r2=rx*rx+ry*ry+rz*rz
                        s2=1.0d0/r2
                        s6=s2*s2*s2
                        a = a*s6*s6
                        b = b*s6
                        e1 = (a - b)
                        df1 = -6.0d0*s2*(a+e1)
			if (dabs(s1p14(3,i)).gt.1.d-12) then
			 s   = dsqrt(s2)
			 e2  = s2p14(3,i)*s
			 df2 = -e2*s2
			else
			 e2  = 0.d0
			 df2 = 0.d0
			end if
                        df = df1 + df2
                        rx = df*rx
                        ry = df*ry
                        rz = df*rz
                        dxi = dxi + rx
                        dyi = dyi + ry
                        dzi = dzi + rz
                        se_vdw2 = se_vdw2 + e1
                        se_el2 = se_el2 + e2
                   if (j .eq. imb1(n)) then
                    factr(j1)=1
                   else if (j .eq. imb2(n)) then
                     factr(j1)=-1
                   else
                     factr(j1)=0
                   end if
                   if (k .eq. imb1(n)) then
                    factr(k1)=1
                   else if (k .eq. imb2(n)) then
                     factr(k1)=-1
                   else
                     factr(k1)=0
                   end if
		   xtemp2(j1) = -(1-f(n))*rx
		   ytemp2(j1) = -(1-f(n))*ry
		   ztemp2(j1) = -(1-f(n))*rz
		   xtemp2(k1) = (1-f(n))*rx
		   ytemp2(k1) = (1-f(n))*ry
		   ztemp2(k1) = (1-f(n))*rz
                   dfx(n)=tempdfx(n)*(e1+e2)
                   dfy(n)=tempdfy(n)*(e1+e2)
                   dfz(n)=tempdfz(n)*(e1+e2)
                   tempdx2(j1) = xtemp2(j1) - factr(j1) * dfx(n)
                   tempdy2(j1) = ytemp2(j1) - factr(j1) * dfy(n)
                   tempdz2(j1) = ztemp2(j1) - factr(j1) * dfz(n)
                   tempdx2(k1) = xtemp2(k1) - factr(k1) * dfx(n)
                   tempdy2(k1) = ytemp2(k1) - factr(k1) * dfy(n)
                   tempdz2(k1) = ztemp2(k1) - factr(k1) * dfz(n)
                   sdx(j1)=sdx(j1)+tempdx2(j1)
                   sdy(j1)=sdy(j1)+tempdy2(j1)
                   sdz(j1)=sdz(j1)+tempdz2(j1)
                   sdx(k1)=sdx(k1)+tempdx2(k1)
                   sdy(k1)=sdy(k1)+tempdy2(k1)
                   sdz(k1)=sdz(k1)+tempdz2(k1)
1300              continue
        if (debug) then
         print*,'se_el1=',se_el1
         print*,'se_el2=',se_el2
         print*,'se_vdw1=',se_vdw1
         print*,'se_vdw2=',se_vdw2
         print*,'spoint11=',spoint11(poitype1(n)-1)
         print*,'spoint21=',spoint21(poitype2(n)-1)
        end if

	return 
	end
