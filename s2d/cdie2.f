	subroutine cdie2()

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
	include 'COMMON/GETVEC.BLOCK'
	include 'COMMON/SCNDRV.BLOCK'

	double precision pick 
	double precision epstmp
	double precision rx,ry,rz,r2,s2,a,b,q,ai,bi
	double precision s,s6,tmp,b6,b12
	double precision s8,s3,abq1,abq2,add
	integer index,jndex,kndex
	double precision qi

	integer i,j,k,jbeg,jend


	if (eelyes) then
		epstmp = 332.0716d0/eps
	else
		epstmp =  0.0d0
	end if

	if (evdyes) then
		pick = 1.0d0
	else
		pick = 0.0d0
	end if


	jbeg=1

	do 100 i=1,npt-1


		jend = point1(i)


		if (jbeg.le.jend) then

		tmp = 1.d0/ptwei(i)
		a   = epsgm12(i)*pick
                if (arith) then
		b   = epsgm6(i)
                else
                b   = epsgm6(i)*pick
                endif
		q   = ptchg(i)*epstmp
		do 200 k=jbeg,jend
			j=list1(k)
			rx = coor(1,i) - coor(1,j)
			ry = coor(2,i) - coor(2,j)
			rz = coor(3,i) - coor(3,j)
			r2=rx*rx+ry*ry+rz*rz
			if (r2.gt.cutvdw2 .and.
     1				(.not.nocut)) then
				ai = 0.d0
				bi = 0.d0
			else 
				ai= a*epsgm12(j)
                           if (arith) then
				bi= b+epsgm6(j)
                           else
                                bi= b*epsgm6(j)
                           endif
			end if
			qi= q*ptchg(j)
			if (lesid(i).ne.0) then
			  if (lesid(i).eq.lesid(j)) then
			   ai = ai*tmp
			   bi = bi*tmp
			   qi = qi*tmp
			  end if
			end if
			s2=1.0d0/r2
			s6=s2*s2*s2
			s8=s2*s6
			s = dsqrt(s2)
			s3 = s2*s
                        if (arith) then
                        b6 = bi*bi*bi*bi*bi*bi
                        b12 = b6*b6
			abq1=168.d0*ai*b12*s8*s8-48.d0*ai*b6*s8*s2
     1				+3.d0*qi*s3*s2
			abq2=-12.d0*ai*b12*s6*s8+6.d0*b6*ai*s8-qi*s3
                        else
                        abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1                          +3.d0*qi*s3*s2
                        abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3
                        endif
c			write(*,*)'#1 abq1 abq2 ',abq1,abq2
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
			add        = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
			add        = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
			add        = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_1(kndex) = -add
200		continue
		end if
		jbeg = jend + 1
100	continue

        jbeg=1
  
c start SECOND loop including particles with cutoff
c cutvdw2 - cutoff appropriate for van der Waals forces
c includes particles with zero charges and therefore
c NO ELECTROSTATIC
c
       do 105 i=1,npt-1


              jend = point2(i)


              if (jbeg.le.jend) then
              tmp = 1.d0/ptwei(i)
              do 205 k=jbeg,jend
                      j=list2(k)
                      rx = coor(1,i) - coor(1,j)
                      ry = coor(2,i) - coor(2,j)
                      rz = coor(3,i) - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutvdw2 .and.
     1                                (.not.nocut)) go to 205
                        a = epsgm12(i)*epsgm12(j)*pick
                        if (arith) then
                        b = epsgm6(i)+epsgm6(j)
                        else
                        b = epsgm6(i)*epsgm6(j)*pick
                        endif
                      if (lesid(i).ne.0) then
                        if (lesid(i).eq.lesid(j)) then
                         a = a*tmp
                         b = b*tmp
                        end if
                      end if
			s2=1.0d0/r2
			s6=s2*s2*s2
			s8=s2*s6
                        if (arith) then
                        b6 = b*b*b*b*b*b
                        b12 = b6*b6
			abq1=168.d0*a*b12*s8*s8-48.d0*a*b6*s8*s2
			abq2=-12.d0*a*b12*s6*s8+6.d0*b6*a*s8
                        else
                        abq1=168.d0*a*s8*s8-48.d0*b*s8*s2
                        abq2=-12.d0*a*s6*s8+6.d0*b*s8
                        endif
c			write(*,*)' i point(i) ',i,point2(i)
c			write(*,*)'#2 abq1 abq2 ',abq1,abq2
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
			add         = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
			add         = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
			add         = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_2(kndex) = -add
205           continue
              end if
              jbeg = jend + 1
105    continue


       jbeg=1

c start THIRD loop including particles with upper cutoff
c cutele2 - cutoff appropriate for electrostics - and lower
c cutoff - cutvdw2 - the lower cutoff electrostatic was already
c calculated using the van der Waals (first) loop for charged
c particles. Includes ONLY electrostic forces
c
       do 110 i=1,npt-1


              jend = point3(i)


              if (jbeg.le.jend) then

              tmp = 1.d0/ptwei(i)
              do 210 k=jbeg,jend
                      j=list3(k)
                      rx = coor(1,i) - coor(1,j)
                      ry = coor(2,i) - coor(2,j)
                      rz = coor(3,i) - coor(3,j)
                      r2=rx*rx+ry*ry+rz*rz
                      if (r2.gt.cutele2 .and.
     1                                (.not.nocut)) go to 210

                      qi=ptchg(i)*ptchg(j)*epstmp
                      if (lesid(i).ne.0 .and.
     1                          lesid(i).eq.lesid(j)) qi = qi*tmp

			s2=1.0d0/r2
			s = dsqrt(s2)
			s3 = s2*s
			abq1=3.d0*qi*s3*s2
			abq2=-qi*s3
c			write(*,*)'#3 abq1 abq2 ',abq1,abq2
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
			add         = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
			add         = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
			add         = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2nb_3(kndex) = -add
210           continue
              end if
              jbeg = jend + 1
110    continue



C	now calculate for the 1-4 special interactions


	do 400 k=1,totspe
		i = spec1(k)
		j = spec2(k)
		rx = coor(1,i) - coor(1,j)
		ry = coor(2,i) - coor(2,j)
		rz = coor(3,i) - coor(3,j)
		r2=rx*rx+ry*ry+rz*rz
		s2=1.0d0/r2
		s6=s2*s2*s2
		s8=s2*s6
		s = dsqrt(s2)
		s3 = s2*s
		if (p14(1,k).gt.1.d-12) then
		 abq1=168.d0*p14(1,k)*s8*s8-48.d0*p14(2,k)*s8*s2
     1				+3.d0*p14(3,k)*s3*s2
		 abq2=-12.d0*p14(1,k)*s6*s8+6.d0*p14(2,k)*s8-p14(3,k)*s3
		else
		 abq1=3.d0*p14(3,k)*s3*s2
		 abq2=-p14(3,k)*s3
		end if
		index = 6*(i-1)+1
		jndex = 6*(j-1)+1
		kndex = 6*(k-1)+1
		add         = abq1*rx*rx + abq2
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
		add         = abq1*ry*rx
		index = index + 1
		jndex = jndex + 1
		kndex = kndex + 1
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
		add         = abq1*rz*rx 
		index = index + 1
		jndex = jndex + 1
		kndex = kndex + 1
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
		add         = abq1*ry*ry + abq2
		index = index + 1
		jndex = jndex + 1
		kndex = kndex + 1
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
		add         = abq1*rz*ry
		index = index + 1
		jndex = jndex + 1
		kndex = kndex + 1
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
		add         = abq1*rz*rz + abq2
		index = index + 1
		jndex = jndex + 1
		kndex = kndex + 1
		diag(index) = diag(index) + add
		diag(jndex) = diag(jndex) + add
		d2spec(kndex) = -add
400	continue

	return 
	end
