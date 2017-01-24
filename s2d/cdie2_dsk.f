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
	include 'COMMON/SCNDRV.BLOCK'

	double precision pick 
	double precision epstmp
	double precision rx,ry,rz,r2,s2,a,b,q,ai,bi
	double precision s,s6,tmp
	double precision s8,s3,abq1,abq2,add
	integer index,jndex,kndex
	double precision qi

	integer i,j,k,jbeg,jend

c	<<****************** added by jing ***********************
c
c	For using two data files act as two big array d2nb_1()
c	and d2nb_2() respectively. All stuff added by jing are
c	included in two "****************" comment lines, and old
c	statement will add a "c*" in each of them.

	integer fid1, fid2

	parameter (fid1=30)
	parameter (fid2=31)


	double precision d1, d2, d3, d4, d5, d6

c	open two data files use fid1, fid2 as units respectively
	open (unit=fid1, file='d2nb_1.data', access='sequential',
     $ 	      form='unformatted', status='unknown')
	open (unit=fid2, file='d2nb_2.data', access='sequential',
     $ 	      form='unformatted', status='unknown')

c	********************************************************>>

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

	e_vdw = 0.0d0
	e_el=0.d0

	jbeg=1

	do 100 i=1,npt-1


		jend = point1(i)


		if (jbeg.le.jend) then

		tmp = 1.d0/ptwei(i)
		a   = epsgm12(i)*pick
		b   = epsgm6(i)*pick
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
				bi= b*epsgm6(j)
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
			abq1=168.d0*ai*s8*s8-48.d0*bi*s8*s2
     1				+3.d0*qi*s3*s2
			abq2=-12.d0*ai*s6*s8+6.d0*bi*s8-qi*s3
c			write(*,*)'#1 abq1 abq2 ',abq1,abq2
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*    			d2nb_1(kndex) = -add

c	<<*******************************************************
	  		d1 = -add
c	*******************************************************>>

			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_1(kndex) = -add

c	<<*******************************************************
			d2 = -add
c	*******************************************************>>

			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_1(kndex) = -add

c	<<*******************************************************
			d3 = -add
c	*******************************************************>>

			add        = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_1(kndex) = -add

c	<<*******************************************************
			d4 = -add
c	*******************************************************>>

			add        = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_1(kndex) = -add

c	<<*******************************************************
			d5 = -add
c	*******************************************************>>

			add        = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_1(kndex) = -add

c	<<*******************************************************
			d6 = -add

c	 		write d1-6 to the data file d2nb_1.data
			write (fid1) d1, d2, d3, d4, d5, d6

c	*******************************************************>>

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
                        b = epsgm6(i)*epsgm6(j)*pick
                      if (lesid(i).ne.0) then
                        if (lesid(i).eq.lesid(j)) then
                         a = a*tmp
                         b = b*tmp
                        end if
                      end if
			s2=1.0d0/r2
			s6=s2*s2*s2
			s8=s2*s6
			abq1=168.d0*a*s8*s8-48.d0*b*s8*s2
			abq2=-12.d0*a*s6*s8+6.d0*b*s8
c			write(*,*)' i point(i) ',i,point2(i)
c			write(*,*)'#2 abq1 abq2 ',abq1,abq2
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d1 = -add
c	***********************************************************>>

			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d2 = -add
c	***********************************************************>>

			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d3 = -add
c	***********************************************************>>

			add         = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d4 = -add
c	***********************************************************>>

			add         = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d5 = -add
c	***********************************************************>>

			add         = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
c*			d2nb_2(kndex) = -add

c	<<***********************************************************
			d6 = -add

c	write d1-d6 into data file d2nb_2.data
			write (fid2) d1, d2, d3, d4, d5, d6

c	***********************************************************>>

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
			add         = abq1*rz*rz
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
		abq1=168.d0*p14(1,k)*s8*s8-48.d0*p14(2,k)*s8*s2
     1				+3.d0*p14(3,k)*s3*s2
		abq2=-12.d0*p14(1,k)*s6*s8+6.d0*p14(2,k)*s8-p14(3,k)*s3
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


