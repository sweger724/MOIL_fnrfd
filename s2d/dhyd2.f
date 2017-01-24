	subroutine dhyd2()
c---Written by M. Kara-Ivanov 10.01.94 ; principle of cdie.f was modified 
c--- for Hydrophobic interactions------------
      	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/NBLIST.BLOCK'
      	include 'COMMON/CONNECT.BLOCK'
      	include 'COMMON/ENERGY.BLOCK'
      	include 'COMMON/COORD.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/SPECL.BLOCK'
	include 'COMMON/GETVEC.BLOCK'
	include 'COMMON/SCNDRV.BLOCK'

        
        double precision les_scale,avg2
        double precision rx,ry,rz,r2,s,s2,h,hydro_th2
	double precision abq1,abq2,add
	integer index,jndex,kndex
	integer i,j,k,l,m

	do 10 i=1,6*ibeta2
	   d2hyd(i)=0
 10	continue
c----- cutoff for hydrophobic interactions------
        hydro_th2 = hydro_th * hydro_th
        avg2      = 2.d0*avg_hydro
c-- l and m are indexes for beta atoms!        
c--- Remember! Written in this way it should work only for rather short
c--- peptides: size: (6*nbeta*(nbeta-1))/2. Alternative way would be to 
c--- to introduce massive of types of beta atoms!!!!
        k=0
        do 100 l=1,nbeta - 1
           i = betap(l)
           les_scale = ptwei(i)
           do 200 m = l + 1,nbeta
                k=k+1
                j = betap(m)
                if ((lesid(i).ne.0) .and.
     *              (lesid(i) .eq.lesid(j)))  then
                    if (cplbl(i) .ne. cplbl(j)) go to 200
                end if
                if(lesid(i).ne.lesid(j))
     *                  les_scale = les_scale*ptwei(j)
c-----This is potential part-------------------------
                h = (cbeta(l)+cbeta(m)+avg2)*les_scale
                rx = coor(1,i) - coor(1,j)
                ry = coor(2,i) - coor(2,j)
                rz = coor(3,i) - coor(3,j)
                r2=rx*rx+ry*ry+rz*rz

                if (r2 .gt. hydro_th2) go to 200
			s2=1.0d0/r2
			s = dsqrt(s2)
                        abq1= -hydro_scale*h*s*s2
                        abq2= hydro_scale*h*s
			index = 6*(i-1)+1
			jndex = 6*(j-1)+1
			kndex = 6*(k-1)+1
			add         = abq1*rx*rx + abq2
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
			add         = abq1*ry*rx
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
			add         = abq1*rz*rx 
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
			add        = abq1*ry*ry + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
			add        = abq1*rz*ry
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
			add        = abq1*rz*rz + abq2
			index = index + 1
			jndex = jndex + 1
			kndex = kndex + 1
			diag(index) = diag(index) + add
			diag(jndex) = diag(jndex) + add
			d2hyd(kndex) = -add
200		continue
100	continue

	return 
	end


