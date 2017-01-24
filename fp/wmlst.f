	subroutine wmlst(ucrd,pp,pn,n,small)

c       Write pp, then write pp + small*pn (path format)
c
c       INPUT:
c       ucrd - file to write to
c       pp - milestone coordinate ("plane point")
c       pn - milestone normal ("plane normal")
c       n - number of atoms per structure
c       small - amount of pn to add to pp

	implicit none
	integer ucrd, n, i, k
	double precision small, e, pn(3,*), pp(3,*)

	e = 999.888
	write(ucrd) e, ((pp(k,i),i=1,n),k=1,3)
	write(ucrd) e, ((pp(k,i)+small*pn(k,i),i=1,n),k=1,3)

	return
	end
