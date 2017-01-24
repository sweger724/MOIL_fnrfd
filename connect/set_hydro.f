	subroutine set_hydro()
c	This subroutine set set the cbeta and betap vector.
c	each betha carbon is given i'ts hydrophobicity acording to
c	Casari & Sippl (JMB 1992 224:725-32). 
	
C	include 'COMMON/PRLL.BLOCK'
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/COORD.BLOCK'
	include 'COMMON/UNITS.BLOCK'
C	include 'COMMON/GENERIC.BLOCK'
	include 'COMMON/ENERGY.BLOCK'
	include 'COMMON/FREEZ.BLOCK'
	
	integer i

	avg_hydro = 0.36d0
	nbeta = 0
c	r0_hydro = 23.d0
	do 1 i = 1,npt
		if ((moname(poimon(i)) .eq. 'GLY1') .or.
     *		    (moname(poimon(i)) .eq. 'NGL2')) then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				cbeta(nbeta) = 0.3d0
				betap(nbeta) = i
			end if
		else if ((moname(poimon(i)) .eq. 'CYS') .or.
     *		    (moname(poimon(i)) .eq. 'CYX')) then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				cbeta(nbeta) = 2.3d0
				betap(nbeta) = i
			end if
		else if (moname(poimon(i)) .eq. 'TRP') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				cbeta(nbeta) = 1.9d0
				betap(nbeta) = i
			end if
		else if (moname(poimon(i)) .eq. 'ILE') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				cbeta(nbeta) = 1.8d0
				betap(nbeta) = i
			end if
		else if (moname(poimon(i)) .eq. 'PHE') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 1.3d0
			end if
		else if (moname(poimon(i)) .eq. 'VAL') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 1.d0
			end if
		else if (moname(poimon(i)) .eq. 'TYR') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 0.9d0
			end if
		else if (moname(poimon(i)) .eq. 'LEU') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 0.9d0
			end if
		else if ((moname(poimon(i)) .eq. 'MET') .or.
     *		    (moname(poimon(i)) .eq. 'CMET')) then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 0.8d0
			end if
		else if ((moname(poimon(i)) .eq. 'HIS') .or.
     *		    (moname(poimon(i)) .eq. 'HIP')) then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 0.8d0
			end if
		else if (moname(poimon(i)) .eq. 'ALA') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = 0.5d0
			end if
		else if (moname(poimon(i)) .eq. 'THR') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.1d0
			end if
		else if (moname(poimon(i)) .eq. 'ASN') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.1d0
			end if
		else if (moname(poimon(i)) .eq. 'ARG') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.3d0
			end if
		else if (moname(poimon(i)) .eq. 'PRO') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.6d0
			end if
		else if (moname(poimon(i)) .eq. 'SER') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.4d0
			end if
		else if (moname(poimon(i)) .eq. 'GLN') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.7d0
			end if
		else if (moname(poimon(i)) .eq. 'GLU') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -0.9d0
			end if
		else if (moname(poimon(i)) .eq. 'ASP') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -1.d0
			end if
		else if (moname(poimon(i)) .eq. 'LYS') then
			if (ptnm(i)(1:4).eq.'CB') then
				nbeta = nbeta + 1
				betap(nbeta) = i
				cbeta(nbeta) = -1.2d0
			end if
		end if

1	continue
	if ((moname(1) .eq. 'NTER') .and. 
     *      (poimon(betap(1)) .eq. 2)) cbeta(1) = cbeta(1) - 1.d0
	if (moname(poimon(betap(nbeta)) + 1) .eq. 'CTER') 
     *		cbeta(nbeta) = cbeta(nbeta) - 1.d0
	return
	end
