	subroutine init_ssbp(ipick)
	include 'COMMON/LENGTH.BLOCK'
	include 'COMMON/CONNECT.BLOCK'
	include 'COMMON/SSBP.BLOCK'
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/UNITS.BLOCK'
	include 'COMMON/PARALLEL.BLOCK'

	integer ipick(*)

c local
	integer i,level,namel
	integer geti
	logical find
	double precision getd
	character*9 name
	data namel/9/
	data name/'init_ssbp'/
Cdeb+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Spherical Solvent Boundary Potential (SSBP)
C
C Adapted from the original code by
C       Dmitrii Beglovd and Benoit Roux  (1994)
C           Department of Chemistry, University of Montreal
C
C Default is include all energy terms for the
C    solvent boundary potential
C----------------------------------------------------------

      QSSBP = .TRUE.
      QFIX  = .FALSE.
      QKIRK = .TRUE.
      QANGU = .TRUE.
      QCAVI = .TRUE.
      QHSR  = .TRUE.
      QHARM = .FALSE.

      nmult = geti('nmul',15)
      diecst = getd('diec',78.4d0)

      write(stdo,'(a,f6.3)')
     & ' Kirkwood dielectric reaction field, dielectric constant ',
     &   DIECST

      drmax1 = getd('drdi',2.8D0)
      drmax2 = getd('drca',2.6D0)

      presi  = getd('pres',1.D0)
      stens  = getd('surt',0.033D0)
C
C-----The following pick is for selecting the TIP3 oxygens
C     The cavity potential is applied on TIP3 oxygens only.
C
      call pick(ipick,i)
      NTSSBP=0
      DO 22 I=1,npt
      IF(ipick(I).EQ.1)THEN
      NTSSBP=NTSSBP+1
      LSTSSBP(NTSSBP)=I
      ENDIF
   22 CONTINUE

	
C  Handle the case with fixed radius

      if (find('fixe')) then
        fxcas  = getd('radi',0.0D0)
        write(stdo,'(A)')
     &      ' Finite system with FIXED radius (must give the RADIUS)'
        if(fxcas.ne.0.d0)then
          write(stdo,'(A,F10.5)') ' FIXED radius of ', FXCAS
        else
          level = 1
          call alert(name,namel,'Value of fixed radius not
     1    defined',17,level)
        endif
      QFIX   = .TRUE.
      if (find('harm')) QHARM = .TRUE.
      endif

         if (qkirk  .and. find('nokirk')) qkirk   = .false.
         if (qcavi .and. find('nocavi')) qcavi  = .false.
         if (qangu .and. find('noag')) qangu  = .false.
         if (qhsr .and. find('nohsr')) qhsr  = .false.

Cdeb+++++++++++++++++++++++++++++++++++++++++++++++++++++++

	return
	end
                                                                              
