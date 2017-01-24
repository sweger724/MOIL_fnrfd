      program extract

C     extract certain first nslc particles of a path file and write to
C     another  path file

      implicit none

      include 'COMMON/LINE.BLOCK'
      include 'COMMON/UNITS.BLOCK'

      character*7 name
      integer geti, namel, i, j, of, nslc, nprt, urcrd, uwcrd
      data urcrd,uwcrd/2*99/
      double precision e0, r(3,10000)
      logical find

      jnkf      = 25
      stdi      = 5
      stdo      = 6
      stderr    = 0
      name      = 'extract'
      namel     = 7

      urcrd     = 0
      uwcrd     = 0

      nprt      = 0
      nslc      = 0

      open ( unit=jnkf, status='scratch' )

 1    call rline(name,namel,stdi)
      if (find('file')) then
         if (find('rcrd')) urcrd = of()
         if (find('wcrd')) uwcrd = of()
      end if
      nprt = geti('#prt',nprt)
      nslc = geti('#slc',nslc)
      if (find('acti')) go to 2
      go to 1


 2    write(*,*) 'urcrd, uwcrd, nprt, nslc: ', urcrd, uwcrd, nprt, nslc

      read(urcrd) e0, ((r(j,i),i=1,nprt),j=1,3)
      write(uwcrd) e0, ((r(j,i),i=1,nslc),j=1,3)

      stop
      end
