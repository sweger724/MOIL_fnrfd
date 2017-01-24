      function of()
c
c Open File: a routine to open a file, does the following
c (i) extract a file name and check that file exists
c (ii) extract a unit number and return the unit number as "of"
c
      integer of
      integer geti
      
      include 'COMMON/LINE.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      character*80 getchar,namef,def
      integer lc,level
      logical exist,find

c get file name
c
      namef = ' '
      lc = 79
      def = 'none'
      namef = getchar('name',def,lc)
      exist = .false.
      if (debug) then
         write(stdo,*)' name of file ',namef(1:lc)
         write(stdo,*)' length of name ',lc
      end if
      inquire(file=namef(1:lc),exist=exist)
c
      of = geti('unit',-1)
      if (of.eq.-1) call getunit(of)
c
      if (find('read')) then
         if (of.eq.5) then
            write(stdo,100)
            return
         end if
         if (of.eq.-1) then
            of = stdi
            return
         else if (of.eq.stdi) then
            write(stdo,100)
 100        format(/,1x,'Input taken from standard file',/)
            return
         else if (find('bina')) then
            if (.not.exist) then
               level = 1
               call alert('of',2,'READ file missing',17,level)
               return
            end if
            open(unit=of,file=namef(1:lc),form='unformatted',
     1           status='old')
            rewind of
            return
         else
            if (.not.exist) then
               level = 1
               call alert('of',2,'READ file missing',17,level)
               return
            end if
            open (unit=of,file=namef(1:lc),status='old')
            rewind of
            return
         end if
c
      else if (find('wovr')) then
         if (of.eq.stdo) then
            level = 1
            call alert('of',2,'Missing unit, set to stdo',25,level)
            return
         end if
         if (find('bina')) then
            open(unit=of,file=namef(1:lc),form='unformatted',
     1           status='unknown')
            rewind of
            return
         else
            open (unit=of,file=namef(1:lc),status='unknown')
            rewind of
            return
         end if
c
      else if (find('writ')) then
         if (exist) then
            write(stdo,101)
 101        format(1x,//,' *** ILLEGAL ATTEMPT TO OPEN FILE ',//,
     1           1x,'*** IF EXISTS MUST SPECIFY read OR wovr ',//)
            level = 1
            call alert('of',2,' Illegal Open File statement',28,level)
         end if
         if (find('bina')) then
            open(unit=of,file=namef(1:lc),form='unformatted',
     1           status='new')
            rewind of
         else
            open (unit=of,file=namef(1:lc),status='new')
            rewind of
         end if
      end if
      
      return
      end



