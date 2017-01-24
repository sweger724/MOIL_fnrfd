	character*(*) function getchar(expr,def,lcpext)
c check if command line (line) includes a key word (expr). If yes
c get getc, if not set getc to def
c Dec, 23 1993 RE: Add null (\0) to end of character for compatability
c with C
c
	integer lcpext
	character*4 expr
	character*1 BLANK
	character*80 def
	include 'COMMON/LINE.BLOCK'
	include 'COMMON/DEBUG.BLOCK'
	include 'COMMON/UNITS.BLOCK'
c
c local
c
	integer i,j,k,lcp
	DATA BLANK/' '/
c
c verify that the line exists ...
c
	if (point(100).le.0 .or. nexp.le.0 ) then
		call alert('getchar',7,'No input line',13,0)
		getchar(1:lcpext)=def(1:lcpext)
		return
	end if
        j=1
	do 2 i=1,nexp
		if(expr(1:4).eq.line(j:j+3)) then
		 do 1 k=j,point(i)
c
c check if this is a substitution command
c
		  if (line(k:k).eq.'=') then
                     getchar(1:point(i)-k-2)=line(k+2:point(i)-1)
		     line(j:point(i))=BLANK
	     	     lcp=point(i)-k-2
		     if (lcp.ne.lcpext) lcpext = lcp
c		     if (lcp.lt.lcpext) getchar(lcp+1:lcpext)=BLANK
c		     if (lcpext.eq.0) lcpext = lcp
		     getchar(lcpext+1:lcpext+1) = '\0'
		     return
                  end if
1		 continue
		 call alert('getchar',7,'Missing =',9,0)
		 getchar(1:lcpext)=def(1:lcpext)
		 getchar(lcpext+1:lcpext+1)='\0'
		 return
		end if
		j=point(i)+1
2	continue
	getchar(1:lcpext)=def(1:lcpext)
	getchar(lcpext+1:lcpext+1)='\0'
	return
	end
