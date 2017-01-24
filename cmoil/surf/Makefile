# Compilation flags
#CC	    = cc
CC	    = gcc
INCLUDE     = -I. 
#LINCLUDE    = -lcurses -ltermcap -lm
LINCLUDE    = -lm
OPT_CFLAGS  = -O2 $(FLAGS) $(INCLUDE)
#CFLAGS	    = -g $(FLAGS) $(INCLUDE)
CFLAGS	    = -O2 $(FLAGS) $(INCLUDE)

# These are the user object files in the application
SRCS        = surf.c io.c compute.c dual.c utils.c lp.c chull.c tessel_cases.c \
	      tessel_patches.c tessel_convex.c tessel_concave.c tessel_torus.c 

# obj files for sequential execution
OBJS        = surf.o io.o compute.o dual.o utils.o lp.o chull.o tessel_cases.o \
	      tessel_patches.o tessel_convex.o tessel_concave.o tessel_torus.o  

# make objects
surf: $(OBJS) Makefile 
	$(CC) $(CFLAGS) $(OBJS) -o surf $(LINCLUDE)
	mv surf ../../exe

lint:
	lint $(INCLUDE) $(SRCS)

clean:
	-rm -f *.o 

release:
	tar -cvf surf.tar README *.[hc] Makefile
	compress surf.tar

.c.o: 
	$(CC) $(CFLAGS) -c $*.c


# make depend makes the proper include file dependencies.  You _could_ run
# it on a sun4, but there's a bug in the SunOS version of sed that causes
# the all-important "cat makedep | sed -e 's/\.o/.d4o/p'" etc. line to
# not print the results to the stdout (but to keep pumping it through the
# edits.  The 'p' flag effectively does nothing). To avoid this I have
# replaced the 'p' flag by the 'w file' flag in sed.

DEPEND = $(SRCS)

depend: $(DEPEND)
	@ echo making dependencies...
	@ echo ' ' > makedep
	@ makedepend $(INCLUDE) -f makedep $(DEPEND)
	@ echo '/^# DO NOT DELETE THIS LINE/+2,$$d' >eddep
	@ echo '$$r makedep' >>eddep
	@ echo 'w' >>eddep
	@ cp Makefile Makefile.bak
	@ ed - Makefile < eddep
	@ rm eddep makedep 
	@ echo '  ' >> Makefile
	@ echo '# DEPENDENCIES MUST END AT END OF FILE' >> Makefile
	@ echo '# IF YOU PUT STUFF HERE IT WILL GO AWAY' >> Makefile
	@ echo '# see make depend above' >> Makefile
	@ echo 'done.'
	@ echo '  ' >> Makefile

# DO NOT DELETE THIS LINE -- make depend depends on it.
  
 
# DO NOT DELETE THIS LINE -- make depend depends on it.

surf.o: surf.h /usr/include/stdio.h /usr/include/math.h /usr/include/stdlib.h
surf.o: /usr/include/string.h /usr/include/sys/time.h
surf.o: /usr/include/sys/types.h 
surf.o: /usr/include/sys/select.h linalg.h
io.o: surf.h /usr/include/stdio.h /usr/include/math.h /usr/include/stdlib.h
io.o: /usr/include/string.h /usr/include/sys/time.h
io.o: /usr/include/sys/types.h 
io.o: /usr/include/sys/select.h linalg.h
compute.o: surf.h /usr/include/stdio.h /usr/include/math.h
compute.o: /usr/include/stdlib.h /usr/include/string.h
compute.o: /usr/include/sys/time.h /usr/include/sys/types.h
compute.o: /usr/include/sys/select.h linalg.h
compute.o: chull.h dual.h
dual.o: surf.h /usr/include/stdio.h /usr/include/math.h /usr/include/stdlib.h
dual.o: /usr/include/string.h /usr/include/sys/time.h
dual.o: /usr/include/sys/types.h 
dual.o: /usr/include/sys/select.h linalg.h dual.h chull.h
utils.o: surf.h /usr/include/stdio.h /usr/include/math.h
utils.o: /usr/include/stdlib.h /usr/include/string.h
utils.o: /usr/include/sys/time.h /usr/include/sys/types.h
utils.o: /usr/include/sys/select.h linalg.h
lp.o: surf.h /usr/include/stdio.h /usr/include/math.h /usr/include/stdlib.h
lp.o: /usr/include/string.h /usr/include/sys/time.h
lp.o: /usr/include/sys/types.h 
lp.o: /usr/include/sys/select.h linalg.h
chull.o: surf.h /usr/include/stdio.h /usr/include/math.h
chull.o: /usr/include/stdlib.h /usr/include/string.h
chull.o: /usr/include/sys/time.h /usr/include/sys/types.h
chull.o: /usr/include/sys/select.h linalg.h
chull.o: chull.h
tessel_cases.o: surf.h /usr/include/stdio.h /usr/include/math.h
tessel_cases.o: /usr/include/stdlib.h 
tessel_cases.o: /usr/include/string.h /usr/include/sys/time.h
tessel_cases.o: /usr/include/sys/types.h 
tessel_cases.o: /usr/include/sys/select.h linalg.h dual.h
tessel_patches.o: surf.h /usr/include/stdio.h /usr/include/math.h
tessel_patches.o: /usr/include/stdlib.h 
tessel_patches.o: /usr/include/string.h /usr/include/sys/time.h
tessel_patches.o: /usr/include/sys/types.h 
tessel_patches.o: /usr/include/sys/select.h linalg.h
tessel_convex.o: surf.h /usr/include/stdio.h /usr/include/math.h
tessel_convex.o: /usr/include/stdlib.h 
tessel_convex.o: /usr/include/string.h /usr/include/sys/time.h
tessel_convex.o: /usr/include/sys/types.h 
tessel_convex.o: /usr/include/sys/select.h linalg.h
tessel_concave.o: surf.h /usr/include/stdio.h /usr/include/math.h
tessel_concave.o: /usr/include/stdlib.h 
tessel_concave.o: /usr/include/string.h /usr/include/sys/time.h
tessel_concave.o: /usr/include/sys/types.h 
tessel_concave.o: /usr/include/sys/select.h linalg.h
tessel_torus.o: surf.h /usr/include/stdio.h /usr/include/math.h
tessel_torus.o: /usr/include/stdlib.h 
tessel_torus.o: /usr/include/string.h /usr/include/sys/time.h
tessel_torus.o: /usr/include/sys/types.h 
tessel_torus.o: /usr/include/sys/select.h linalg.h
  
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
  
