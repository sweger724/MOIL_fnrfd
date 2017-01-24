# directory where modules will be installed

bindir = $(prefix)/modules_@ARCH@

# names of module loadables

SUFFIXES = .in .c

bin_PROGRAMS = Surf

# module description files

bin_DATA = Surf.mdf

# source code for individual modules

Surf_SOURCES   = dxsurf.c Surf.mdf.in io.c compute.c \
		dual.c utils.c lp.c \
		 chull.c tessel_cases.c tessel_patches.c \
                 tessel_convex.c tessel_concave.c tessel_torus.c

BUILT_SOURCES = Surf.mdf.c

CLEANFILES = all.mdf Surf.mdf Surf.dll Surf.exp Surf.ilk Surf.lib Surf.pdb 
DISTCLEANFILES = ${BUILT_SOURCES}

if SGI_LINK
XTM_LINK = ld -n32 $(LDFLAGS) -o $@
else
XTM_LINK = $(CCLD) $(AM_FLAGS) $(CFLAGS) $(LDFLAGS) -o $@
endif

#CLEANFILES all.mdf Surf.mdf Surf.dll Surf.exp \
#           Surf.ilk Surf.lib Surf.pdb $(BUILT_SOURCES)
#
# Build Entry point functions for loadable modules
#

.in.c:
	awk -v dynamic=1 -v config=$(CONFIG_HEADER) -f ../../mdf2c.awk $* > $@

EXTRA_DIST = LICENSE

install-data-local:
	@modules=`ls $(prefix)/modules_@ARCH@/*.mdf` ;	\
	if test -r all.mdf ; then rm -f all.mdf ; fi ;		\
	for i in $$modules ; do					\
		if test $$i != "$(prefix)/modules_@ARCH@/all.mdf" ; then \
			echo "   Adding $$i to all.mdf" ;		\
			echo "" >> all.mdf ;			\
			cat $$i >> all.mdf ;			\
		else						\
			rm -f $$i ; 				\
		fi ;						\
	done ;							\
	echo " $(INSTALL_DATA) all.mdf $(DESTDIR)$(bindir)" ;	\
	$(INSTALL_DATA) all.mdf $(DESTDIR)$(bindir) ;		\
	echo "" ; echo "" ; echo "" ;				\
	echo "---------------------------------------------------" ; \
	echo "Module now installed. To use this loadable module," ; \
	echo "type the following when starting dx:" ;		\
	echo "" ;						\
	echo "dx -mdf $(prefix)/modules_@ARCH@/all.mdf" ;	\
	echo "" ;						\
	echo "It may useful to set this up as an alias." ; echo ""


clean-local:
	if test -r all.mdf ; then rm -f all.mdf ; fi
	rm -f *.mdf

