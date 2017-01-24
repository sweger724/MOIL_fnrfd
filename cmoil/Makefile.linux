#
# This is the Makefile for CMOIL on linux .
#
# Cmoil integrates GLUT, Magick++ library from ImageMagick, Mark Podlipec's gifmerge, 
# Amitabh Varshney's surf, and mpeg2encode from MPEG Software Simulation Group. 
#
# Before compiling CMOIL, please 
#
#  1. (you can skip it if GLUT has been installed)
#     install GLUT binary or build from cmoil by following:
#     1) download GLUT (http://www.opengl.org/resources/libraries/glut/glut-3.7.tar.gz) and
#     2) unpack it into $GLUT_HOME with "gzip -dc glut-3.7.tar.gz | tar xvf -" 
#
#  2. (you can skip it if ImageMagick has been installed)
#     install ImageMagick binary or build from cmoil by following:
#     1) download ImageMagick (ftp://ftp.imagemagick.org/pub/ImageMagick/ImageMagick-6.2.3-5.tar.gz) and
#     2) unpack it into $MAGICK_HOME with "gzip -dc ImageMagick-6.2.3-5.tar.gz | tar xvf -"
#
#  3. download MPEG (ftp://ftp.mpeg.org/pub/mpeg/mssg/mpeg2vidcodec_v12.tar.gz) and
#     unpack it into $MPEG_HOME with "gzip -dc mpeg2vidcodec_v12.tar.gz | tar xvf -" 
#
#  4. modify following if necessary:
# 	MAGICK_HOME GLUT_HOME MPEG_HOME CXX all 
#  
#  Now issue "make all" in this directory to build CMOIL and related executables 	
#
#  Notes: 
#    1. GLUT must be installed for compiling cmoil;  
#    2. Missing ImageMagick and/or MPEG packages will only effect some image saving functions;
#    3. surf and gifmerge packages come with the cmoil.
#
.PHONY:  cmoil cmoil_magick glut magick surf mpeg glut_install magick_install

#glut package
# NOTE: I have removed the glut header/lib from the cmoil distribution in preference
# of users having a stable/recent version installed on their machine.  (tfb)
GLUT_HOME=$HOME/thomas/freeglut-2.4.0/
#GL_LIBS =  -L/usr/lib -L/usr/X11R6/lib -lGL -lGLU -lglut -lm -lX11 -lXi -lXmu  -lXm
GL_LIBS =  -L /usr/local/lib -L/usr/lib -L/usr/X11R6/lib -lGL -lGLU -lglut -lm -lX11
GLUT_INCDIR = /usr/local/include
#-lXext  
        
# install ImageMagick into /usr from
# ftp://ftp.imagemagick.org/pub/ImageMagick/ImageMagick-6.2.3-5.tar.gz
MAGICK_HOME=/home/bwang/ImageMagick-6.2.3/
MAGICK_LIBS=  -lMagick++ -lMagick -lWand 

# MPEG coding
# ftp://ftp.mpeg.org/pub/mpeg/mssg/mpeg2vidcodec_v12.tar.gz 
MPEG_HOME=/home/bwang/mpeg2/
MPEG_DIR=$(MPEG_HOME)src/mpeg2enc/

#
#
# Use the following to force 32-bit
#CXX = g++ -m32 
CXX = g++
#
DEL = rm 
COPY= cp -f 
MOVE= mv -f 
LINK= g++
CD  = cd 
MKDIR = mkdir

CMOIL_DIR=./

SRC_DIR = $(CMOIL_DIR)
OBJ_DIR = $(CMOIL_DIR)obj/
OUT_DIR = $(CMOIL_DIR)../exe/

#######
CFLAGS = -I$(MAGICK_INCDIR) -I$(GLUT_INCDIR) -g
LFLAGS = $(GL_LIBS) $(MAGICK_LIBS)  

OBJS = $(OBJ_DIR)cmoil_bspline.o    $(OBJ_DIR)cmoil_atom.o     \
       $(OBJ_DIR)cmoil_pick.o       $(OBJ_DIR)cmoil_rms.o      \
       $(OBJ_DIR)cmoil_parm.o       $(OBJ_DIR)cmoil_util.o     \
       $(OBJ_DIR)cmoil_dib.o        $(OBJ_DIR)cmoil_msgboard.o \
       $(OBJ_DIR)trackball.o        \
       $(OBJ_DIR)cmoil_movieprint.o $(OBJ_DIR)gifmerge.o       \
       $(OBJ_DIR)cmoil_aminoacids.o $(OBJ_DIR)cmoil_ce.o       \
       $(OBJ_DIR)cmoil.o 	    $(OBJ_DIR)cmoil_menu.o     \
       $(OBJ_DIR)cmoil_drawshape.o

SRCS = $(OBJS: %.o %.cpp)

HDRS = cmoil.h cmoil_const.h cmoil_dib.h cmoil_movieprint.h   \
       cmoil_msgboard.h cmoil_parm.h cmoil_pick.h trackball.h

cmoil_all :  $(OUT_DIR) cmoil cmoil_magick 
all       :  glut cmoil_all magick mpeg surf cmoil_magick

$(OUT_DIR):
	-$(MKDIR) $(OUT_DIR)
	-$(MKDIR) $(OBJ_DIR)

cmoil : $(OBJS) $(HDRS)
	$(LINK) $(GL_LIBS) $(OBJS) -o $(OUT_DIR)$@  

cmoil_magick: $(OBJ_DIR)cmoil_magick.o
	$(LINK) $(MAGICK_LIBS) $^  -o $(OUT_DIR)$@

#compile surf
surf : 
	-@make -C $(SRC_DIR)surf
	-@$(COPY) $(SRC_DIR)surf/surf  $(OUT_DIR)

#Magick should be installed to default path
magick :
	-@cd $(MAGICK_HOME); ./configure
	-@make -C $(MAGICK_HOME)
	-@su - root -c "cd $(MAGICK_HOME); make install"

magick_install:
	-@su - root -c "cd $(MAGICK_HOME); make install"

# compile MPEG 
mpeg:
	-@make -C $(MPEG_DIR)
	-@$(COPY) $(MPEG_DIR)mpeg2encode  $(OUT_DIR)mpeg2enc
	-@$(COPY) $(SRC_DIR)mpeg2/mpeg2enc.par  $(OUT_DIR)

#compile glut
glut :
	@$(COPY) $(GLUT_HOME)linux/Glut.cf $(GLUT_HOME)
	@cd $(GLUT_HOME); mkmkfiles.imake
	@$(COPY) $(GLUT_HOME)linux/Makefile $(GLUT_HOME)lib/glut
	@make -C $(GLUT_HOME)lib/glut
	@echo Now you need root priviledge to install glut library to /usr/local/lib .........
	-@su - root -c "cd /usr/local/lib; $(COPY) $(GLUT_HOME)lib/glut/libglut.so.3.7 . ; ln -s libglut.so.3.7 libglut.so.3; ln -s libglut.so.3.7 libglut.so" 
glut_install: 
	-@su - root -c "cd /usr/local/lib; $(COPY) $(GLUT_HOME)lib/glut/libglut.so.3.7 . ; ln -s libglut.so.3.7 libglut.so.3; ln -s libglut.so.3.7 libglut.so" 
swap4 : $(OBJ_DIR)swap4.o
	$(LINK) $^ -o $@

swap8p : $(OBJ_DIR)swap8p.o
	$(LINK) $^ -o $@

clean :
	-$(DEL) $(OBJ_DIR)*.o

clobber: clean
	-$(DEL) $(OUT_DIR)cmoil


$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	$(CXX) $(CFLAGS) $< -c -o $@ 
