#ifndef _CMOIL_CONST_H
#define _CMOIL_CONST_H
//****************************************************************************************************
//*  Filename:   cmoil_const.h
//*
//*  Description: constant definitions for cmoil
//*
//*  History:
//*  Date		Developer	Description
//*  ------------ ----------- ------------------------------------------------------------------------
//*  June 2001	B. Wang	Initial Development
//*  Sept 2007  Thomas Blom   Update for OSX/Linux
//****************************************************************************************************
#ifdef __APPLE__
	#include <GLUT/glut.h>
#elif defined(WIN32)
	#include "GLwin32/glut.h"
#else
	#include <GL/glut.h>
#endif
#include <math.h>

namespace CMOIL {
enum DisplayMode 				// for both overall & single atom's display mode
{
  SKIP=0, SPACEBALL=0x1, STICK=0x2, RIBBON=0x4, BACKBONE=0x8, UNITBALL=0x10, 
  STICKBALL=0x20, STRUCT2nd=0x40, SURF=0x80, SUBSURF=0x100, 
  INSPHERE=0x200,                   // in sphere selection, temporary use

  REDISPLAY=0x80000000,

  NONE_STEREO = 0x10000,           // for StereoMode
  HARDWARE_STEREO = NONE_STEREO+1,  //  hardware stereo is available
  ANAGLYPH_STEREO = NONE_STEREO+2
}; 	

typedef char SecStructType;

enum SecStructTypeList               // list of secondary structure
{
  STRUCT_INIT=0x80,                  // no structure exists
  STRUCT_UNKNOWN=0,
  STRUCT_PHI=0x1,                    // phi<=[-170, 0] for helix, other for [0,190]
  STRUCT_PSI=0x2,                    // psi<=0
  STRUCT_END=0x4,                    // end of beta sheet
  STRUCT_SHEET=STRUCT_PHI,
  STRUCT_LASTOFSHEET=STRUCT_SHEET|STRUCT_END,   // last residure of a sheet structure
  STRUCT_HELIX=STRUCT_PHI|STRUCT_PSI,
  STRUCT_CHAINEND=0x10,              // end of a chain
  STRUCT_OFOUND=0x20                 // found O
};

#define sIsHelix(xx)    (((xx)&0xF)==STRUCT_HELIX)
#define sIsSheet(xx)    (((xx)&0xF)==STRUCT_SHEET)
#define sIsUnknown(xx)  (((xx)&STRUCT_PHI)==0)
#define sIsChainEnd(xx) (((xx)&STRUCT_CHAINEND)==STRUCT_CHAINEND)
#define sFoundO(xx)     (((xx)&STRUCT_OFOUND)==STRUCT_OFOUND)
#define sFoundHBond(xx) (((xx)&STRUCT_4I_HBOND)==STRUCT_4I_HBOND)
#define sLowType(xx)     ((xx)&0x0F)
   									
const unsigned int FILENAME_MAXLEN  =	256;		// max length for a filename
const unsigned int LINE_MAXLEN  	=	300;		// max chars per line
const unsigned int MONO_MAXNUM  	= 	20000;	// max number of monomers(residuels)
const unsigned int ATOM_MAXNUM  	=   	100000;	// max number of atoms,for ReadXYZ only
const unsigned short int MAXLEN_NAMES  = 	5;		// max chars for an atom or a residue
typedef char NameString[MAXLEN_NAMES];			// atom or residue name

const GLubyte DELTA_COLOR	=	30; 	      // for rgb 

const unsigned short int ATOM_MAX_PICK 	  =	4;	// max number of picked atoms in select mode
const unsigned short int ATOM_MAX_NEIGHBORS =	14;	// for atom neighbor counting, consider alternatives

const unsigned short int ARRAY_MAXNUM	  =	20;	// maximum number of input pictures

const double PI = 3.14159265359;					// Pi

const float PI10th=PI/10.;
const float PI2DEGREE=180.0/PI;

// true constants
//
const int TubeListStart=100000;
const int NoSelectLoadName=2000000;                   // Name not for GL_SELECT
const double AngleToRadius=PI/180;
const int BsplineInterPoints = 5 ;      
const char CmoilAppWinName[]="Molecular Viewer 1.0";

//#define Color2Mono(r,g,b)  ((GLubyte)(0.30*(r)+0.59*(g)+0.11*(b)))
#define Color2Mono(r,g,b)  (0.5*(r)+0.3*(g)+0.2*(b))
class tRGBA                                 //RGBA color  0-255
{
public:
  GLubyte r;
  GLubyte g;
  GLubyte b;
  GLubyte w;       // grey scale

  tRGBA(void) { set(0,0,0); };
  tRGBA(GLubyte rm,  GLubyte gm, GLubyte bm )  { set(rm,gm,bm); };
  void set(GLubyte rm, GLubyte gm, GLubyte bm) { r=rm;g=gm;b=bm;setGrey(); };
  void setGrey(void) { w=Color2Mono(r,g,b); } ;
  friend tRGBA operator*(const float &f, const tRGBA &clr) { tRGBA c(int(clr.r*f), int(clr.g*f), int(clr.b*f)); return c; };
};
// cp : pointer to tRGBA structure

typedef unsigned int tLEN;

typedef enum {SURF_COLOR_UNI=0, SURF_COLOR_RES=0x1, SURF_COLOR_ATOM=0x2 } tSurfColorType;    

enum { NoClip=0, MinusClipDelta=NoClip+1, SameClipDelta=NoClip+2, PlusClipDelta=NoClip+3, ResetClip=NoClip-1, ClipOnOff=NoClip-2 };
enum { ColorFilterNone=0, ColorFilterRed=1, ColorFilterCyan=2};

enum {SurfCavity=0, SurfOuterSelf=1, SurfOuterNeighborCheck=2}; 
enum {X=0, Y=1, Z=2 };       // crd index in an array


//*************************************************************
// anaglyph colors
//
// ColorAnaglyph;  // effected only when filter!=ColorFilterNone

// for mono color anaglyph
#define PureAnaglyphR(filter,cp)   	((filter)==ColorFilterCyan?0:(cp)->w)   
#define PureAnaglyphG(filter,cp)   	((filter)==ColorFilterCyan?(cp)->w:0) 
#define PureAnaglyphB(filter,cp)   	((filter)==ColorFilterCyan?(cp)->w:0)  

// for color ganglyph
#define ColorAnaglyphR(filter,cp)   ((filter)==ColorFilterCyan?0:((cp)->w+(cp)->r)*0.5)  
#define ColorAnaglyphG(filter,cp)   ((filter)==ColorFilterCyan?  ((cp)->w+(cp)->g)*0.5:0)
#define ColorAnaglyphB(filter,cp)   ((filter)==ColorFilterCyan?  ((cp)->w+(cp)->b)*0.5:0) 

// cp : pointer to tRGBA
#define Color2R(filter,cp)  ((unsigned char)((filter)==ColorFilterNone?(cp)->r:((ColorAnaglyph)?ColorAnaglyphR(filter,cp):PureAnaglyphR(filter,cp))))
#define Color2G(filter,cp)  ((unsigned char)((filter)==ColorFilterNone?(cp)->g:((ColorAnaglyph)?ColorAnaglyphG(filter,cp):PureAnaglyphG(filter,cp))))
#define Color2B(filter,cp)  ((unsigned char)((filter)==ColorFilterNone?(cp)->b:((ColorAnaglyph)?ColorAnaglyphB(filter,cp):PureAnaglyphB(filter,cp))))
//*************************************************************

#ifdef _WIN32
#include <direct.h>
#include <process.h>
#define UNLINK	_unlink
#define POPEN	_popen
#define PCLOSE	_pclose
#define MKDIR	_mkdir
#define RMDIR	_rmdir
#define CHDIR	_chdir
#define GETCWD	_getcwd
#define SLEEP(x)	Sleep(x)
#define STAT    _stat
#define DIRMASK	_S_IFDIR
#define DIR1CMD  "DIR /W /O:GN /B"
#define DIRDELIM '\\'
#define DIRDELIMalt '/'
#define BINSURFFIX ".exe"
#define WRAP_CE    "ce.bat"
#define WRAP_EXE   "ce.bat"
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/select.h>
#define UNLINK	unlink
#define POPEN	popen
#define PCLOSE	pclose
#define MKDIR(pname)  mkdir(pname,S_IRWXU)
#define RMDIR	rmdir
#define CHDIR	chdir
#define GETCWD	getcwd
//#define SLEEP(x)	sleep(x/100)     x is in 1/1000 second
#define SLEEP(x) { struct timeval tv;          \
                   tv.tv_sec=0; tv.tv_usec=(x)*1000;  \
                   select ( 0, NULL, NULL, NULL, &tv); }
#define STAT    stat
#ifdef __APPLE__
	#define DIRMASK	S_IFDIR
#else
	#define DIRMASK	__S_IFDIR
#endif

#define DIR1CMD  "ls -1"
#define DIRDELIM '/'
#define DIRDELIMalt '\\'
#define BINSURFFIX " "
#define WRAP_CE    "ce.sh"
#define WRAP_EXE   "exec.sh"
#endif
}

#endif

