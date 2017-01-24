#ifndef _CMOIL_H
#define _CMOIL_H
//*************************************************************************************************
//*  Filename:   cmoil.hpp
//*
//*  Description: head file for Main class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 02 2005	Baohua Wang	Reshape
//*  Sep. 07 2007   Thomas Blom Update for OSX/Linux
//*************************************************************************************************

#include <fstream>
#include "trackball.h"
#include "cmoil_bspline.h"
#include "cmoil_globals.h"
#include "zvec.h"

#define INIT_TIMERS( count ) ZHashTable tHash; static int nTimed=10; static int nToTime=count
	// set nTimed to 0 if you want timing results printed out -- it prints them until nTimed is 5
#ifdef USE_TIMERS
#define START_TIMER( name ) ZTime name; if( nTimed<nToTime ) name=zTimeNow()
#define END_TIMER( name ) if( nTimed<nToTime ) { name=zTimeNow() - name; tHash.putD( #name, name ); }
#define TRACE_TIMERS() if( nTimed++ && 0 && nTimed<nToTime ) { \
							 ZTime accum = 0.0; \
							 for( ZHashWalk f( tHash ); f.next(); ) { \
								char   *key = f.key; \
								double *val = (double*)f.val; accum += *val;\
								trace( "  (%9.8lf secs) %s\n", *val, key ); \
							 }\
							 trace( "    total: %9.8lf\n", accum ); \
							 nTimed++; \
					   }
#else
#define START_TIMER( name )
#define END_TIMER( name )
#define TRACE_TIMERS() nTimed++
#endif


#define  EFFECT_AAptr  (NO_OVERLAP?CMOIL::Main::AAptrs[CMOIL::Main::ImageIndex]:CMOIL::Main::AAptrs[CMOIL::Main::Overlap1])

namespace CMOIL { 

extern tRGBA orange;
extern tRGBA purple;
extern tRGBA green;
extern tRGBA blue;
extern tRGBA red;
extern tRGBA cyan;
extern tRGBA yellow;
extern tRGBA white;
extern tRGBA ltGreen;
extern tRGBA ltYellow;
extern tRGBA zmoilColors[12];

// tRGBA zmoilColors[11] = { red, green, blue, yellow, orange, cyan, purple, ltRed, ltGreen, ltBlue, ltYellow, white };

#define tRGBA_to_Int( x ) ( x.r << 24 | x.g << 16 | x.b << 8 | 0xFF )
enum ZmoilColor {
	ZC_Red=0,
	ZC_Green,
	ZC_Blue,
	ZC_Yellow,
	ZC_Orange,
	ZC_Cyan,
	ZC_Purple,
	ZC_LtRed,
	ZC_LtGreen,
	ZC_LtBlue,
	ZC_LtYellow,
	ZC_White,
	ZC_NumColors,
};
				

bool initPath();

#define NO_OVERLAP     (CMOIL::Main::ImageIndex<CMOIL::Main::NumImages)

class Main : public State
{
public:
  static int MovieMode;

  static int Spinning;
  static int Moving;
  static int Rocking;
  static int Translating;
  static int BeginX, BeginY;
  static int Scaling;
  //static float Curquat[4];
  //static float Lastquat[4];
  static FQuat CurQuat;
  static FQuat LastQuat;
  static float ScaleFactor;
  static FVec3 Translate;
  static DVec3 BoundsMax;
  static DVec3 BoundsMin;
	// defines a bounding box that contains all atoms from all structures

  static float ax[3];
  static float angle;
  static int   ZclipFront;
  static int   ZclipBack;

  static void IconMain();

  static void recalcmodelView( AtomArray *model=0 );
  
  static void animate(void);
  static void animate_zclip(void);
  static void animate_rot();

  static void display_clip();
  static void display_draw(int windowSize);

  static void InitLightSource();
  static void SphereShiness();
  static void SetOverlapRMSDmsg();

  static void mouse(int button, int state, int x, int y);
  static void motion(int x, int y, float W, float H );


public:
  static int main(int argc, char **argv);
  static GLUquadricObj  *qobj;

public:
  static void display( int w=400, int h=400 );
  static void keyboard(unsigned char key, int x, int y);
  static void ZClip(int fron, int back);
  static void Rotates(float x, float y, float z);
  static int ReadFileParams(char *filename, char *filetype, char *cmoilexe);
  static void AddSphereToList();

  static int MainWindow;


};

}

#endif
