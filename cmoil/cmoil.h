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
#ifdef __APPLE__
	#include <GLUT/glut.h>
#elif defined(WIN32)
	#include "GLwin32/glut.h"
	//#include <GL/glu.h>
#else 
	#include <GL/glut.h>
#endif
#include <fstream>

#include "trackball.h"
#include "cmoil_bspline.h"
#include "cmoil_globals.h"

namespace CMOIL { 

#define NO_OVERLAP     (ImageIndex<NumImages)

class Main : public State
{
private:
  static int MovieMode;

  static int Spinning;
  static int Moving;
  static int Rocking;
  static int Translating;
  static int BeginX, BeginY;
  static int Scaling;
  static float Curquat[4];
  static float Lastquat[4];
  static float ScaleFactor;
  static float TranslateX;
  static float TranslateY;
  static double MaxZ, MinZ;   // object space

  static float ax[3];
  static float angle;
  static int   ZclipFront;
  static int   ZclipBack;

  static void IconMain();

  static void DisplayModeInit();
  static void recalcmodelView();
  
  static char* keyboard_checkcmd(unsigned char key);
  static void  keyboard_chrcmd(char key);
  static void  keyboard_pickcmd(char key);
  static void  keyboard_zcmd(char *cmd);

  static void animate(void);
  static void animate_zclip(void);
  static void animate_rot();

  static void display_clip();
  static void display_draw(int windowSize);

  static void ApplyBasicMatrices(double maxZ);
  static void InitLightSource();
  static void SphereShiness();
  static void SetOverlapRMSDmsg();

  static void mouse(int button, int state, int x, int y);
  static void motion(int x, int y);


public:
  static int main(int argc, char **argv);
  static GLUquadricObj  *qobj;
  static void lightposition();

protected:
  static void display(void);
  static void keyboard(unsigned char key, int x, int y);
  static void ZClip(int fron, int back);
  static void Rotates(float x, float y, float z);
  static void ReadFileParams(char *filename, char *filetype, char *cmoilexe);
  static void showMolImg(int value, char overlapType);
  static void AddSphereToList();

  static int MainWindow;


};

}

#endif
