//*************************************************************************************************
//*  Filename:   cmoil.cpp
//*
//*  Description: 
//*    Main loop of graphic display in cmoil program.
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. ?? 2000	??		Initial Development
//*  Oct. 09 2000	Baohua Wang	Seperated header file cmoil.hpp
//*  Mar.    2001	Baohua Wang	Add text display
//*  Jul. 18 2001	Baohua Wang	Reshape to support mutiple structure display and overlap
//*  Jul. 28 2005	B. Wang 	Linux support
//*************************************************************************************************
#
#include "cmoil.h"
#include "cmoil_rms.h"
#include "cmoil_globals.h"
#include "cmoil_menu.h"
#include "cmoil_pick.h"
#include "cmoil_camera.h"
#include "cmoil_bspline.h"

#ifdef _WIN32
#include <windows.h>
#endif

#define  EFFECT_AAptr  (NO_OVERLAP?AAptrs[ImageIndex]:AAptrs[Overlap1])

namespace CMOIL {

/**************************************************************************************************/

// static variables of class State
OverlapAlign State::Ov;			// overlap
AtomArray *State::AAptrs[ARRAY_MAXNUM];	// All input Structures
PickedAtoms State::PickQ;	 // most recently selected 4 atoms
MessageBoard  State::MsgBrd;
MessageBoard  State::MsgBrd2;
Camera  State::camera;
char *State::ExeDir=NULL;          // set this up with calls to getenv()
int State::ImageIndex=0;	// index of current displayed/processed AtomArray, [0, NumImages), & NumImages(overlap)
int State::NumImages=0;	// total number of images for display
int State::Overlap1=0;	// remember 2 recent access structure indexes for overlap purpose
int State::Overlap2=0;
int State::PickMode=0;          // global Pick Mode (0/1)
int State::StereoMode=NONE_STEREO;
GLfloat State::Shine=0.0;       //0-1.0
bool State::ColorAnaglyph=false;
bool State::HasHardwareStereoSupport=false;  

// draw backbone/loop in either line or tube mode: start,stop,draw
bool  State::DrawBondAsLine=true;
bool  State::DrawBackboneAsCurve=true;

// static variables of class Main
int Main::MovieMode=0;
int Main::Spinning = 0;
int Main::Moving = 0;
int Main::Rocking = 0;
int Main::Translating = 0;
int Main::BeginX, Main::BeginY;
int Main::Scaling;
float Main::Curquat[4];
float Main::Lastquat[4];
float Main::ScaleFactor = 1.0;
float Main::TranslateX = 0.0;
float Main::TranslateY = 0.0;
double Main::MaxZ=0, Main::MinZ=0;   // object space
int Main::MainWindow=-1;


void initPath() {
	// Set the static State::ExeDir to the location of cmoil executables.
	// Added by Thomas Blom Sep2007 to fix problem when cmoil tries to 
	// launch cmoil_magick and doesn't know where to find it (because
	// popen apparently doesn't read the environment PATH?)
	 
	static char exePath[256];

	char *moilHome = getenv( "MOIL_HOME" );
	if( moilHome ) {
		strcpy( exePath, moilHome );
		int lastChar = strlen( exePath ) - 1;
		if( exePath[ lastChar  ] == '/' || exePath[ lastChar ] == '\\' ) {
			exePath[ lastChar ] = 0;
		}
		#ifdef WIN32
		strcat( exePath, "/moil.exe/" );	
		#else
		strcat( exePath, "/moil.source/exe/" );	
		#endif
		State::ExeDir = exePath;
	}
	else {
		printf( "CMOIL still relies on MOIL_HOME to be set to find the location\n");
		printf( "of cmoil_magick.  If this is not set, you might not be able to\n" );
		printf( "capture screen shots.  Please consider using ZMOIL instead, which\n"); 
		printf( "is the much-enhanced successor to CMOIL.\n\n"); 
	}
}

//**** Main Loop *************
int Main::main(int argc, char **argv)
{  
  AtomArray *aaptr=NULL;
  char      filename[FILENAME_MAXLEN], *filetype;
  int       topWindow;

  initPath();

  IconMain();

  glutInit(&argc, argv);
 
  NumImages = 0;
  for (int inFile=1; inFile<argc;)
  {
    if ( argv[inFile][0] == '-' )
      filetype = argv[inFile++];
    else
      filetype = "";
    strncpy(filename, argv[inFile++], FILENAME_MAXLEN);
    SetValidPath(filename);
    ReadFileParams(filename, filetype, argv[0]);
  }	
  ImageIndex = (NumImages>0)?0:-1;				// display first image by default 
  trackball(Curquat, 0.0, 0.0, 0.0, 0.0);

  DisplayModeInit();
  #ifdef WIN32
    // for some reason this call chokes on OSX and some linux variants
	glFinish();
  #endif

  topWindow=glutCreateWindow(CmoilAppWinName);
  MainWindow=topWindow;

  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);

  // Window init
  InitLightSource();

  //sphere list 
  for (int j=0; j<NumImages; j++)
    AddSphereToList();
  glFinish();
 
  MenuBar::mainMenu(topWindow);
  
  if ( NumImages>0 && AAptrs[ImageIndex]->InParm.structend > 0 )
    keyboard('m', 0, 0 );

  glutMainLoop();

  return 0;
}

#ifdef _WIN32
void Main::IconMain() {
  HWND w=GetForegroundWindow();
  if ( w != NULL)
    CloseWindow(w);
}
#else
void Main::IconMain() {;}
#endif

//---------------------------------------------------------------------------------------//
void Main::DisplayModeInit()
{
  HasHardwareStereoSupport=false;
  if ( StereoMode==HARDWARE_STEREO )
  {
    try {
       glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGB|GLUT_STEREO);
       HasHardwareStereoSupport=true;
    } 
    catch (char *str) {
       str=NULL;
       StereoMode=NONE_STEREO;
       fprintf(stderr, "CMOIL: unable to use hardware stereo!\n");
    } 
  }
  if (StereoMode != HARDWARE_STEREO )
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_ACCUM|GLUT_RGB);
}

void Main::lightposition(void)
{
  //GLfloat light_position[4] = {100.0,100.0,100.0,1.0};
  GLfloat light_position[4] = {1.0,1.0,1.0,0.0};

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glPopMatrix();
}

void Main::recalcmodelView(void)
{
  GLfloat m[4][4];
  build_rotmatrix(m, Curquat);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPushMatrix();
  glMultMatrixf(&m[0][0]);
  glTranslatef(TranslateX, TranslateY, 0.0);
  glScalef(ScaleFactor, ScaleFactor, ScaleFactor);
}

// the keyboard simulator to send out command as zCMDCMDCMDz
// assume zCMDCMDCMDz is atomic
//
// get command 
char* Main::keyboard_checkcmd(unsigned char key)
{
  static char cmd[250];
  static bool in_zcmd=false;
  if (!in_zcmd)
  {  
    cmd[0]=key;           //get a command
    cmd[1]='\0';     
    if (key=='z')
      in_zcmd=true;        // z-cmd starting
  } 
  else if (key == 'z')
   in_zcmd=false;            // z-cmd ending
  else
  {
    int cmdlen=strlen(cmd);
    if (cmdlen<249)
    {
      cmd[cmdlen++]=key;        // append a char (to z-cmd)
      cmd[cmdlen]='\0';
    } 
    else
      in_zcmd=false;
  }

  if (in_zcmd)
    return NULL;
  else
    return cmd;
}

//single char command
void Main::keyboard_chrcmd(char key) 
{
  AtomArray *aaptr;
  switch (key)
  {
      case 'p':                         // prepare for pick mode
        if ( PickMode ==0 )
        {
          PickMode = 1;
          if (NO_OVERLAP)
            PickQ.reset(ImageIndex);	// erase selected atom for different image
          else
            PickQ.reset(Overlap1, Overlap2);
        } 
        else
        {
          PickMode = 0;
          SetOverlapRMSDmsg();
        }
        break;           
      case 'm':                         // prepare for movie
        if (MovieMode == 1) {
          MovieMode = 0;
        } else {
          MovieMode = 1;
          glutIdleFunc(animate);
        }
        glutPostRedisplay();
        break;
      //        
      // other direct commands:
      case 'a' :                        // advance a frame, or next model of a PDB
        if ( PickMode == 1)
          break;
        MsgBrd.erase();
        int i;   //VC++ problem
        for (i=(NO_OVERLAP)?ImageIndex:Overlap1;;)
        {							// read at most two structures
          aaptr = AAptrs[i];
          if ( aaptr->InParm.filetag != fPDB )
          {
            aaptr->ReadNextStruct();
          } 
          else if (NO_OVERLAP)
          {
            aaptr->DisplayNextPdbModel();      
          }

          if ( NO_OVERLAP || i == Overlap2 )
            break;
          else 
            i = Overlap2;
        }

        glutPostRedisplay();
        break;
      case 'c':                    // advanced to next Chain
        if ( NO_OVERLAP )
        {   
          MsgBrd.erase();
          aaptr = AAptrs[ImageIndex];
          if ( aaptr!= NULL && aaptr->InParm.filetag == fPDB )
          {
            aaptr->DisplayNextPdbChain();
            glutPostRedisplay();
          }
        }
        break;
      case 's':      // save image into file
        aaptr = EFFECT_AAptr;
	  // use first structure's parameters for printing overlap picture
        if (!NO_OVERLAP)		
        {							// temp change bmp filename for AAptrs[0]
            aaptr->InParm.bmpName[strlen(aaptr->InParm.bmpName)-4] = '\0';
            strcat(aaptr->InParm.bmpName, "_ovlp.bmp");
        }
        aaptr->mprint_start(1^MovieMode);      // input: oneFrameOnly

        if (!NO_OVERLAP)
        {							// restore bmp filename
            aaptr->InParm.bmpName[strlen(aaptr->InParm.bmpName)-9] = '\0';
            strcat(aaptr->InParm.bmpName, ".bmp");
        }
        break;
      case 'r':      // rocking
        if ( Rocking > 0  )
        {
           Rocking=0;
           glutIdleFunc(NULL);
        }
        else
        {  
           MovieMode=0;
           Rocking=1; 
           glutIdleFunc(animate);
        }
        glutPostRedisplay();
        break;
  }  // case(key)
}  

//processing z-cmd
void Main::keyboard_zcmd(char *imagecmd) 
{
  AtomArray *aaptr;

  int dispBit=0;          
  int polymode=0;
  if (!strncmp(imagecmd, "zstick", 9))
          dispBit = STICK;
  else if (!strncmp(imagecmd, "zbackbone", 9))
          dispBit = BACKBONE;
  else if ( !strncmp(imagecmd, "zribbon1",  8) )
  {
          dispBit = RIBBON;
          polymode=GL_LINE;
  }
  else if ( !strncmp(imagecmd, "zribbon2",  8) )
  {
          dispBit = RIBBON;
          polymode=GL_FILL;
  }
  else if ( !strncmp(imagecmd, "zribbonwiden",  8) )
  {
         AAptrs[ImageIndex]->InParm.ribbonWidthFactor *= 1.1;      // increase 10%
         display();
  }
  else if ( !strncmp(imagecmd, "zribbonnarrow",  8) )
  {
         AAptrs[ImageIndex]->InParm.ribbonWidthFactor *= 0.9;      // decrease 10%
         display();

  }
  else if ( !strncmp(imagecmd, "zstickwiden",  8) )
  {
         AAptrs[ImageIndex]->InParm.stickWidthFactor += 1. ;      // increase 10%
         display();
  }
  else if ( !strncmp(imagecmd, "zsticknarrow",  8) )
  {
         AAptrs[ImageIndex]->InParm.stickWidthFactor -= 1.;      // decrease 10%
         if (AAptrs[ImageIndex]->InParm.stickWidthFactor < 1.0 )
           AAptrs[ImageIndex]->InParm.stickWidthFactor = 1.0; 
         display();
  }
  else if ( !strncmp(imagecmd, "zstructure",  8) )
  {
          dispBit = STRUCT2nd;
          polymode=GL_FILL;
  }
  else if ( !strncmp(imagecmd, "zsurf_",  6) )
  {
        dispBit = REDISPLAY;          // trun on surface first
        //polymode=GL_FILL;

        for (int i=(NO_OVERLAP)?ImageIndex:Overlap1;;)
        {							// read at most two structures
          aaptr = AAptrs[i];

          if ( !strncmp(imagecmd+6, "color", 5) )
          {
            if ( !strncmp(imagecmd+12, "res", 3) ) 
               aaptr->InParm.AAcolorGroups = ((aaptr->InParm.AAcolorGroups^SURF_COLOR_RES)&SURF_COLOR_RES);
            else
               aaptr->InParm.AAcolorGroups = ((aaptr->InParm.AAcolorGroups^SURF_COLOR_ATOM)&SURF_COLOR_ATOM);
          } 
          if ( !strncmp(imagecmd+6, "axe", 3) )
          {
            aaptr->setSurfCutoffPlane(1);
          } 
          else if ( !strncmp(imagecmd+6, "subarea", 4) )
          {
            aaptr->InParm.showSubSurf = !(aaptr->InParm.showSubSurf);  // highlight the area
          }
          else if ( !strncmp(imagecmd+6, "trans", 5) )
          {
            if ( aaptr->InParm.transparent == 255 )
              aaptr->InParm.transparent = 200;    // partial transparent
            else
              aaptr->InParm.transparent = 255;    
          }
          else if ( !strncmp(imagecmd+6, "cavity", 6) )
          {
             aaptr->InParm.showCavity = !(aaptr->InParm.showCavity);  
             if ( aaptr->InParm.showCavity )
               aaptr->InParm.cutoffval = 1.e+20;                    // no cutoff by default
          }
          else if ( !strncmp(imagecmd+6, "mesh", 4) )
          {
             aaptr->InParm.showSurfInMesh = !(aaptr->InParm.showSurfInMesh);  
          }
          else if ( !strncmp(imagecmd+6, "far", 3) )
          {
            aaptr->InParm.cutoffval -= aaptr->InParm.cutoffstep; 
          }
          else if ( !strncmp(imagecmd+6, "near", 4) )
          {
            aaptr->InParm.cutoffval += aaptr->InParm.cutoffstep; 
          }
          else if ( !strncmp(imagecmd+6, "nocut", 4) )
          {
            aaptr->InParm.cutoffval = 1.e+20; 
          }
          else if ( !strncmp(imagecmd+6, "pick", 4) )
          {
            aaptr->setSurfCutoffPlane(0);
          }
          else if ( !strncmp(imagecmd+6, "flip", 4) )
          {
            aaptr->flipSurfCutoffPlane();
          }
          else if ( !strncmp(imagecmd+6, "step", 4) )
          {
            aaptr->InParm.cutoffstep=atof(imagecmd+11);
            if (aaptr->InParm.cutoffstep<0)
              aaptr->InParm.cutoffstep = 0.5;
          }
          if ( NO_OVERLAP || i == Overlap2 )
            break;
          else 
            i = Overlap2;
        }
        glutPostRedisplay();
  }
  else if ( !strncmp(imagecmd, "zsurf",  5) )
  {
        dispBit = SURF;
        polymode=GL_FILL;
  }
  else if ( !strncmp(imagecmd, "zstereo",  7) )
  {
        if (!strncmp(imagecmd+8, "off",3) )
          StereoMode = NONE_STEREO;
        else if (!strncmp(imagecmd+8, "hardware",3) )
        {
          if (HasHardwareStereoSupport)
            StereoMode = HARDWARE_STEREO;
        }
        else
        {
          StereoMode = ANAGLYPH_STEREO;
          if ( !strncmp(imagecmd+8, "color", 5))
             ColorAnaglyph=true;
          else
             ColorAnaglyph=false;
        }
        glutPostRedisplay();
  }

  if ( dispBit != 0 )
  {
        for (int i=(NO_OVERLAP)?ImageIndex:Overlap1;;)
        {							// read at most two structures
          aaptr = AAptrs[i];
          if ( polymode != 0 ) 
            aaptr->InParm.ribbonpolymode=polymode;
          aaptr->InParm.dispmode ^= dispBit;   // switch on<->off :  0<->1
          if ( dispBit&RIBBON )                // solve ribbon/struct2nd conflict 
            aaptr->InParm.dispmode &= (0xffffffff^STRUCT2nd);
          else if ( dispBit&STRUCT2nd )
            aaptr->InParm.dispmode &= (0xffffffff^RIBBON);
          aaptr->ReProcessData();
          if ( NO_OVERLAP || i == Overlap2 )
            break;
          else 
            i = Overlap2;
        }
        display();
        //glutPostRedisplay();
  }
}

void Main::keyboard_pickcmd(char key)
{
    switch (key)
    {
      case 'd': 
          PickQ.print_distance();
          break;
      case 'a': 
          PickQ.print_angle();
          break;
      case 't': 
          PickQ.print_torsional();
          break;
      case '+':   		 // lighter
      case '=':
          PickQ.highlight('a');	//all
          break;
      case '-':    		// blacker
      case '_':
          PickQ.darken('a');	//all
          break;
      case 'r':     
          PickQ.highlight('r');	//red
          break;
      case 'b':    
          PickQ.highlight('b');	//blue
          break;
      case 'g':    
          PickQ.highlight('g');	//green
          break;
      case 'R':    
          PickQ.darken('r');	//red
          break;
      case 'B':    
          PickQ.darken('b');	//blue
          break;
      case 'G':    		
          PickQ.darken('g');	//green
          break;
      case '^':    		
          PickQ.reverse();		
          break;
    } //swith(key)
}

extern void spoofClick( int, int );

//process a key stroke 
void Main::keyboard(unsigned char key, int x, int y) 
{
	// debug
	//printf( "keyboard received %c\n", key);
	//if( key == 'd' ) {
	//	spoofClick( 5, 5 );
	//	return;
	//} 

  if (NumImages <= 0 )
    return;

  ////glFlush();
  if (MainWindow <0 )  return;
  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);
  
	
  char *cmd=keyboard_checkcmd(key); 
  if ( cmd != NULL)
  {
    // execute a keyboard command
    if ( *cmd=='z') 
      keyboard_zcmd(cmd);
    else
    {
      keyboard_chrcmd(*cmd);
      if (PickMode == 1) // pick mode 
        keyboard_pickcmd(*cmd);
    }
  }
}

void Main::mouse(int button, int state, int x, int y)
{
	// printf( "Mouse handler: button %d at (%d,%d), state=%d\n", button, x,y,state);
  ////glFlush();
  if ( NumImages <= 0 )
    return;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if (PickMode == 1) {
      PickQ.put(x, y);
    } else {
      Spinning = 0;
      Moving = 1;
      BeginX = x;
      BeginY = y;
    }
  }
  if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
    Moving = 0;
  }
  if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    Translating = 1;
    BeginX = x;
    BeginY = y;
  }

  if (button == GLUT_MIDDLE_BUTTON && state == GLUT_UP) {
    glutIdleFunc(NULL);
    Translating = 0;
  }
  if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
    Scaling = 1;
    BeginX = x;
    BeginY = y;
  }
  if (Spinning == 0 && Moving == 0) 
    glutIdleFunc(NULL);
  if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP) {
    Scaling = 0;
  }
}


void Main::animate(void)
{

  AtomArray *aaptr;

  if (NumImages <=0 )
    return;
  
  aaptr = EFFECT_AAptr;

  if (Rocking>0)
  {
    ax[0]=0;
    ax[1]=(Rocking==1)?1:-1;
    ax[2]=0;
    angle=aaptr->InParm.rotateDegree; 
    axis_to_quat(ax, angle, Lastquat);
    add_quats(Lastquat, Curquat, Curquat);
    Rocking=Rocking%2+1;
    SLEEP((int)300);
    glutPostRedisplay();
    return;
  }
  else
  {
    if (MovieMode) 
    {
      for (int k=(NO_OVERLAP)?ImageIndex:Overlap1;;)
      {
        aaptr=AAptrs[k];
       
        if ( aaptr->InParm.filetag == fPTH || aaptr->InParm.filetag == fDCD ) 
          aaptr->ReadNextStruct();     // only path and dcd support movie

        if (NO_OVERLAP)
          break;
        else if ( k != Overlap2)
          k = Overlap2;
        else
          break;
      }
      SLEEP((int)aaptr->InParm.sleepseconds);
    }
    if (Spinning == 1)
    {
      add_quats(Lastquat, Curquat, Curquat);
    }
  }
  //NewModel = 1;

  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);

  if ( aaptr->InParm.structend <= 0 )
  {
    glutPostRedisplay();
  }
  else
  {
    if (  (aaptr->InParm.filetag == fPTH || aaptr->InParm.filetag == fDCD ) &&
          aaptr->currentmol == aaptr->InParm.structend )
      glutPostRedisplay();
  }

  if (MovieMode)
  {
    aaptr = EFFECT_AAptr;
    aaptr->mprint_frame();
  }
  glFlush();
}

void Main::motion(int x, int y)
{

  ////glFlush();
  if (NumImages <= 0 )
     return;

  float W=glutGet(GLUT_WINDOW_WIDTH);
  float H=glutGet(GLUT_WINDOW_HEIGHT);
  if (Scaling) 
  {
    ScaleFactor = ScaleFactor * (1.0 + (((float) (BeginY - y)) / H));
    BeginX = x;
    BeginY = y;
    //NewModel = 1;
    glutPostRedisplay();
    return;
  }
  if (Moving) 
  {
    W *= 0.5;
    H *= 0.5;
    trackball(Lastquat,(BeginX-W)/W,(H-BeginY)/H,(x-W)/W,(H-y)/H);
    BeginX = x;
    BeginY = y;
    Spinning = 1;
    glutIdleFunc(animate);
  }

  if (Translating) {
    TranslateY -=   5.0*(y - BeginY)/float(H);
    TranslateX +=   5.0*(x - BeginX)/float(W);
    BeginX = x;  // save last location
    BeginY = y;  // save last location

    //NewModel = 1;
    glutPostRedisplay();
    return ;
  }
}

float Main::ax[3]; 
float Main::angle=0;

void Main::animate_rot()
{
  if (MainWindow < 0 || ImageIndex<0)  return;
  if (glutGetWindow() != MainWindow )  glutSetWindow(MainWindow);
  axis_to_quat(ax, angle, Lastquat);
  add_quats(Lastquat, Curquat, Curquat);
  SLEEP((int)300);
  glutPostRedisplay();
}

// (x,y,z)Axis= 0, 1, -1
void Main::Rotates(float xAxis, float yAxis, float zAxis) // deltaX=0, int deltaY=0, int deltaZ=0)
{
  if (MainWindow < 0 || ImageIndex<0)  return;
  if (glutGetWindow() != MainWindow )  glutSetWindow(MainWindow);
  
  AtomArray *aaptr=AAptrs[NO_OVERLAP?ImageIndex:Overlap1];
  //NewModel=1;
  ax[0]=-xAxis;
  ax[1]=-yAxis;
  ax[2]=-zAxis;
  if ( ax[0]==0 && ax[1]==0 && ax[2]==0)
    glutIdleFunc(NULL);
  else 
  {
    angle =  aaptr->InParm.rotateDegree; 
    axis_to_quat(ax, angle, Lastquat);
    add_quats(Lastquat, Curquat, Curquat);
    glutIdleFunc(animate_rot);
    glutPostRedisplay();
  }
}

int Main::ZclipFront=0, Main::ZclipBack=0;
void Main::animate_zclip()
{
  if (MainWindow < 0 || ImageIndex<0)  return;
  if (glutGetWindow() != MainWindow )  glutSetWindow(MainWindow);
  camera.Clip(ZclipFront, ZclipBack);
  SLEEP((int)300);
  glutPostRedisplay();
}

// front, back = -1,0,1,   -2:reset,  2: on|off
void Main::ZClip(int front, int back)
{
  if (front==back)
  {
    if ( front==ClipOnOff )
    { 
      if (ZclipFront!=NoClip || ZclipBack!=NoClip)
      {
        ZclipFront = NoClip;
        ZclipBack = NoClip;
      } else {
        ZclipFront=ClipOnOff;
        ZclipBack=ClipOnOff;
      }
    }
    else if (front==ResetClip)
    {
      camera.Clip(ResetClip, ResetClip);
      ZclipFront = NoClip;
      ZclipBack = NoClip;
    }
    glutIdleFunc(NULL);
  } 
  else
  {
    ZclipFront = front; 
    ZclipBack = back;
    camera.Clip(front, back);
    glutIdleFunc(animate_zclip);
  }
  if (MainWindow < 0 || ImageIndex<0)  return;
  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);
  glutPostRedisplay();
}

void Main::display_clip()
{
  // clip front and back of the screen picture
  if (ZclipFront!=0 || ZclipBack!=0 )
  {
    glMatrixMode(GL_MODELVIEW);
    GLdouble m[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, m);
    double z = MaxZ * m[10];     // convert maxz to eye coordinate
    double stepDelta = z*0.02; 
    if (camera.clipBack!=0 )  // back cut
    { 
      GLdouble eqn[4]={0,0,1,z-camera.clipBack*stepDelta}; 
      glClipPlane(GL_CLIP_PLANE1, eqn);
      glEnable(GL_CLIP_PLANE1);
    } 
    if (camera.clipFront>0 )  // front cut
    { 
      GLdouble eqn[4]={0,0,-1,z-camera.clipFront*stepDelta}; 
      glClipPlane(GL_CLIP_PLANE0, eqn );
      glEnable(GL_CLIP_PLANE0);
    } 
  }
}

//-------------------------------------------------------------------------
void Main::display_draw(int windowSize)
{
  AtomArray *aaptr;
  int colorFilter=(StereoMode!=ANAGLYPH_STEREO)?ColorFilterNone:((camera.eye==camera.LEFT)?ColorFilterRed:ColorFilterCyan);
 
  MsgBrd.print(1);		// 1:  bg = black, do not apply clip

  display_clip();
  Main::recalcmodelView();
  SphereShiness();

  Main::lightposition();
  
  int ii=0;
  int displayedStructs[8];
  int displayedStructCount=0;
  for ( int k=(NO_OVERLAP?ImageIndex:Overlap1); ; )
  {
    aaptr = AAptrs[k];
    if (aaptr == NULL)
      break;

    aaptr->DisplayStructure(ii, colorFilter, windowSize);
	displayedStructs[ displayedStructCount++ ] = aaptr->currentmol;

    if ((aaptr->dispmodeComb&BACKBONE)!=0 )          // global Backbone mode
      aaptr->backbone(windowSize, colorFilter);

    if ( (aaptr->dispmodeComb&STRUCT2nd) !=0 )
      aaptr->struct2nd(windowSize,colorFilter);
    else if ( (aaptr->dispmodeComb&RIBBON)!=0 ) 
      aaptr->ribbon(windowSize,colorFilter);

    if ( (aaptr->dispmodeComb&SURF) !=0 )
      aaptr->surf(colorFilter);

    // draw sphere net
    //if ( (aaptr->pickSphereCenter)>=0 && aaptr->InParm.showSphereNet )
    //{
    //  aaptr->SphereNet(colorFilter);
    //}

    if ( NO_OVERLAP || k == Overlap2 )
      break;				// only display one picture
    else 
      k = Overlap2;                 // display second structure

  } // for(k)


  if (PickMode==1  || (aaptr->InParm.dispmode&SURF )) 
    PickQ.display_pick(windowSize, colorFilter);

  glDisable(GL_CLIP_PLANE1);  // if any clip on, 
  glDisable(GL_CLIP_PLANE0);

  // Second MessageBoard added for lower-righth-hand display: currently always displays index of displayed struct
  char buffer[128];
  strcpy( buffer, "Index: ");
  for( int n=0; n<displayedStructCount; n++ ) {
	char buf[8];
	sprintf( buf, "%d", displayedStructs[ n ] );
	strcat( buffer, buf );
	if( n+1 < displayedStructCount ) {
		strcat( buffer, ", " );
	}
  }
  MsgBrd2.erase();
  if( displayedStructCount ) {
	MsgBrd2.set( buffer );
	MsgBrd2.print( 1, 1 );
		// print lower right corner
  }



}  // display_draw()


void Main::display(void)
{
  int wWidth, wHeight, windowSize;
 
  if (NumImages <= 0 )    // for empty command line filename
  {
    glClearColor(0.,0.,0.,1.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
    glutSwapBuffers();
    return;
  }
  if (MainWindow<0) return;
  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);

  glClearColor(AAptrs[0]->InParm.bg.r/255., 
               AAptrs[0]->InParm.bg.g/255., AAptrs[0]->InParm.bg.b/255., 1);

   /* Clear the buffers */
   glDrawBuffer(GL_BACK_LEFT);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   if ( HasHardwareStereoSupport)
   {
     glDrawBuffer(GL_BACK_RIGHT);
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
   }
   if ( StereoMode == ANAGLYPH_STEREO)
   { 
     glClearAccum(0.0,0.0,0.0,0.0);
     glClear(GL_ACCUM_BUFFER_BIT);
   }

  // classify window to small(1), middle(2), large(3)
  wWidth=glutGet(GLUT_WINDOW_WIDTH);
  wHeight=glutGet(GLUT_WINDOW_HEIGHT); 
  if ( wWidth<350 || wHeight<350)
    windowSize=1;    
  else if  ( wWidth<600 || wHeight<600) 
    windowSize=2;
  else
    windowSize=3;

  camera.SetSizes(wWidth, wHeight, MaxZ, MinZ); 

  camera.LookAt(camera.LEFT, StereoMode);

  display_draw(windowSize); 

  if ( StereoMode != NONE_STEREO)
  {
    if ( StereoMode == ANAGLYPH_STEREO )
    {
      glAccum(GL_LOAD, 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    }

    camera.LookAt(camera.RIGHT, StereoMode);
    display_draw(windowSize); 
  
    if ( StereoMode == ANAGLYPH_STEREO )
    {
      glAccum(GL_ACCUM, 1.0);
      glAccum(GL_RETURN,1.0);
    }
  } 
  if (StereoMode != HARDWARE_STEREO && HasHardwareStereoSupport )
  {
    glAccum(GL_LOAD, 1.0);
    glDrawBuffer(GL_BACK_RIGHT);
    glAccum(GL_RETURN,1.0);
  }
  glutSwapBuffers();
}  // display()

GLUquadricObj *Main::qobj=NULL;

void Main::ApplyBasicMatrices(double maxZ)
{
  glMatrixMode(GL_PROJECTION);
  //glPopMatrix();
  glLoadIdentity(); 

  //glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS, 128. /*85.0 */ );
  //glMaterialf(GL_BACK, GL_SHININESS, 35.0);
  Set_gluPerspective();
  glPushMatrix();

  glMatrixMode(GL_MODELVIEW);
  //glPopMatrix();
  glLoadIdentity();    
  gluLookAt(0.0, 0.0, maxZ+10,  /* Eye is at (0,0,maxZ) */
	      0.0, 0.0, 0.0,      /* Center is at (0,0,0) */
	      0.0, 1.0, 0.);      /* Up is in positive Y direction. */

  glPushMatrix();
}

/* Add to sphere List */
void Main::AddSphereToList() 
{
  if ( NumImages <= 0 )            /* check last image to add missing sphere to list */
    return;
  
  AtomArray *aaptr = AAptrs[NumImages-1];
  Atom *listi = aaptr->atomList;
  int i, l;
  Coordinate  crd;

  for (i=0;i<aaptr->numAtom;i++,listi++)
  {
    l = int(listi->radius*100+0.5);
    if (!glIsList(l) ) 
    {
      glNewList(l, GL_COMPILE);         // Create sphere display list.
	gluSphere(qobj, listi->radius, 20, 20);
	glEndList();
    }
  }
  if ( !glIsList(1) )              // add unit sphere, redius differ from Stick and backbone Tube
  {
    glNewList(1, GL_COMPILE); 
    gluSphere(qobj, /*Radius*/ Stick::StickRadius*3., /*Slices*/ 20,  /*Stacks*/ 20); 
    glEndList();
  }
  Tube::AddToQuadricList(qobj); 

  double maxZ = AAptrs[0]->maxZ, minZ=AAptrs[0]->minZ ;
  for (i=1; i<NumImages; i++)
  {
    if (maxZ<AAptrs[i]->maxZ)
      maxZ=AAptrs[i]->maxZ;
    if (minZ>AAptrs[i]->minZ)
      minZ=AAptrs[i]->minZ;
  }

  // Handle ModelView and Projection Matrixes 
  ApplyBasicMatrices(maxZ);
  MaxZ = maxZ;
  MinZ = minZ;
}

void Main::SphereShiness()
{
  // specular color
  GLfloat sp[4]={Shine, Shine, Shine, 1.0}; 
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,sp);
  glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,100.);
}

void Main::InitLightSource () 
{
  qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj, GLU_FILL);

  // diffuse color 0.75
  GLfloat light_diffuse[4]  = {0.75, 0.75, 0.75, 1.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);

  lightposition(); 
  glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 5.0*ScaleFactor);

  GLfloat sp[4]={1.0, 1.0, 1.0, 1.0};
  glLightfv(GL_LIGHT0, GL_SPECULAR, sp);

  // ambient color
  GLfloat tv[4] = {0.5,0.5,0.5,1.0};    // grey color 
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,tv);   
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,0.0);

  SphereShiness(); 

  //glMaterialf(GL_FRONT,GL_SHININESS,50); 
  //glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,100.); 

  GLfloat clr[]={0.1, 0.5, 0.8, 1.0};
  //glMaterialfv(GL_BACK,GL_AMBIENT_AND_DIFFUSE,clr);
  glMaterialfv(GL_BACK,GL_AMBIENT_AND_DIFFUSE,clr);

  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_DEPTH_TEST);
}


// handle mutiple Molecular display,  value = imageOrder = imageIndex+1
void Main::SetOverlapRMSDmsg()
{        
  if (NO_OVERLAP)
    MsgBrd.erase();
  else 
  {
    char msg[120];
    switch (Ov.ovtype)
    {
      case 'A': 
          sprintf(msg, "%s,%s aligned with CE: length=%d, Rmsd=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.autoAlign.numAlignedMonos, Ov.rmsd);
          break;

      case 'B': 
          sprintf(msg, "%s,%s aligned by backbone: %d residues, %d atoms, RMSD=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.autoAlign.numAlignedMonos, Ov.numAlignedAtoms, Ov.rmsd);
          break;
      case 'C': 
          sprintf(msg, "%s,%s aligned by CA: %d residues, RMSD=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.autoAlign.numAlignedMonos, Ov.rmsd);
          break;
      case 'D':
          sprintf(msg, "%s,%s aligned by selected particles: %d atoms, RMSD=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.numAlignedAtoms, Ov.rmsd);
          break;
      case 'E':
          sprintf(msg, "%s,%s aligned by selected backbone: %d atoms, RMSD=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.numAlignedAtoms, Ov.rmsd);
          break;
      case 'F':
          sprintf(msg, "%s,%s aligned by selected CA: %d atoms, RMSD=%g", 
                       Ov.structs[0], Ov.structs[1], Ov.numAlignedAtoms, Ov.rmsd);
          break;
    }
    MsgBrd.set(msg);
  }
}

void Main::showMolImg(int value, char overlapType)
{
  if ( MainWindow<0 )  return;
  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);

  for (int i=0; i<4; i++)
    Curquat[i]=0.0;
  trackball(Curquat, 0.0,0.0,0.0,0.0);
  TranslateX=TranslateY=0.;
  //NewModel=1;          // reset model view

  if ( value > 0 && value <= NumImages )
  {
    if (NO_OVERLAP)
      MsgBrd.erase();           // keep the possible save crd message

    ImageIndex = value-1;       // set img to display
    if (AAptrs[ImageIndex]->orgCrdModified)
      AAptrs[ImageIndex]->ReProcessData();     
  
    if (Overlap2 != ImageIndex)
    {
      Overlap1 = Overlap2;      // turn around the img index
      Overlap2 = ImageIndex;    
    }
  }
  else if ( (value==NumImages+1) && (overlapType!= Ov.ovtype || ImageIndex!=NumImages ) )   // ignore if already overlapped
  {
    Ov.ovtype=overlapType;
    if (Overlap1 != Overlap2)
    {
      if (AAptrs[Overlap1]->orgCrdModified)
        AAptrs[Overlap1]->ReProcessData();     
      if (AAptrs[Overlap2]->orgCrdModified)
        AAptrs[Overlap2]->ReProcessData();     

      if ( Ov.overlap(*(AAptrs[Overlap1]), *(AAptrs[Overlap2]) ) >0 )
      {
        ImageIndex = NumImages;			// this is a overlap condition
        SetOverlapRMSDmsg();
      }
      else
      {
        PickMode = 0;
        MsgBrd.set(Ov.ErrMsg);
      }
    }
  }
 
  if ( PickMode == 1)
  {
    if ( NO_OVERLAP )
      PickQ.reset(ImageIndex);
    else
      PickQ.reset(Overlap1, Overlap2);
  }
  display();
}

// set globals: ImageIndex  NumImages AAptrs[ImageIndex] 
void Main::ReadFileParams(char *fileName, char *fileType, char *cmoiexe)
{
  AtomArray *aaptr=NULL;
  int        currCrd=0;

  if (MainWindow>=0)
  if ( glutGetWindow() != MainWindow ) glutSetWindow(MainWindow);

  if (fileName == NULL)
    return;

  aaptr = new AtomArray ;
  if (aaptr == NULL)
  {
    fprintf(stderr, "Error: Unable to allocate memory allocation at ReadAnInstrFile().\n");
    exit(-1);
  }

  if ( strcmp(fileType, "-pdb" ) == 0 )
  {
    strcpy(aaptr->InParm.crdName, fileName);      // no instruction file
    aaptr->InParm.filetag=fPDB;
  }
  else if ( strcmp(fileType, "-xyz") == 0 )
  {
    strcpy(aaptr->InParm.crdName, fileName);   // no instruction file
    aaptr->InParm.filetag=fXYZ;
  }
  else
    currCrd = aaptr->GetInputParameters(fileName, currCrd, cmoiexe);

  if ( currCrd >= 0)	
  {
    ImageIndex = NumImages;
    AAptrs[ImageIndex] = aaptr;  
    NumImages++; 
    aaptr->ProcessData();
  }
  else
    delete aaptr;
}

}


int main(int argc, char **argv)
{
  return CMOIL::Main::main(argc, argv);
}


