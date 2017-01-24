//*************************************************************************************************
//*  Filename:   cmoil_menu.cpp
//*
//*  Description: 
//*    Main menu for cmoil.  Export function:  CmoilMainMenu(int topWindow)
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Jun. 25 2003	Baohua Wang	Initial Development
//*  Sep. 07 2007   Thomas Blom Update for OSX/Linux
//*************************************************************************************************
// to manager Menus
#include "cmoil_menu.h"
#include "cmoil_camera.h"

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

#ifdef _WIN32
	#include <direct.h>
	#include <windows.h>
#endif

#if !defined( __APPLE__ ) && !defined( WIN32 )
	#include <GL/glx.h>
#endif

#ifndef WIN32	
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#endif

namespace CMOIL {


#ifndef _WIN32
int g_windowsFound=0;
Window FindWindow(Display *display, Window pwin, const char winname[])
{
  char *wname=NULL;
  int fetchStatus = XFetchName(display, pwin, &wname);
  if( fetchStatus ) {
  	//printf( "FindWindow called with parent = %s", wname );
  	XFree( wname );
  }	
  else {
  	//printf( "FindWindow called with pwin = %d, but no name found.", pwin );
  }
    
  Window r, p, *kids; 
  unsigned int numkids;
  if ( XQueryTree(display, pwin, &r, &p, &kids, &numkids) == 0 ) {
  	 printf( "XQueryTree Failed.  Returning. \n" );
     return 0;	// error occured
  }
  else {
  	//printf( "   ( %d children found )\n", numkids );
  }
 
  XWindowAttributes attr;
  int foundID = 0;
  for (int i=0; foundID==0 && i<numkids; i++)
  {
    XGetWindowAttributes(display, kids[i], &attr);
    g_windowsFound++;
    //if ( attr.map_state!=IsViewable )
    //  continue;
    if ( XFetchName(display, kids[i], &wname) != 0 && strcmp(winname, wname)==0 )   
    {
       foundID = kids[i];
       printf("(%d) +win=%s (%d): (%d,%d) %d %d \n", g_windowsFound, wname, foundID, attr.x, attr.y, attr.width, attr.height);
       XFree(wname);
       continue;
    }
    //printf("(%d) -win=%s: (%d,%d) %d %d \n", g_windowsFound, wname, attr.x, attr.y, attr.width, attr.height);
    if ( wname != NULL )
      XFree(wname);

    foundID = FindWindow( display, kids[i], winname );
  }
  return foundID;
}
#endif

//*********************************************************
//*  MenuButton
//*********************************************************
void MenuButton::set(int mywindowId, char* mytext, int mygroup, int mynameId, int pWindow)
{
  windowId=mywindowId;
  text=mytext;
  group=mygroup;
  nameId=mynameId;
  parentId=pWindow;
}

void MenuButton::draw(const char *dstr)
{
   glPushMatrix();
   gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );   // (left, right, bottom, top)
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
   glBegin(GL_POLYGON);
   if ( *dstr == '|' || *dstr == '>' )
     glColor3f(0,0,0);
   else
     glColor3f(0.8,0.8,0.8);
   glVertex3f(0,0,-1);
   glVertex3f(0,1,-1);
   glVertex3f(1,1,-1);
   glVertex3f(1,0,-1);
   glEnd();

   // button style
   if ( !strncmp(dstr," X",2) || !strncmp(dstr, " Y",2)   || !strncmp(dstr, " Z",2) ||
        !strncmp(dstr,"Front",5) || !strncmp(dstr,"Back",4)  )
   {             // button shape
     glLineWidth(2.0);
     glBegin(GL_LINE_STRIP);
     glColor3f(0.5,0.5,0.5);
     glVertex2f(0,0);
     glVertex2f(0,0.9);
     glVertex2f(0.9,0.9);
     glVertex2f(0.9,0);
     glVertex2f(0,0);
     glEnd();
   }
   else if ( strncmp(dstr,"Rot",3) && strncmp(dstr, "Clip",4) && *dstr != '>' && *dstr != '|' )
   {             // button shape
     glLineWidth(2.0);
     glBegin(GL_LINE_STRIP);
     glColor3f(2.0,2.0,2.0);     // light lines
     glVertex2f(0,0);
     glVertex2f(0,1);
     glVertex2f(1,1);
     glEnd();
     glLineWidth(3.0);
     glBegin(GL_LINE_STRIP);     // dark lines
     glColor3f(0.5,0.5,0.5);
     glVertex2f(1,1);
     glVertex2f(1,0);
     glVertex2f(0,0);
     glEnd();
   }
   glColor3f ( 0.0, 0.0, 0.0 );		      // text color
   glRasterPos2f(0.05,0.25);
     for ( int i=0; i< (int)strlen(dstr); i++ )
       glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, dstr[i]);
   // draw bitmap for rotations
   // glBitmap(24,13,0,0,0,0, (GLubyte*)xRotate);
  
   glPopMatrix();
   glutSwapBuffers();
}

void MenuButton::display()
{
  glutSetWindow(windowId);
  if (group == MB_NOTSUB || (group==MB_PICK && PickMode>0) )
  {
    draw(text);
    if (group==MB_PICK)
      glutShowWindow();
  } 
  else if ( group==MB_PICK )
    glutHideWindow();
}


//*********************************************************
//*    Menu Bar
//*********************************************************

// init static members
int MenuBar::windowId=-1;
int MenuBar::parentId=-1;
int MenuBar::X=0;
int MenuBar::optionX=0;
int MenuBar::Y=0;
int MenuBar::nButtons=0;
MenuButton MenuBar::buttons[MAX_NUM_MBUTTONS];
int MenuBar::button2_menuid=-1;
int MenuBar::button2_pulldownOverlap=-1;
int MenuBar::filelist_menuid=-1;
char * MenuBar::openfileType="";

// callback functions
void MenuBar::displayNop()    // for function in CreateMenuButton()
{
  ;
}

void MenuBar::RotateX(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
  {
    if (button == GLUT_LEFT_BUTTON)
      Rotates(-1,0,0);
    else
      Rotates(1,0,0);
  }
  else 
    Rotates(0,0,0);
}
void MenuBar::RotateY(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
  {
    if (button == GLUT_LEFT_BUTTON)
      Rotates(0,-1,0);
    else
      Rotates(0,1,0);
  }
  else 
    Rotates(0,0,0);
}

void MenuBar::RotateZ(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
  {
    if (button == GLUT_LEFT_BUTTON)
      Rotates(0,0,-1);
    else
      Rotates(0,0,1);
  }
  else 
    Rotates(0,0,0);
}

void MenuBar::ClipFront(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
  {
    if (button == GLUT_LEFT_BUTTON)
      ZClip(PlusClipDelta,SameClipDelta);
    else
      ZClip(MinusClipDelta,SameClipDelta);
  }
  else
    ZClip(SameClipDelta,SameClipDelta);
}

void MenuBar::ClipBack(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
  {
    if (button == GLUT_LEFT_BUTTON)
      ZClip(SameClipDelta,PlusClipDelta);
    else
      ZClip(SameClipDelta,MinusClipDelta);
  } 
  else
    ZClip(SameClipDelta,SameClipDelta);
}

void MenuBar::ClipNoClip(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN)
    ZClip(ClipOnOff,ClipOnOff);   // on/off
}

// menu commands split Menu1 into 2 types: Menu1_1, Menu2_2
//
void MenuBar::Menu1_1(int value )      // menu function part1 for MenuButton#1
{
  char cmd[50]="";
  unsigned int i;

  switch (value)
  { 
  case M1_BACKBONE:
    strcpy(cmd, "zbackbonez");
    break;
  case M1_SURF:				  //  surf menu entries 
    strcpy(cmd, "zsurfz");
    break;
  case M1_SURF_COLOR:
    strcpy(cmd, "zsurf_color_resz");
    break;
  case M1_SURF_COLOR_ATOM:
    strcpy(cmd, "zsurf_color_atomz");
    break;
  case M1_SURF_FAR:
    strcpy(cmd, "zsurf_farz");       // submenu
    break;
  case M1_SURF_NEAR:
    strcpy(cmd, "zsurf_nearz");       // submenu
    break;
  case M1_SURF_NOCUT:
    strcpy(cmd, "zsurf_nocutz");       // submenu
    break;
  case M1_SURF_AREA:
    strcpy(cmd, "zsurf_subareaz");       // submenu
    break;
  case M1_SURF_TRANS:
    strcpy(cmd, "zsurf_transparentz");       // submenu
    break;
  case M1_SURF_CAVITY_OFF:
    strcpy(cmd, "zsurf_cavityz");       // submenu
    break;
  case M1_SURF_MESH :
    strcpy(cmd, "zsurf_meshz");
    break;
  case M1_SURF_SWITCH_AXE:
    strcpy(cmd, "zsurf_axez");       // submenu
    break;
  case M1_SURF_PICK:
    strcpy(cmd, "zsurf_pickz");       // submenu
    break;
  case M1_SURF_FLIP:
    strcpy(cmd, "zsurf_flipz");       // submenu
    break;
  case M1_STICK:
    strcpy(cmd, "zstickz");
    break;
  case M1_STICK_WIDEN:                   // share width parameter b/w ribbon & 2nd_struct
    strcpy(cmd, "zstickwidenz");
    break;
  case M1_STICK_NARROW:
    strcpy(cmd, "zsticknarrowz");
    break;
  case M1_WHITEBG:
    AAptrs[0]->InParm.bg.r ^= 0xFF;
    AAptrs[0]->InParm.bg.g ^= 0xFF;
    AAptrs[0]->InParm.bg.b ^= 0xFF;
    display();
    break;
  case M1_SHINE_NO:
    Main::Shine=0.0;
    display();
    break;
  case M1_SHINE_LOW:
    Main::Shine=0.1;
    display();
    break;
  case M1_SHINE_MEDIUM:
    Main::Shine=0.2;
    display();
    break;
  case M1_SHINE_HIGH:
    Main::Shine=0.4;
    display();
    break;
  case M1_BACKBONE_QUALITY:
    DrawBackboneAsCurve = !DrawBackboneAsCurve; 
    display();
    break;
  case M1_STICK_QUALITY:
    DrawBondAsLine = !DrawBondAsLine;
    display();
    break;
  case M1_QUALITY:
    DrawBondAsLine = !DrawBondAsLine;
    if ( DrawBondAsLine )
      DrawBackboneAsCurve = true;
    else
      DrawBackboneAsCurve = false;
    display();
    break;
  case M1_SPHERE_NET:
    AAptrs[0]->InParm.showSphereNet= !(AAptrs[0]->InParm.showSphereNet);
    display();
    break;
  case M1_SAVE:       // saveImage
    strcpy(cmd, "s");
    break;
  case M1_MOVIE_MODE:       //movie
    strcpy(cmd, "m");
    break;
  case M1_PICK_MODE:       //pick
    strcpy(cmd, "p");
    break;
  case M1_ROCK:       //pick
    strcpy(cmd, "r");
    break;
  case M1_RIBBON:                      // ribbon menu entries
    strcpy(cmd, "zribbon2z");
    break;
  case M1_STRUCT2ND:                   // structure menu entries
    strcpy(cmd, "zstructurez");
    break;
  case M1_RIB_WIDEN:                   // share width parameter b/w ribbon & 2nd_struct
    strcpy(cmd, "zribbonwidenz");
    break;
  case M1_RIB_NARROW:
    strcpy(cmd, "zribbonnarrowz");
    break;
  case M1_OPEN_INSTRUCT:
    openfileType="";
    FileBrowser();
    break;
  case M1_OPEN_PDB:
    openfileType="-pdb";
    FileBrowser();
    break;
  case M1_STEREO_OFF:
    strcpy(cmd, "zstereo_offz");
    break;
  case M1_ANAGLYPH_STEREO:
    strcpy(cmd, "zstereo_anaglyphz");
    break;
  case M1_COLOR_ANAGLYPH_STEREO:
    strcpy(cmd, "zstereo_color_anaglyphz");
    break;
  case M1_HARDWARE_STEREO:
    strcpy(cmd, "zstereo_hardwarez");
    break;
  case M1_CLIP_RESET:
    ZClip(ResetClip,ResetClip);
    break;
  case M1_CLIP_WIN:
    camera.clipSpace=camera.ScreenSpace;
    glutPostRedisplay();
    break;
  case M1_CLIP_OBJ:
    camera.clipSpace=camera.ObjectSpace;
    glutPostRedisplay();
    break;
  case M1_EXIT:
    exit(0);
    break;
  }
  for ( i=0; i<strlen(cmd); i++ )
    keyboard(cmd[i], 0, 0);

  if (value== M1_PICK_MODE)       //pick
  {
    displayMenubar();    
    glutSetWindow(windowId);   // refresh menubar
    glutPostRedisplay();
  }
}

// for surface step distance
void MenuBar::Menu1_2(int value)        // menu function part2 for MenuButton#1 : surface cut jump distance
{
  char cmd[50]="";
  sprintf(cmd, "zsurf_step %dz", value);
  for (unsigned int i=0; i<strlen(cmd); i++ )
    keyboard(cmd[i], 0, 0);
}

//
//extern int Overlap2;

// existing #item=NumImages-1, items in them menu is  NumImages-1+4
#define M2Item_Overlap  (NumImages+3)    
#define M2Item_PrintCrd (NumImages+6)

void MenuBar::Menu2(int value)
{
  if ( value>=M2_OVERLAP_A && value<=M2_OVERLAP_F)
  {
    // must have ( NumImages > 1 ) now 
    char ovtype='A'+ value-M2_OVERLAP_A;       // 'A', 'B', 'C'   // M2_OVERLAP_A|B|C are continouse
    showMolImg(NumImages+1, ovtype);
    glutSetMenu(button2_menuid);   // add "print overlapped crd" item to menu 
    if ( glutGet(GLUT_MENU_NUM_ITEMS) < M2Item_PrintCrd )   /* the "Print.." item does not exit */
      glutAddMenuEntry("Print Overlapped CRD", M2_PRINTCRD);
  }
  else 
  {
    if (!NO_OVERLAP)         // previous overlapped? 
    {
      // Overlap2 is the seconde overlapped structure, recover the coordinate first 
      //AAptrs[Overlap2]->currentmol=0;
      //AAptrs[Overlap2]->ReProcessData();
      glutSetMenu(button2_menuid);           // remove "print overlapped crd" item to menu 
      glutRemoveMenuItem(M2Item_PrintCrd);
    }

    switch (value)
    {
      case M2_NEXTSTRUCT:        //Next Structure
        keyboard('a', 0, 0);
        break;
      case M2_NEXTCHAIN:        //Next Chain(PDB only)
        keyboard('c', 0, 0);
        break;
      case M2_PRINTCRD:          //PrintOverlappedCrd();
        if (!NO_OVERLAP)
          AAptrs[Overlap2]->WriteOrgCrd(MsgBrd.GetText());
        value=M2_LAST+Overlap2+1;   // go through to display the image of printed crd
      default:
        showMolImg(value-M2_LAST, ' ');          
        break;
    }
  }
}

//
void MenuBar::subButtonCmd(int button, int state, int x, int y)
{
  int myWindow = glutGetWindow();
  for ( int i=0; i<nButtons; i++ )
  if ( myWindow==buttons[i].windowId )
  {
     if ( buttons[i].group == MB_PICK )
     {
       switch (buttons[i].nameId)
       {
         case MBP_QUIT: 
           Menu1_1(M1_PICK_MODE);
           break;
         case MBP_DISTANCE: 
           keyboard('d', 0, 0);
           break;
         case MBP_ANGLE: 
           keyboard('a', 0, 0);
           break;
         case MBP_TORANGLE: 
           keyboard('t', 0, 0);
           break;
       }
     }
  }
}


void  MenuBar::displayMenubar()
{
  glutSetWindow(windowId);
  glClearColor(0.8,0.8,0.8,1);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glutSwapBuffers(); 
  for ( int i=0; i<nButtons; i++ )
    buttons[i].display(); 
  glutSwapBuffers(); 
}

//
// nameId==1, First OptionButton 
// subGroup == 0, not a subgroup
// assume:  avg each letter occupies 6 pix
void MenuBar::addButton(char *buttonText, int subGroup, int nameId)
{
  int x, w, mywId;

  if ( windowId ==-1 )
    return;

  if ( *buttonText == '|' )
    w = 2;
  else
    w=6+strlen(buttonText)*6;
  if (subGroup==MB_NOTSUB)
  {
    x= X;
    X +=w;      // X == next button left edge
  }
  else 
  {
    if ( nameId==1)
      optionX=X;
    x=optionX;
    optionX += w;
  }
  mywId = glutCreateSubWindow(windowId, x, Y, w, MenuButtonHeight);
  glutDisplayFunc(displayNop);   // display
  // binding mouse click function to rotation command
  if ( *(buttonText+1) == 'X') 
    glutMouseFunc (RotateX);
  else if ( *(buttonText+1) == 'Y') 
    glutMouseFunc (RotateY);
  else if ( *(buttonText+1) == 'Z') 
    glutMouseFunc (RotateZ );
  else if (strstr(buttonText, "Front")!=NULL)
    glutMouseFunc (ClipFront);
  else if (strstr(buttonText, "Back")!=NULL)
    glutMouseFunc (ClipBack);
  else if (strstr(buttonText, "NoClip")!=NULL)
    glutMouseFunc(ClipNoClip);
  buttons[nButtons].set(mywId, buttonText, subGroup, nameId, windowId);
  nButtons++;

  menuinit();
} //addButton()

void MenuBar::menuinit(void)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();             
  glPushMatrix();
}  


//*********************************************************
//*    File Browser Menu
//*********************************************************

void MenuBar::Menu_FileBrowser (int value)    // menu function for FileBrowser
{
  if ( value <=1 )
  { 
    CHDIR("..");              // goto uplever
    FileBrowser();
  } 
  else
  {
    char dirbuf[LINE_MAXLEN];  
    struct STAT statbuf;

    FILE  *mydir=POPEN(DIR1CMD, "r");   // add list of current directory to menu
    for (int myindex=2; myindex < value; myindex++)
      fscanf(mydir, "%*s");
    fscanf(mydir, "%s", dirbuf);
    PCLOSE(mydir);

    if ( STAT (dirbuf, &statbuf) == 0 )
    {
      if ( statbuf.st_mode & DIRMASK )
      {
        CHDIR(dirbuf);             // goto uplever
        FileBrowser();
      } else {
        if (filelist_menuid>= 0 )       // destroy file browser when reading a file
        {
          glutDestroyMenu(filelist_menuid);    // destroy existing menu
          filelist_menuid=-1;
        }
        ReadFileParams(dirbuf,openfileType, NULL);
        AddStructToMenu();
        AddSphereToList();
        Overlap1=Overlap2;
        Overlap2=ImageIndex;
        display();
      }
    }
  }
} // Menu_FileBrowser()

// delete existing file sub-submenu and create a new file sub-submenu from current directory
// and attach it to the right button.
// return the menu id.
int MenuBar::submenu_filelist()
{
  int taglen1 = 3, taglen2=5;

  if (filelist_menuid>= 0 )
     glutDestroyMenu(filelist_menuid);    // destroy existing menu 

  filelist_menuid = glutCreateMenu(Menu_FileBrowser);      // create the menu
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  //browser
  char dirbuf[LINE_MAXLEN];  
  strncpy(dirbuf, "<--", taglen1);         // mark is a parent directory "<--"
  if ( GETCWD(dirbuf+taglen1, LINE_MAXLEN-taglen1) != NULL )
  { 
    dirbuf[LINE_MAXLEN-1]='\0';
    glutAddMenuEntry(dirbuf, 1);           // post new entries  == 1
  }

  struct STAT statbuf;

  FILE  *mydir=POPEN(DIR1CMD, "r");   // add list of current directory to menu

  int    nitems=1;
  while ( fscanf(mydir, "%s", dirbuf+taglen2) > 0 )
  {
    if ( STAT(dirbuf+taglen2, &statbuf) == 0 )
    {
      if ( statbuf.st_mode & DIRMASK )
        strncpy(dirbuf, "  -->", taglen2);     // is sub directory         "  -->"
      else 
        strncpy(dirbuf, "     ", taglen2);     // is regular file          "    "
      glutAddMenuEntry(dirbuf, ++nitems);
    }
  }
  PCLOSE(mydir);
  return filelist_menuid;

} // submenu_filelist()


void spoofClick( int x, int y )
{
#ifndef WIN32
	static int offset = 0;
	x += offset;
	y += offset;
	offset += 2;
	printf( "Spoofing mouse click at %d,%d\n", x,y);
  // XQueryTree() doesn't work well with Exceed 
  Display *display=XOpenDisplay(XDisplayName(NULL)); 
  Window  screen=DefaultScreen(display);
  Window  root=RootWindow(display, screen);
  Window  menu=FindWindow(display, root, CmoilAppWinName);
  
  XButtonEvent xbe;
  xbe.serial = 0; 
  xbe.send_event = True; 
  xbe.display = display; 
  xbe.window = menu; 
  xbe.root = root; 
  xbe.subwindow = menu; 
  xbe.time = 0; 
  xbe.x = x ; 
  xbe.y = y ; 
  xbe.x_root = 0; 
  xbe.y_root = 0; 
  xbe.state = Button1Mask; 
  xbe.button = Button1; 
  xbe.same_screen = True;
  
  xbe.type = ButtonPress;
  int status = XSendEvent(display, menu, false, ButtonPressMask, ( XEvent *) &xbe ); 
  printf( "XSendEvent (press) status = %d\n", status);

  xbe.type = ButtonRelease;
  xbe.time = 20;
  status = XSendEvent(display, menu, false, ButtonReleaseMask, ( XEvent *) &xbe ); 
  printf( "XSendEvent (release) status = %d\n", status);
  #endif
}




// #include "X11/extensions/XTest.h"
// uncomment this for the XTest code below
void MenuBar::ClickRightButton()
{
  //it should be in the leftmost menubutton window already
  //simulate right mouse button up and down to open filelist

#ifdef  _WIN32
  HWND frtW=GetWindow(GetWindow(GetActiveWindow(), GW_CHILD), GW_CHILD);  // get HWND of the "Display"
  PostMessage(frtW, WM_RBUTTONDOWN, 0, 0);
  PostMessage(frtW, WM_RBUTTONUP, 0, 0);
#else


	// @TODO (tfb): this strange mechanism of spoofing a right-mouse click works under win32 but
	// not under linux.  FIXME.
	
	
  // XQueryTree() doesn't work well with Exceed 
  Display *display=XOpenDisplay(XDisplayName(NULL)); 
  Window  screen=DefaultScreen(display);
  Window  root=RootWindow(display, screen);
  
  /*
  // I had thought to use these test extensions for X11 to simplify the spoof; but this
  // doesn't work either. (Tfb)
  
  int event_base, error_base, major_version, minor_version;
  XTestQueryExtension(display, &event_base, &error_base, &major_version, &minor_version);
  printf( "XTest results: event_base[%d] error_base[%d] major_v[%d] minor_v[%d]\n",
  				event_base,error_base,major_version,minor_version); 
  
  int button = 1;
  int is_press = 1;
  int delay = 0;
  XTestFakeButtonEvent(display, button, is_press, delay);
  is_press = 0;
  XTestFakeButtonEvent(display, button, is_press, delay);
  
  return;
  */

Window menu = FindWindow(display, root, CmoilAppWinName);
  
  printf( "[ClickRightButton]: found cmoil window: %d\n", menu );

  XButtonEvent xbe;
  xbe.serial = 0; 
  xbe.send_event = True; 
  xbe.display = display; 
  xbe.window = menu; 
  xbe.root = root; 
  xbe.subwindow = menu; 
  xbe.time = 0; 
  xbe.x = 5 ; 
  xbe.y = 5 ; 
  xbe.x_root = 10; 
  xbe.y_root = 39; 
  xbe.state = Button3Mask; 
  xbe.button = Button3; 
  xbe.same_screen = True;
  
  xbe.type = ButtonPress;
  int status = XSendEvent(display, menu, false, ButtonPressMask, ( XEvent *) &xbe ); 
  printf( "XSendEvent (press) status = %d\n", status);

  xbe.type = ButtonRelease;
  xbe.time = 20;
  status = XSendEvent(display, menu, false, ButtonReleaseMask, ( XEvent *) &xbe ); 
  printf( "XSendEvent (release) status = %d\n", status);

#endif

}

void MenuBar::FileBrowser()
{
  glutSetWindow(buttons[0].windowId);     // set active window to "Display" button

  submenu_filelist();         // delete old and create new dir sub-submenu, and attach it to the right button
    
  //simulate right mouse button up and down to open filelist
  ClickRightButton();
}

int MenuBar::submenu_surf () 
{
  int menu1_cutID = glutCreateMenu(Menu1_2);       // surf cut distance sub-submenu of surf submenu
  char itext[3];
  glutAddMenuEntry("0.5", -1);
  for (int i=1; i<=64; i+=i )
  {
    sprintf(itext, "%d", i);
    glutAddMenuEntry(itext,  i);
  } 

  int menuid  = glutCreateMenu(Menu1_1);      // surf submenu
  glutAddMenuEntry("On/Off",  M1_SURF);
  glutAddMenuEntry("----------------",  M1_SEPARATOR );
  glutAddMenuEntry("Transparent (On/Off) ",  M1_SURF_TRANS);
  glutAddMenuEntry("Mesh (On/Off) ",  M1_SURF_MESH);
  glutAddMenuEntry("Show selected area (On/Off)",  M1_SURF_AREA);
  glutAddMenuEntry("Show cavity only (On/Off)",  M1_SURF_CAVITY_OFF);
  glutAddMenuEntry("----------------",  M1_SEPARATOR );
  glutAddMenuEntry("No Cutoff",  M1_SURF_NOCUT);
  glutAddMenuEntry("Set surface cutoff to picked plane",  M1_SURF_PICK);
  glutAddMenuEntry("Switch surface cutoff axe Z->X->Y->",  M1_SURF_SWITCH_AXE);
  glutAddMenuEntry(" -> Cut flip over",  M1_SURF_FLIP);
  glutAddMenuEntry(" -> Cut farther", M1_SURF_FAR);
  glutAddMenuEntry(" -> Cut nearer",  M1_SURF_NEAR);
  glutAddMenuEntry("----------------",  M1_SEPARATOR );
  glutAddMenuEntry("Color Residues: RED(chg+), BLU(chn-), GRN(pol), GREY(hyd), PUR(unknow)",  M1_SURF_COLOR);
  glutAddMenuEntry("Color Atoms: Red(chg+) <--> Grey(0) <--> Blue(chg-)",  M1_SURF_COLOR_ATOM);
  glutAddMenuEntry("----------------",  M1_SEPARATOR );
  glutAddSubMenu("Cutoff Jump distance",  menu1_cutID);
  return menuid;
}

int MenuBar::submenu_ribbon() {
  int menuid = glutCreateMenu(Menu1_1);   // ribbon submenu
  glutAddMenuEntry("On/Off", M1_RIBBON);
  glutAddMenuEntry("Width Increase",  M1_RIB_WIDEN);
  glutAddMenuEntry("Width Decrease",  M1_RIB_NARROW);
  return menuid;
}

int MenuBar::submenu_stick()
{
  int menuid = glutCreateMenu(Menu1_1);   // 2nd struct submenu
  glutAddMenuEntry("On/Off", M1_STICK);
  glutAddMenuEntry("Width Increase",  M1_STICK_WIDEN);
  glutAddMenuEntry("Width Decrease",  M1_STICK_NARROW);
  glutAddMenuEntry("Quality high/low",  M1_STICK_QUALITY);
  return menuid;
}

int MenuBar::submenu_backbone()
{
  int menuid = glutCreateMenu(Menu1_1);   // 2nd struct submenu
  glutAddMenuEntry("On/Off", M1_BACKBONE);
  glutAddMenuEntry("Quality high/low",  M1_BACKBONE_QUALITY);
  return menuid;
}

int MenuBar::submenu_struct()
{
  int menuid = glutCreateMenu(Menu1_1);   // 2nd struct submenu
  glutAddMenuEntry("On/Off", M1_STRUCT2ND);
  glutAddMenuEntry("Width Increase",  M1_RIB_WIDEN);
  glutAddMenuEntry("Width Decrease",  M1_RIB_NARROW);
  return menuid;
}

int MenuBar::submenu_stereo()
{
  int menuid=glutCreateMenu(Menu1_1);    // stereo submenu
  glutAddMenuEntry("off",  M1_STEREO_OFF);
  glutAddMenuEntry("Anaglyph stereo (Red/Cyan)",  M1_ANAGLYPH_STEREO);
  glutAddMenuEntry("Color Anaglyph stereo (Red/Cyan)",  M1_COLOR_ANAGLYPH_STEREO);
  if (HasHardwareStereoSupport) 
    glutAddMenuEntry("Hardware stereo",  M1_HARDWARE_STEREO);
  return menuid;
}

int MenuBar::submenu_shine()
{
  int menuid=glutCreateMenu(Menu1_1);    // stereo submenu
  glutAddMenuEntry("no",  M1_SHINE_NO);
  glutAddMenuEntry("low", M1_SHINE_LOW);
  glutAddMenuEntry("medium", M1_SHINE_MEDIUM);
  glutAddMenuEntry("high",M1_SHINE_HIGH);
  return menuid;
}

int MenuBar::submenu_openfile()
{
  int menuid = glutCreateMenu(Menu1_1);    // file submenu
  glutAddMenuEntry("Open instruction file", M1_OPEN_INSTRUCT);
  glutAddMenuEntry("Open PDB file",         M1_OPEN_PDB);
  return menuid;
}

void MenuBar::menubutton1(void)
{
  // create menu for first menubar button 
  int menu1_surfID  =submenu_surf();      // surf submenu
  int menu1_ribbonID=submenu_ribbon();
  int menu1_structID=submenu_struct();
  int menu1_openID  =submenu_openfile();
  int menu1_stereoID=submenu_stereo();
  int menu1_shineID =submenu_shine();
  int menu1_stickID =submenu_stick();
  int menu1_backboneID=submenu_backbone();

  addButton("Display", MB_NOTSUB, MB_NOTSUB);       // menu botton #1
  glutCreateMenu(Menu1_1);
  glutAttachMenu(GLUT_LEFT_BUTTON);
    glutAddMenuEntry("Set Movie Mode               [m]",M1_MOVIE_MODE);
    glutAddMenuEntry("Set Pick Mode                  [p]", M1_PICK_MODE);
    glutAddMenuEntry("Set Rocking Mode              [r]", M1_ROCK);

    glutAddMenuEntry("----------------",  M1_SEPARATOR );
    glutAddSubMenu("Stereo Switch",     menu1_stereoID);
    glutAddMenuEntry("Save Image(s) to File(s)     [s]",  M1_SAVE);

    glutAddMenuEntry("----------------",  M1_SEPARATOR );
    glutAddSubMenu("Show Bonds As Sticks",menu1_stickID );
    glutAddSubMenu("Show Backbon",        menu1_backboneID);
    glutAddSubMenu("Show Ribbon",         menu1_ribbonID );
    glutAddSubMenu("Show 2nd Structure",  menu1_structID );
    glutAddSubMenu("Show Surface",        menu1_surfID);    

    glutAddMenuEntry("----------------",    M1_SEPARATOR );
    glutAddMenuEntry("Set White Backgroud", M1_WHITEBG);
    glutAddMenuEntry("Show Net for PickSphere (on|off)", M1_SPHERE_NET);
    glutAddSubMenu("Enhanced Display", menu1_shineID); 
    glutAddMenuEntry("High/Low Quality Shapes", M1_QUALITY);

    glutAddMenuEntry("----------------",  M1_SEPARATOR );
    glutAddSubMenu("Open File",           menu1_openID);    
    glutAddMenuEntry("Exit", M1_EXIT);
}

void MenuBar::menubutton2(void)
{
  // create menu for second menubar button 
  addButton( "Structure", MB_NOTSUB, MB_NOTSUB);  // menu botton #2
  button2_pulldownOverlap= glutCreateMenu(Menu2);       // pull-down submenu of overlap 
  glutAddMenuEntry("Auto-align with CE",                   M2_OVERLAP_A);       // default is NO_OVERLAP
  glutAddMenuEntry("Auto-align by backbone(N,CA,C,O)",     M2_OVERLAP_B);       // default is NO_OVERLAP
  glutAddMenuEntry("Auto-align by CA",                     M2_OVERLAP_C);       // default is NO_OVERLAP
  glutAddMenuEntry("----------------", M1_SEPARATOR );
  glutAddMenuEntry("Align by all selected atoms",          M2_OVERLAP_D);       // default is NO_OVERLAP
  glutAddMenuEntry("Align by backbone of selected atoms",  M2_OVERLAP_E);       // default is NO_OVERLAP
  glutAddMenuEntry("Align by CA of selected atoms",        M2_OVERLAP_F);       // default is NO_OVERLAP


  button2_menuid=glutCreateMenu(Menu2);
  glutAttachMenu(GLUT_LEFT_BUTTON);
  if ( NumImages > 0 )      // init menu2,  set init entry of images
  {
    char tmpstr[25];  
    int  currCrd;

    glutAddMenuEntry("Next Structure              [a]", M2_NEXTSTRUCT);    // #1
    glutAddMenuEntry("Next Chain(PDB only)   [c]", M2_NEXTCHAIN);     // #2
    glutAddMenuEntry("---------",                M2_SEPARATOR );    // #3

    for (currCrd = 1; currCrd <= NumImages; currCrd++)
    {
      sprintf(tmpstr, "MOLC #%d : %4s", currCrd, AAptrs[currCrd-1]->moleculeName);
      glutAddMenuEntry(tmpstr, currCrd+M2_LAST);   // the image index 
    }
    if (NumImages > 1)
    {
      glutAddMenuEntry("---------",                M2_SEPARATOR );    // #3
      glutAddSubMenu("Overlap Two Structures",  button2_pulldownOverlap);
    }
  }
}

void MenuBar::mainMenu(int pWindow)
{
  glFinish();
  // create menubar first
  parentId = pWindow;
  windowId=glutCreateSubWindow(pWindow, 0, 0, 2000, MenuButtonHeight+1);   // use reshape callback
  glutDisplayFunc(displayMenubar);

  // left most buttons
  menubutton1();
  menubutton2();

  // othere buttons
  //
  addButton("Rot", MB_NOTSUB, 0);
  addButton(" X", MB_NOTSUB, 0);
  addButton(" Y", MB_NOTSUB, 0);
  addButton(" Z", MB_NOTSUB, 0);
  addButton("Clip", MB_NOTSUB, 0);
    glutCreateMenu(Menu1_1);
    glutAttachMenu(GLUT_LEFT_BUTTON);
    glutAddMenuEntry("Reset",M1_CLIP_RESET);
  addButton("Front", MB_NOTSUB, 0);
  addButton("Back ", MB_NOTSUB, 0);
  addButton("NoClip", MB_NOTSUB, 0);

  addButton("Distance[d]",        MB_PICK, MBP_DISTANCE);
  glutMouseFunc(subButtonCmd);

  addButton("Angle[a]",           MB_PICK, MBP_ANGLE);
  glutMouseFunc(subButtonCmd);

  addButton("Torsional Angle[t]", MB_PICK, MBP_TORANGLE);
  glutMouseFunc(subButtonCmd);
  glFinish();
}

void MenuBar::AddStructToMenu()
{
  char tmpstr[25];  
  glutSetMenu(button2_menuid);

  if (NumImages == 1 ) {           // add Next item when adding structure is first
    glutAddMenuEntry("Next Structure  [a]",      M2_NEXTSTRUCT);    // #1
    glutAddMenuEntry("Next Chain(PDB only) [c]", M2_NEXTCHAIN);     // #2
    glutAddMenuEntry("---------",                M2_SEPARATOR );    // #3
  } 
  else if ( NumImages > 2 )      // existing overlap #item=NumImages-1
  {  
    while ( glutGet(GLUT_MENU_NUM_ITEMS) >= M2Item_Overlap)  // check existing menu items
    {
      glutRemoveMenuItem(M2Item_Overlap);   //  remove menu item after "Overlap...."
    }
  }

  sprintf(tmpstr, "MOLC #%d : %4s", NumImages, AAptrs[NumImages-1]->moleculeName);
  glutAddMenuEntry(tmpstr, NumImages+M2_LAST);
  if (NumImages > 1) 
  {
    glutAddMenuEntry("---------",             M2_SEPARATOR );    // #3
    glutAddSubMenu("Overlap Two Structures",  button2_pulldownOverlap);
    if (!NO_OVERLAP)
    {
      glutAddMenuEntry("Print Overlapped CRD", M2_PRINTCRD);
    }
  }
}

}
