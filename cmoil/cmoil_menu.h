//*************************************************************************************************
//*  Filename:   cmoil_menu.h
//*
//*  Description: 
//*    Main menu for cmoil.  
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Jun. 25 2003	Baohua Wang	Initial Development
//*  Aug. 02 2005	Baohua Wang	Reshape
//*
//*************************************************************************************************
// to manager Menus
#include "cmoil.h"
#include "cmoil_globals.h"
#include "cmoil_atom.h"
#include "cmoil_msgboard.h"

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

namespace CMOIL {

//button group ID
enum MopButtonGroup { MB_NOTSUB=0, MB_PICK, MB_MOVIE } ;

class MenuButton : public State
{
private:
  char *text;
  int  parentId;

public:
  int  windowId;
  int  nameId;      // for pick name
  int  group;

  void set(int mywindowId, char* mytext, int mygroup, int mynameId, int pWindow);
  void display(void);
  void draw(const char *dstr);
};

class MenuBar : public Main
{
public:
  static void mainMenu(int parentWindowId);   // create menubar 

  enum {MAX_NUM_MBUTTONS=20, MenuButtonHeight=20 } ;
private:

  static int  windowId;
  static int  parentId;    // parent window id
  static int  X;           // for button auto layout
  static int  optionX;     //
  static int  Y;           //
  static int  nButtons;
  static MenuButton buttons[MAX_NUM_MBUTTONS];

  // special menuids
  static int button2_menuid;
  static int button2_pulldownOverlap;
  static int filelist_menuid;

  static char *openfileType;

  static void displayMenubar();
  static void addButton(char *buttonText, int subGroup, int nameId) ;

  static void menuinit(void);

  static void menubutton1(void);   // leftmost fixed menu buttons
  static void menubutton2(void);

  static int submenu_openfile(void);   // return menu id
  static int submenu_stereo(void);
  static int submenu_struct(void);
  static int submenu_ribbon(void);
  static int submenu_surf(void);
  static int submenu_shine(void);  
  static int submenu_backbone(void);  
  static int submenu_stick(void);  

  static int submenu_filelist(void);  

  static void ClickRightButton (void);  
  static void FileBrowser(void);
  static void AddStructToMenu(void);


  static void displayNop();

  // Mouse callbacks   
  static void RotateX(int button, int state, int x, int y);
  static void RotateY(int button, int state, int x, int y);
  static void RotateZ(int button, int state, int x, int y);
  static void ClipFront (int button, int state, int x, int y);
  static void ClipBack  (int button, int state, int x, int y);
  static void ClipNoClip(int button, int state, int x, int y);
  static void subButtonCmd(int button, int state, int x, int y);

  // Menu callbacks
  static void Menu_FileBrowser(int value); 
  static void Menu1_1(int value); 
  static void Menu1_2(int value); 
  static void Menu2  (int value); 

  enum MBpickItems { MBP_QUIT=-1, MBP_DISTANCE=1, MBP_ANGLE, MBP_TORANGLE } ; // first item in subMenu must be 1

  enum M1items { M1_SEPARATOR=-999, M1_MOVIE_MODE=1, M1_PICK_MODE, M1_ROCK, M1_SAVE, 
               M1_STICK,
               M1_BACKBONE, M1_RIBBON, M1_STRUCT2ND, 
               M1_SURF, M1_SURF_COLOR, M1_SURF_COLOR_ATOM, M1_SURF_NOCUT, 
               M1_SURF_SWITCH_AXE, M1_SURF_FAR, M1_SURF_NEAR,
               M1_SURF_PICK, M1_SURF_FLIP, M1_SURF_AREA, M1_SURF_TRANS,
               M1_SURF_CAVITY_OFF, M1_SURF_MESH,
               M1_WHITEBG, M1_SHINE_NO, M1_SHINE_LOW, M1_SHINE_MEDIUM, M1_SHINE_HIGH,
               M1_BACKBONE_QUALITY, M1_STICK_QUALITY, M1_QUALITY,
               M1_SPHERE_NET,
               M1_RIB_WIDEN, M1_RIB_NARROW, M1_STICK_WIDEN, M1_STICK_NARROW,
               M1_OPEN_INSTRUCT, M1_OPEN_PDB, 
               M1_STEREO_OFF, M1_ANAGLYPH_STEREO, M1_COLOR_ANAGLYPH_STEREO, M1_HARDWARE_STEREO,
               M1_SPHERE_OFF, M1_SPHERE_RADIUS_INCR, M1_SPHERE_RADIUS_DECR,
               M1_CLIP_RESET, M1_CLIP_OBJ, M1_CLIP_WIN,
               M1_EXIT } ;
  enum M2items { M2_SEPARATOR=-999, M2_NEXTSTRUCT=1, M2_NEXTCHAIN=2,  
               M2_OVERLAP_A=51, M2_OVERLAP_B=52, M2_OVERLAP_C=53, M2_OVERLAP_D=54, 
               M2_OVERLAP_E=55, M2_OVERLAP_F=56,
               M2_PRINTCRD=90, M2_LAST=100 };
}; // MenuBar

}


