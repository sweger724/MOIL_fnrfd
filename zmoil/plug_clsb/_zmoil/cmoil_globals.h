#ifndef _CMOIL_GLOBALS_H
#define  _CMOIL_GLOBALS_H

#include "cmoil_const.h"

namespace CMOIL {

class AtomArray;
class PickedAtoms;
class MessageBoard;
class Camera;
class SplineShape;
class Tube;
class Curve;

class State 
{
public:
  static AtomArray *AAptrs[ARRAY_MAXNUM];	// All input Structures
  static PickedAtoms PickQ;	 				// most recently selected 4 atoms
  static MessageBoard  MsgBrd;
  static MessageBoard  MsgBrd2;				// added for additional displays, lower right-hand corner
  static Camera  camera;

  static char *ExeDir;						// point to argv[0] 

  static int ImageIndex;					// index of current displayed/processed AtomArray, [0, NumImages), & NumImages(overlap)
  static int NumImages;						// total number of images for display
											// overlap display condition:	ImageIndex >= NumImages

  static int Overlap1;	// remember 2 recent access structure indexes for overlap purpose
  static int Overlap2;

  static int PickMode;						// global Pick Mode (0/1)
  enum PickInfoType { piNone, piAngle, piDistance, piTorsion };
  static PickInfoType PickModeInfo;
	// tracks the currently selected type of information for display.
	// this is used to auto-refresh the displayed values when viewing new struct
  static int StereoMode;
  static bool ColorAnaglyph;
  static GLfloat Shine;						// shine effect

  static bool HasHardwareStereoSupport;  
  static bool DrawBondAsLine;				// draw atom bond as line/stick
  static bool DrawBackboneAsCurve;			// draw atom bond as line/stick

};

}

#endif
