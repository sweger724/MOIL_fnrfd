#ifndef _CMOIL_PICK_H
#define _CMOIL_PICK_H
//*************************************************************************************************
//*  Filename:   cmoil_pick.hpp
//*
//*  Description: head file for cmoil.cpp
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 2000   	??         	Initial Development
//*  Oct. 09 2000	Baohua Wang	Extract from cmoil.cpp and modified
//*  Sep. 07 2007   Thomas Blom Update for OSX/Linux
//*************************************************************************************************
#ifdef __APPLE__
	#include <GLUT/glut.h>
#elif defined(WIN32)
	#include "GLwin32/glut.h"
#else
	#include <GL/glut.h>
#endif
#include "cmoil_atom.h"
#include "cmoil_msgboard.h"

namespace CMOIL {
class PickedAtoms : public State
{
private:
  int 	 q[ATOM_MAX_PICK+1]; 		// picked atom index 
  AtomArray  *aaptr[ATOM_MAX_PICK+1];	// pointed Images
  Atom	 *aptr[ATOM_MAX_PICK+1];      // picked atom array
  int		 numentries;
  
  int		 currImage;				// check imageIndex for reset purpose

  int        image1;
  int        image2;				// image2 = -1  if not overlap
  int 	 get_atomnum(int x, int y);	// get atom index from the input position

public:

  PickedAtoms();
  void reset(int imageIndex1=-1, int imageIndex2=-1);
  void put(int x, int y);		// put picked atom to the q array
  void print_distance();
  void print_angle();
  void print_torsional();
  void highlight(char rgbtype);	// color operation
  void darken(char rgbtype);		// color operation
  void reverse();				// color reverse

  void footer_print(int pickmode);
  void footer_erase();
  void display_pick(int windowSize=2, int colorFilter=ColorFilterNone, GLubyte ballR=255, GLubyte ballG=255, GLubyte ballB=0);    // == 1(small), 2(middle), 3(large)

  int  getPickedCrds(Coordinate *crds);

  char *atom_name( int atomIndex,  char fileTag=fPTH ) ; // not for PDB by default

  int isPicked(Atom *atom);
};

void display();

}

#endif
