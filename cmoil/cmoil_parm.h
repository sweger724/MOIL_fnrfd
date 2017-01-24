#ifndef _CMOIL_PARM_H
#define _CMOIL_PARM_H

//*************************************************************************************************
//*  Filename:   cmoil_parm.h
//*
//*  Description: head file for InputParameters class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. ?? 2000	??         	Initial Development
//*  Oct. 09 2000	Baohua Wang	Extract from cmoil.cpp and modified to support crd, dcd coordinates
//*
//*************************************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmoil_const.h"
#include "cmoil_aminoacids.h"
#include "cmoil_globals.h"

namespace CMOIL {

const float MOIL_MAXCUTOFFVAL = 1000000.0;            // for surface cutoff

typedef enum {fNONE='\0', fPTH='p', fCRD='c', fDCD='d', fXYZ='x', fPDB='b' } tFileTag;

#define IsBinaryCrd(crdFileTag)    ( crdFileTag==fPTH || crdFileTag==fDCD)

class InputParameters : public State
{
public:
  InputParameters(void);
  ~InputParameters(void);

  char fileName[FILENAME_MAXLEN];
  int ReadInput(char *filename, int crdIndex=0, char *cmoil=NULL);	// return current CRD index
  char connName[FILENAME_MAXLEN];
  char filetag;							// p_pth, c_crd, d_dcd file, x_xyz, b_pdb file
  union {
     char pthName[FILENAME_MAXLEN];
     char crdName[FILENAME_MAXLEN];
     char dcdName[FILENAME_MAXLEN];
     char pdbName[FILENAME_MAXLEN];
  };
  char bmpName[FILENAME_MAXLEN];
  char surfName[FILENAME_MAXLEN];
  int dispmode;               // DisplayMode
  int bondmode;
  int swapmode;
  void SwapBite(void *booger, int varsize=8);
  NameString *skipRes;		// indexed by numRes [MONO_MAXNUM];
  NameString *stickRes;	      // [MONO_MAXNUM];
  NameString *spaceRes;		// [MONO_MAXNUM];
  NameString *ribbonRes;	// [MONO_MAXNUM];
  int structno;
  int structend;              // structure index from 0 to stop at, default=-1, no stop
  int skipresno;
  int stickresno;
  int spaceresno;
  int ribbonresno;
                             // text strings with pick syntax, lengthes are privates
  char *skipsurf;            // single line .. 
  char *picksphere;
  float pickradius;          // radius of picksphere
  float rotateDegree;        // rotate degree for each button click

  char *subsurf;
  char *skip;
  char *stick;
  char *space;
  char *ribbon;
  char *alignsettings;       // single line .. end
  char *colorsettings;       // multi-lines color pick list, separated by ", "
  char *ribboncolor;         // multi-lines color pick list for coloring ribbon

  int   ribbonpolymode;
  tRGBA ribbonColorBase;  // basic color for the ribbon, backbone etc, default green {0,255,0}
  tRGBA bg;               // background color   0-1
  float ribbonWidthFactor;   //  default=1.0;  >1.0 get wider;  <1.0 narrow
  float stickWidthFactor;    //  default=1.0;  >1.0 get wider;  <1.0 narrow

  float cutoffval;        //  surface display cutoff: only display cutoffaxe<=cutoffval, default z<=CoMZ
  char  cutoffaxe;        // ='X','Y','Z'
  float cutoffstep;       // far-near distance moved in one step, default=2
  float surfprobe;        // input protein surface probe redius
  bool  useExistingVtxFile;    // using an existing vertexes file from previous run if any

  int AAcolor[NUM_AA_TYPE][3]; // amino acid colors, indexed by AA_TYPE, RGB
  int AAcolorGroups;      //  defined by tSurfColorType 

  bool  showSubSurf;        // hightlight selected sub-surf 
  bool  showCavity;         // display cavity only
  bool  showSurfInMesh;       // display surface in mesh instead of color filling

  bool  showSphereNet;      // diaplay net for pickSphere
  GLubyte  transparent;     // 0-255

  float sleepseconds;
  int   maxframes; 	  // save max 500 frames for movie in default
  char  imagetype[4];     // supported image type gif, mpg, bmp, png
  int   gifloops;	
  int   framedistance;	  // frame distance in #structures
  int   structpstart;     // the first structure to be save into file    

  char exedir[FILENAME_MAXLEN];         // cmoil, mpeg2enc directory

  tLEN alignsettingsLen;

private:
  tLEN colorsettingsLen;               // memory lengthes passed into AddText
  tLEN ribboncolorLen;   
  tLEN skipsurfLen;                    // single line ..
  tLEN picksphereLen;
  tLEN pickcavityLen;
  tLEN subsurfLen;
  tLEN skipLen;
  tLEN stickLen;
  tLEN spaceLen;
  tLEN ribbonLen;
             
  char *AddText(char **allocStr, char *newStr, unsigned int *lenVar, char *separator);
  void DropText(char **allocStr, unsigned int *lenVar);

};

void SetValidPath(char *pathname);

}

#endif
