#ifndef _CMOIL_BSPLINE_H
#define _CMOIL_BSPLINE_H

//*************************************************************************************************
//*  Filename:   cmoil_bspline.h
//*
//*  Description: head file for Bspline class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 02 2005	Baohua Wang	Reshape
//*
//*************************************************************************************************
#include "cmoil_atom.h"

namespace CMOIL { 

class Bspline : public State
{
protected:

private:
  int N;	              // number of interpolated points per segment
  int J;                  // number of points inputed
  int numIn;              // count of current points in a_in
  Atom psudoAtom;

  enum { MaxNumSegments=10, NP=4, /*NumInputPoints=4*/ }; 
  float M[NP][NP],S[NP][NP], B[NP][NP], G[NP][NP], BG[NP][NP];

  void StepFit(bool isLastInput=false);

public:
  Atom      *a_in[NP];       // circular array of input points
  Coordinate a_out[MaxNumSegments];     // a_out= pointer to (ouput) array
  Atom      *out_atom[MaxNumSegments];  // a_out corresponding atom 
  int        currIn;        // current index  of a_in;
  int        numOut;        // count of output points in a_out
  AtomArray *aaptr;         // the structure current working on

  Bspline(AtomArray *inArray=NULL, int N=BsplineInterPoints);
  void init(AtomArray *inArray, int N=BsplineInterPoints);
  void NextPoint(Atom *aptr);
  void reset() { currIn=-1; numIn=0; numOut=0; };

};

}

#include "cmoil_drawshape.h"

#endif
