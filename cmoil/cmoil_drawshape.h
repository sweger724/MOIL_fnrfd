#ifndef _CMOIL_DRAWSHAPE_H
#define _CMOIL_DRAWSHAPE_H

//*************************************************************************************************
//*  Filename:   cmoil_drawshape.h
//*
//*  Description: head file for display structural information 
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Nov. 02 2005	Baohua Wang	Reshape
//*
//*************************************************************************************************
#include "cmoil.h"
#include "cmoil_crd.h"
#include "cmoil_const.h"
#include "cmoil_bspline.h"

namespace CMOIL {

class Bspline;

class SplineShape  : public Bspline 
{
private:
  int numCA;        			// used by quad()  :
  Coordinate v[4];  			// for quad() return value
  Coordinate preCA, preO, preD;

protected:
  AtomArray  *aaPtr;  

  Coordinate* quad(Coordinate CA, Coordinate O, int reset); 
  int colorFilter;
  int windowSize;

  static const float StdRibbonWidth;    // use with amplify factor
  double ribbonWidth;   // actual width between bs1 & bs2, if ribbonWidth=0, bs1&bs2 merged

  const int numBss;              // number of bspline curves: 1, 2
  Bspline    bss[2];             // max dim: 1 for backbone, 2 for others

  Atom       input[2][4];        // input array for two spline
  int        currCaIndex;        // current CA index of input[0] for splines of ribbon&secStruct 

  Coordinate lastDraw[2];        // previous pair of points, for ribbon & 2nd struct numBss=2;
  Coordinate lastNorm1;
  Coordinate *nvp;

  tRGBA  *colorbase;                // base ribbon color
  tRGBA  *colorsetting;             // ribbon color input array
  //tRGBA  *cp;                       // current color

public:
  SplineShape () : numBss(1) { ; };
  SplineShape  (AtomArray *aaptr, int windowsize, int colorfilter, int numbss) ;
  
  virtual void ChainBegin() ;
  virtual void ChainEnd();

  virtual void DrawBegin(SecStructType type=STRUCT_UNKNOWN) {;}; 
  virtual void DrawEnd()  { ; };

  virtual void Draw(Atom *myatom) {;}  ;   // single point bslpline, eg. backbone    
  virtual void Draw(Coordinate newCA, Coordinate newO, int rn, SecStructType secStruct=STRUCT_UNKNOWN) {;} ;
                                           // peptide plan bspline,  eg. ribbon & 2nd struct 
  void setColor(int rn, tRGBA *dflColor); 
};

//=======================================================================================

// for backbones
class SplineBelt1D: public SplineShape
{
public:
   SplineBelt1D(AtomArray *aaptr,int windowSize,int colorfilter):SplineShape(aaptr,windowSize,colorfilter,1){ribbonWidth=0;} ;
};

// ribbons draw as 2D belts
class SplineBelt2D: public SplineShape
{
public:
   SplineBelt2D(AtomArray *aaptr,int windowSize,int colorfilter):SplineShape(aaptr,windowSize,colorfilter,2){;};
};

//ribbons draw as 3D belts
const int NumEdges=4;   // four faces for the 3D ribbon
class SplineBelt3D : public SplineShape
{
protected:
   static const float thickness;

   Coordinate lastDraw1, lastDraw2;
   double structFactor;

   SecStructType secStruct;
   SecStructType oldStructLow;
   SecStructType prevStruct;
   SecStructType bsplineOutStruct;  //prev to prevStruct, struct type corresponding to the bspline output segment

   // draw different shape upon bsplineOut[]
   void DrawQuadstrip(int startIndex, int stopIndex, int edge1, int edge2, bool SetDrawDelta, tRGBA *dflColor);
   void Draw4Quadstrips(int startIndex,int outIndex, tRGBA *dflColor);
   void DrawQuadstripEnd(int outIndex, tRGBA *dflColor);

public:
   SplineBelt3D(AtomArray *aaptr,int windowSize,int colorfilter);
   ~SplineBelt3D();

   void Draw(Coordinate newCA, Coordinate newO, int rn, SecStructType newStruct=STRUCT_UNKNOWN);
   //void DoSpline(Coordinate newCA, Coordinate newO, int rn, SecStructType newStruct=STRUCT_UNKNOWN);

   int         bsplineOutCount;
   Coordinate *bsplineOut[NumEdges];       // indexed by the output bsplineOutCount
   int        *resNum;                     // indexed by the output bsplineOutCount
   SecStructType *structTypes;

   void ChainBegin() { SplineShape::ChainBegin(); 
                       bsplineOutStruct=STRUCT_UNKNOWN; 
                       oldStructLow=STRUCT_UNKNOWN; 
                       bsplineOutCount=0; }; 
};

//=======================================================================================

class Curve : public SplineBelt1D
{
private:
  void Draw(Coordinate center) { glVertex3d(center.x, center.y, center.z); };

public:
  Curve(AtomArray *aaptr,int windowSize,int colorfilter):SplineBelt1D(aaptr,windowSize,colorfilter) {;} ;
  
  void Draw(Atom *myatom) ;
  void DrawBegin(SecStructType type=STRUCT_UNKNOWN) { glBegin(GL_LINE_STRIP); glNormal3d(0,0,1); };    
  void DrawEnd()    { Draw(NULL); glEnd(); }; 
     
  void ChainBegin() { SplineShape::ChainBegin(); DrawBegin(); };
  void ChainEnd()   { SplineShape::ChainEnd();   DrawEnd();   }; 
                
};


enum { TubeEndListNum=2, TubeDiskSlice=20 }; 
typedef float DiskArray[TubeDiskSlice][3];

class Tube : public SplineBelt1D
{
    
private:

  static const DiskArray UnitCircle;   // unit circle for speed

  int count;             // existing points
  DiskArray preDisk;     // translated disk
  Coordinate preCntr;
  Coordinate preDirect;

public:
  Tube(AtomArray *aaptr,int windowSize,int colorfilter):SplineBelt1D(aaptr,windowSize,colorfilter) {count=0;} ;

  void Draw(Atom *myatom) ;
  void Draw(Coordinate center, SecStructType type) ;
  void DrawBegin(SecStructType type=STRUCT_UNKNOWN)  { count=0; };   // reset for drawing next tube
  void DrawEnd()    { Draw(NULL); if (count>0) DrawEndBall(preCntr); };

  void ChainBegin() { SplineShape::ChainBegin(); DrawBegin(); };
  void ChainEnd()   { SplineShape::ChainEnd();  DrawEnd();  }; 

  static void AddToQuadricList(GLUquadricObj *qobj);

  static void  RotTransDisk(Coordinate center, Coordinate normal, DiskArray &disk);            
  static void  DrawEndBall(Coordinate &center);
  static const float TubeRadius;       // tube radius
  static void  DrawTubeBtw2Centers(Coordinate center, Coordinate &preCntr, DiskArray &preDisk, bool isFirst); 

};


class Ribbon : public SplineBelt2D  
{
 public:
   Ribbon (AtomArray *aaptr, int windowSize, int colorfilter):SplineBelt2D(aaptr,windowSize,colorfilter){;} ;
   void Draw(Coordinate center, SecStructType secStruct=STRUCT_UNKNOWN) {;} ; 
   void Draw (Coordinate newCA, Coordinate newO, int rn, SecStructType secStruct=STRUCT_UNKNOWN) ;   // draw ribbon from CA,O and rn
   void DrawDo () ; 

   void DrawBegin(SecStructType type=STRUCT_UNKNOWN) { glBegin(GL_QUAD_STRIP); }  ;    
   void DrawEnd() {DrawDo(); glEnd();  };
};


class Ribbon3D : public SplineBelt3D
{

 public:
   Ribbon3D(AtomArray *aaptr,int windowSize,int colorfilter):SplineBelt3D(aaptr,windowSize,colorfilter){;}

   void Draw(Coordinate center, SecStructType secStruct=STRUCT_UNKNOWN) {;} ; 
   void DrawDo () ;   // draw all points after spline

   void DrawBegin   ( SecStructType type=STRUCT_UNKNOWN) { ; };
   void DrawEnd()   { DrawDo(); };

  // void ChainBegin() { SplineShape::ChainBegin(); }; 
};

class SecStruct : public SplineBelt2D  
{
 private:
   SecStructType secStruct;
   Coordinate *nvp, lastDraw1, lastDraw2;
   double structFactor;
   SecStructType oldStructLow;
   SecStructType prevStruct;
   SecStructType bsplineOutStruct;  //prev to prevStruct, struct type corresponding to the bspline output segment

   void Draw ();
   void StructEndBegin(SecStructType oldStruct, SecStructType newStruct);

 public:
   SecStruct(AtomArray *aaptr, int windowSize, int colorfilter);

   void Draw(Coordinate center, SecStructType secStruct=STRUCT_UNKNOWN) {;} ; 
   void Draw(Coordinate newCA, Coordinate newO, int rn, SecStructType secStruct=STRUCT_UNKNOWN) ;

   void DrawBegin(SecStructType type=STRUCT_UNKNOWN) ;    
   void DrawEnd()    ;

   void ChainBegin();

};

class SecStruct3D : public SplineBelt3D
{
 private:
   void Draw() {};
   void StructEndBegin(SecStructType oldStruct, SecStructType newStruct);

   void DrawHelix(int startOutIndex, int stopOutIndex);
   void DrawSheet(int startIndex, int stopIndex);
   void DrawSheetArrow(int startIndex, int stopIndex);
   void DrawLoop(int startIndex, int stopIndex);

   void GetLoopCenters(int startIndex, int stopIndex);

 public:
   SecStruct3D(AtomArray *aaptr,int windowSize,int colorfilter);

   void DrawBegin(SecStructType type=STRUCT_UNKNOWN) ;    
   void DrawEnd()    ;

   void DrawDo () ;   // draw all points after spline

};



class Stick 
{
private :
  static int listNum(double atomDistance) { return int(atomDistance*100.+StickListStart); }; 

public:
  static void Draw(Atom *ptr1, Atom *ptr2);
  static void AddToQuadricList(Atom *ptr1, Atom *ptr2);

  static const int    StickListStart;
  static const double StickRadius;
};

}

#endif
